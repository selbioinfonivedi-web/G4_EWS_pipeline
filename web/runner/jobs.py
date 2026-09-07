"""Job registry: spawn a whitelisted command, stream its output, keep the record.

Processes are spawned with :func:`asyncio.create_subprocess_exec` on an
argv list -- never through a shell -- so nothing a caller supplies can be
interpreted as a shell metacharacter.

Only one job runs at a time. The pipeline stages write shared artifacts
(the Atlas, the QC FASTA, the append-only ledger); letting two run
concurrently would let them interleave writes to the same files, so
requests queue instead.
"""

from __future__ import annotations

import asyncio
import os
import signal
import time
import uuid
from collections import deque
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from . import audit
from .commands import GATE_CLOSED_EXIT, REPO_ROOT, Command, resolve_executable

MAX_LOG_LINES = 4000
MAX_JOBS_KEPT = 60


class JobState(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    PAUSED = "paused"
    SUCCEEDED = "succeeded"
    #: Exit code 3 from a gate-aware command. A correct outcome, not a failure.
    GATE_CLOSED = "gate_closed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class LogLine:
    seq: int
    stream: str  # "stdout" | "stderr" | "meta"
    text: str


@dataclass
class Job:
    id: str
    command_key: str
    title: str
    argv: list[str]
    pathogen: str | None
    created: float
    state: JobState = JobState.QUEUED
    started: float | None = None
    finished: float | None = None
    exit_code: int | None = None
    lines: deque[LogLine] = field(default_factory=lambda: deque(maxlen=MAX_LOG_LINES))
    _seq: int = 0
    paused_seconds: float = 0.0
    _paused_at: float | None = None
    _process: asyncio.subprocess.Process | None = None
    _subscribers: list[asyncio.Queue] = field(default_factory=list)

    # ── log plumbing ────────────────────────────────────────────────
    def emit(self, stream: str, text: str) -> None:
        self._seq += 1
        line = LogLine(seq=self._seq, stream=stream, text=text)
        self.lines.append(line)
        for queue in list(self._subscribers):
            try:
                queue.put_nowait(line)
            except asyncio.QueueFull:
                pass

    def subscribe(self) -> asyncio.Queue:
        queue: asyncio.Queue = asyncio.Queue(maxsize=1000)
        self._subscribers.append(queue)
        return queue

    def unsubscribe(self, queue: asyncio.Queue) -> None:
        if queue in self._subscribers:
            self._subscribers.remove(queue)

    @property
    def terminal(self) -> bool:
        return self.state in {
            JobState.SUCCEEDED,
            JobState.GATE_CLOSED,
            JobState.FAILED,
            JobState.CANCELLED,
        }

    @property
    def duration(self) -> float | None:
        """Wall time actually spent computing, excluding paused intervals."""
        if self.started is None:
            return None
        paused = self.paused_seconds
        if self._paused_at is not None:
            paused += time.time() - self._paused_at
        return (self.finished or time.time()) - self.started - paused

    def summary(self) -> dict:
        return {
            "id": self.id,
            "command": self.command_key,
            "title": self.title,
            "pathogen": self.pathogen,
            "argv": self.argv,
            "state": self.state.value,
            "exit_code": self.exit_code,
            "created": self.created,
            "started": self.started,
            "finished": self.finished,
            "duration": self.duration,
            "paused_seconds": round(self.paused_seconds, 2),
            "line_count": len(self.lines),
        }

    def detail(self) -> dict:
        return {
            **self.summary(),
            "lines": [{"seq": ln.seq, "stream": ln.stream, "text": ln.text} for ln in self.lines],
        }


class JobRunner:
    """Serialised runner for whitelisted pipeline commands."""

    def __init__(self, repo_root: Path = REPO_ROOT) -> None:
        self.repo_root = repo_root
        self.jobs: dict[str, Job] = {}
        self.order: deque[str] = deque()
        self._queue: asyncio.Queue[str] = asyncio.Queue()
        self._worker: asyncio.Task | None = None
        self._loop: asyncio.AbstractEventLoop | None = None
        self.current: str | None = None

    # ── lifecycle ───────────────────────────────────────────────────
    def start(self) -> None:
        """Start the worker. Must be called from the event loop (lifespan)."""
        self._loop = asyncio.get_running_loop()
        if self._worker is None or self._worker.done():
            self._worker = asyncio.create_task(self._run_forever())

    def _enqueue(self, job_id: str) -> None:
        """Hand a job to the worker from whichever thread we are on.

        FastAPI runs synchronous endpoints in a threadpool, where
        ``asyncio.Queue`` is not safe to touch directly: ``put_nowait``
        would append to the deque but resolve the waiting future from the
        wrong thread, and the worker would never wake. Route through the
        loop instead whenever we are not already on it.
        """
        loop = self._loop
        if loop is None:
            self._queue.put_nowait(job_id)
            return
        try:
            running = asyncio.get_running_loop()
        except RuntimeError:
            running = None
        if running is loop:
            self._queue.put_nowait(job_id)
        else:
            loop.call_soon_threadsafe(self._queue.put_nowait, job_id)

    async def stop(self) -> None:
        if self._worker:
            self._worker.cancel()
            try:
                await self._worker
            except asyncio.CancelledError:
                pass
            self._worker = None

    # ── submission ──────────────────────────────────────────────────
    def submit(self, command: Command, pathogen: str | None, argv: list[str]) -> Job:
        job = Job(
            id=uuid.uuid4().hex[:12],
            command_key=command.key,
            title=command.title,
            argv=argv,
            pathogen=pathogen,
            created=time.time(),
        )
        self.jobs[job.id] = job
        self.order.append(job.id)
        self._evict()
        job.emit("meta", f"$ {' '.join(argv)}")
        # Recorded at submission, not only at completion: a job that is
        # queued and then lost to a crash still happened, and an audit
        # that only logs finished work cannot show that.
        audit.record("queued", job)
        if self.current is not None:
            job.emit("meta", "Queued — another job is running. Stages share artifacts, so they run one at a time.")
        self._enqueue(job.id)
        return job

    def _evict(self) -> None:
        while len(self.order) > MAX_JOBS_KEPT:
            stale = self.order.popleft()
            job = self.jobs.get(stale)
            if job and not job.terminal:
                self.order.appendleft(stale)
                return
            self.jobs.pop(stale, None)

    def pause(self, job_id: str) -> bool:
        """Suspend the whole process group with SIGSTOP."""
        job = self.jobs.get(job_id)
        if job is None or job.state is not JobState.RUNNING or job._process is None:
            return False
        try:
            os.killpg(os.getpgid(job._process.pid), signal.SIGSTOP)
        except (ProcessLookupError, PermissionError):
            return False
        job.state = JobState.PAUSED
        job._paused_at = time.time()
        job.emit("meta", "Paused — SIGSTOP sent to the process group. Elapsed time excludes the pause.")
        return True

    def resume(self, job_id: str) -> bool:
        job = self.jobs.get(job_id)
        if job is None or job.state is not JobState.PAUSED or job._process is None:
            return False
        try:
            os.killpg(os.getpgid(job._process.pid), signal.SIGCONT)
        except (ProcessLookupError, PermissionError):
            return False
        if job._paused_at is not None:
            job.paused_seconds += time.time() - job._paused_at
            job._paused_at = None
        job.state = JobState.RUNNING
        job.emit("meta", "Resumed — SIGCONT sent to the process group.")
        return True

    def cancel(self, job_id: str, *, force: bool = False) -> bool:
        job = self.jobs.get(job_id)
        if job is None or job.terminal:
            return False
        if job.state is JobState.QUEUED:
            job.state = JobState.CANCELLED
            job.finished = time.time()
            job.emit("meta", "Cancelled before it started.")
            self._close(job)
            return True
        if job._process is not None and job._process.returncode is None:
            sig = signal.SIGKILL if force else signal.SIGTERM
            job.emit("meta", f"{'Force stopping' if force else 'Stopping'} — sending {sig.name}.")
            try:
                pgid = os.getpgid(job._process.pid)
                # A stopped process cannot act on SIGTERM. Continue it
                # first, or a paused job would sit there holding the queue.
                if job.state is JobState.PAUSED:
                    os.killpg(pgid, signal.SIGCONT)
                    if job._paused_at is not None:
                        job.paused_seconds += time.time() - job._paused_at
                        job._paused_at = None
                    job.state = JobState.RUNNING
                os.killpg(pgid, sig)
            except (ProcessLookupError, PermissionError):
                job._process.kill() if force else job._process.terminate()
            return True
        return False

    # ── worker ──────────────────────────────────────────────────────
    async def _run_forever(self) -> None:
        while True:
            job_id = await self._queue.get()
            job = self.jobs.get(job_id)
            if job is None or job.state is not JobState.QUEUED:
                continue
            self.current = job_id
            try:
                await self._execute(job)
            except Exception as exc:  # noqa: BLE001 - surfaced to the UI, never swallowed
                job.state = JobState.FAILED
                job.finished = time.time()
                job.emit("meta", f"Runner error: {exc}")
                audit.record("failed", job)
                self._close(job)
            finally:
                self.current = None

    async def _execute(self, job: Job) -> None:
        resolved = resolve_executable(job.argv[0])
        if resolved is None:
            job.state = JobState.FAILED
            job.finished = time.time()
            job.emit("meta", f"{job.argv[0]!r} was not found. See docs/installation.md.")
            audit.record("failed", job)
            self._close(job)
            return
        argv = [resolved, *job.argv[1:]]

        job.state = JobState.RUNNING
        job.started = time.time()
        audit.record("started", job)

        env = dict(os.environ)
        env["PYTHONUNBUFFERED"] = "1"
        env["COLUMNS"] = "100"

        process = await asyncio.create_subprocess_exec(
            *argv,
            cwd=str(self.repo_root),
            env=env,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            start_new_session=True,
        )
        job._process = process

        await asyncio.gather(
            self._pump(job, process.stdout, "stdout"),
            self._pump(job, process.stderr, "stderr"),
        )
        code = await process.wait()

        job.exit_code = code
        job.finished = time.time()
        command = _command_for(job.command_key)

        if job.state is JobState.CANCELLED:
            pass
        elif code == 0:
            job.state = JobState.SUCCEEDED
        elif code == GATE_CLOSED_EXIT and command is not None and command.gate_aware:
            job.state = JobState.GATE_CLOSED
            job.emit(
                "meta",
                "Exit 3 — the D.H1 gate is closed. This is a correct outcome, not a failure: "
                "scoring stays blocked until the gate returns SUPPORTED.",
            )
        else:
            job.state = JobState.FAILED

        job.emit("meta", f"Exited with code {code} after {job.duration:.1f}s.")
        # One row per terminal outcome. GATE_CLOSED is recorded as a
        # finished run, not a failure: exit 3 is a correct scientific
        # result and an audit that called it a failure would misreport it.
        audit.record("cancelled" if job.state is JobState.CANCELLED else "finished", job)
        self._close(job)

    async def _pump(self, job: Job, stream: asyncio.StreamReader | None, name: str) -> None:
        if stream is None:
            return
        while True:
            try:
                raw = await stream.readline()
            except ValueError:
                # A single line longer than the reader's limit; keep going.
                job.emit("meta", "(a very long output line was truncated)")
                continue
            if not raw:
                return
            job.emit(name, raw.decode("utf-8", errors="replace").rstrip("\n"))

    def _close(self, job: Job) -> None:
        for queue in list(job._subscribers):
            try:
                queue.put_nowait(None)
            except asyncio.QueueFull:
                pass

    # ── queries ─────────────────────────────────────────────────────
    def list_jobs(self) -> list[dict]:
        return [self.jobs[i].summary() for i in reversed(self.order) if i in self.jobs]

    def get(self, job_id: str) -> Job | None:
        return self.jobs.get(job_id)


def _command_for(key: str) -> Command | None:
    from .commands import BY_KEY

    return BY_KEY.get(key)
