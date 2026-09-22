"""Section 12 two-layer enforcement: the D.H1 gate as executable code.

Build Architecture Section 6 makes Stage 5 (scoring) contingent on a
passed D.H1 test, per pathogen, and Section 16 makes that a procedural
gate as well as an architectural one. This module is the mechanism. It
answers exactly one question — *may scoring run on real data for this
pathogen?* — from the persisted testing ledger, and it answers it
conservatively:

* the ledger is the authority, not a flag someone set;
* a pathogen with no ledger rows is **not** permitted (absence of a
  negative result is not a positive one);
* ``INSUFFICIENT_DATA`` is kept strictly distinct from ``NOT_SUPPORTED``
  — the first means the test never ran, the second means it ran and the
  hypothesis failed, and neither permits scoring;
* both layers must hold: ``operational_mode: true`` in the config **and**
  a ``SUPPORTED`` D.H1 verdict in the ledger. Flipping the config flag
  alone does nothing.

Stage 5 and Stage 6 entry points call :func:`assert_scoring_permitted`
before touching real data. The scoring modules themselves
(``g4watch/scoring/``) remain importable and unit-testable against
synthetic fixtures — that is Phase 1 infrastructure and is validated by
the non-circular ground-truth suite — but they are not wired to a real
corpus until this function stops raising.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from .validation.dh1_gate import Dh1Verdict

#: Ledger verdict written when the Appendix C floor stopped a locus before
#: the D.H1 test could run. Deliberately not a member of ``Dh1Verdict``:
#: it is the absence of a test, not one of its outcomes.
INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


class ScoringPermission(str, Enum):
    PERMITTED = "PERMITTED"
    BLOCKED_NO_LEDGER_ENTRY = "BLOCKED_NO_LEDGER_ENTRY"
    BLOCKED_INSUFFICIENT_DATA = "BLOCKED_INSUFFICIENT_DATA"
    BLOCKED_NOT_SUPPORTED = "BLOCKED_NOT_SUPPORTED"
    BLOCKED_SIGNAL_EXPLAINED_BY_GC = "BLOCKED_SIGNAL_EXPLAINED_BY_GC"
    #: A real, GC-adjusted effect exists but runs AGAINST D.H1 — the locus
    #: is more disrupted than its matched control. Blocked like any other
    #: non-SUPPORTED verdict, but named separately because it is evidence
    #: against the hypothesis rather than an absence of evidence for it.
    BLOCKED_SIGNAL_OPPOSITE_DIRECTION = "BLOCKED_SIGNAL_OPPOSITE_DIRECTION"
    BLOCKED_OPERATIONAL_MODE_OFF = "BLOCKED_OPERATIONAL_MODE_OFF"


@dataclass(frozen=True)
class GateStatus:
    """The full, reportable state of the D.H1 gate for one pathogen.

    Carries enough detail for the dashboard requirement in Section 17
    (gate status must be prominently displayed, never blank) without the
    caller having to re-read the ledger.
    """

    pathogen: str
    permission: ScoringPermission
    n_ledger_rows: int
    n_supported_loci: int
    supported_loci: tuple[str, ...]
    latest_timestamp: str | None
    failing_checks: tuple[str, ...]
    operational_mode: bool

    @property
    def permitted(self) -> bool:
        return self.permission is ScoringPermission.PERMITTED

    def explain(self) -> str:
        if self.permitted:
            return (
                f"{self.pathogen}: D.H1 SUPPORTED at {self.n_supported_loci} locus/loci "
                f"({', '.join(self.supported_loci)}) as of {self.latest_timestamp}, and "
                f"operational_mode is on. Stage 5 scoring is permitted."
            )
        reasons = {
            ScoringPermission.BLOCKED_NO_LEDGER_ENTRY: (
                "the testing ledger records no D.H1 result for this pathogen. The gate has not been "
                "run, and an unrun gate is not a passed gate."
            ),
            ScoringPermission.BLOCKED_INSUFFICIENT_DATA: (
                "every locus was halted at the Appendix C minimum-data floor before D.H1 could run "
                f"(failing checks: {', '.join(self.failing_checks) or 'see ledger'}). This is "
                "INSUFFICIENT_DATA, which is distinct from a NOT_SUPPORTED biological null: the test "
                "was never invoked."
            ),
            ScoringPermission.BLOCKED_NOT_SUPPORTED: (
                "D.H1 ran and returned NOT_SUPPORTED at every locus. Per Concept Paper v2 Section 6.6 "
                "this is a reportable negative result, and scoring stays unwired."
            ),
            ScoringPermission.BLOCKED_SIGNAL_EXPLAINED_BY_GC: (
                "D.H1's apparent signal did not survive the GC-confound gate "
                "(SIGNAL_EXPLAINED_BY_GC). Scoring stays unwired."
            ),
            ScoringPermission.BLOCKED_SIGNAL_OPPOSITE_DIRECTION: (
                "at least one locus shows a real, GC-adjustment-surviving effect that runs AGAINST "
                "D.H1 (SIGNAL_OPPOSITE_DIRECTION): it is MORE disrupted than its matched control, "
                "not less. This is evidence against the hypothesis, not an absence of evidence for "
                "it, and it is reported separately for that reason. Scoring stays unwired."
            ),
            ScoringPermission.BLOCKED_OPERATIONAL_MODE_OFF: (
                "D.H1 is SUPPORTED, but operational_mode is false in the pathogen config. Both "
                "conditions must hold before any deployment claim is made."
            ),
        }
        return f"{self.pathogen}: Stage 5 scoring is BLOCKED — {reasons[self.permission]}"


class ScoringNotPermittedError(RuntimeError):
    """Raised when Stage 5/6 is invoked on real data before the gate passes."""


def read_ledger(ledger_path: str | Path, pathogen: str, test: str = "D.H1") -> list[dict]:
    """Rows of the append-only testing ledger for one pathogen and test.

    A missing ledger file yields no rows rather than an error: that is
    the legitimate state of a study that has not run the gate yet, and
    it is handled as ``BLOCKED_NO_LEDGER_ENTRY`` downstream.
    """
    path = Path(ledger_path)
    if not path.exists():
        return []
    with open(path, newline="") as handle:
        return [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("pathogen") == pathogen and row.get("test") == test
        ]


def evaluate_gate(
    ledger_path: str | Path,
    pathogen: str,
    *,
    operational_mode: bool,
) -> GateStatus:
    """Decide, from the ledger alone, whether Stage 5 may run.

    Only the most recent run is considered. The ledger is append-only, so
    earlier runs stay on the record for auditability, but a stale
    ``SUPPORTED`` from a superseded corpus must not authorise scoring
    after a later run downgraded the verdict.
    """
    rows = read_ledger(ledger_path, pathogen)
    if not rows:
        return GateStatus(
            pathogen=pathogen,
            permission=ScoringPermission.BLOCKED_NO_LEDGER_ENTRY,
            n_ledger_rows=0,
            n_supported_loci=0,
            supported_loci=(),
            latest_timestamp=None,
            failing_checks=(),
            operational_mode=operational_mode,
        )

    latest_timestamp = max(row.get("timestamp", "") for row in rows)
    latest = [row for row in rows if row.get("timestamp", "") == latest_timestamp]

    verdicts = [row.get("verdict", "") for row in latest]
    supported = tuple(row["atlas_id"] for row in latest if row.get("verdict") == Dh1Verdict.SUPPORTED.value)
    failing = tuple(
        sorted({check for row in latest for check in (row.get("minimum_data_failing_checks") or "").split(";") if check})
    )

    if supported:
        permission = (
            ScoringPermission.PERMITTED if operational_mode else ScoringPermission.BLOCKED_OPERATIONAL_MODE_OFF
        )
    elif all(v == INSUFFICIENT_DATA for v in verdicts):
        permission = ScoringPermission.BLOCKED_INSUFFICIENT_DATA
    elif Dh1Verdict.SIGNAL_OPPOSITE_DIRECTION.value in verdicts:
        # Ranked first among the blocking verdicts for the same reason the
        # locus-level roll-up ranks it first: an effect pointing away from
        # the hypothesis is the most specific finding in the run, and the
        # one most easily lost if reported as a bare NOT_SUPPORTED.
        permission = ScoringPermission.BLOCKED_SIGNAL_OPPOSITE_DIRECTION
    elif Dh1Verdict.SIGNAL_EXPLAINED_BY_GC.value in verdicts:
        # A GC-explained signal is the more specific and more important
        # finding to report, so it wins over a bare NOT_SUPPORTED when
        # both appear among the loci of one run.
        permission = ScoringPermission.BLOCKED_SIGNAL_EXPLAINED_BY_GC
    else:
        permission = ScoringPermission.BLOCKED_NOT_SUPPORTED

    return GateStatus(
        pathogen=pathogen,
        permission=permission,
        n_ledger_rows=len(rows),
        n_supported_loci=len(supported),
        supported_loci=supported,
        latest_timestamp=latest_timestamp or None,
        failing_checks=failing,
        operational_mode=operational_mode,
    )


#: What a closed gate does.
#:
#: ``refuse``   -- raise, produce nothing. The default, and the behaviour
#:                every pathogen has unless its config says otherwise.
#: ``annotate`` -- let scoring proceed, and mark the whole result as
#:                non-authoritative with the gate's own reason attached.
#:
#: ``annotate`` exists because refusing produces no artifact at all, and a
#: reader then has nothing to look at and no sense of what the machinery
#: would have said. It is NOT a way to obtain a usable score from a corpus
#: that has not earned one: the reason travels with every result, the
#: warning classifier is capped, and ``authoritative`` stays False. A
#: number carrying its own disclaimer is safer than a number with no
#: disclaimer, and both are less safe than a refusal -- which is why
#: ``refuse`` remains the default and has to be opted out of per pathogen.
ON_BLOCK_REFUSE = "refuse"
ON_BLOCK_ANNOTATE = "annotate"
DEFAULT_ON_BLOCK = ON_BLOCK_REFUSE


def assert_scoring_permitted(
    ledger_path: str | Path,
    pathogen: str,
    *,
    operational_mode: bool,
    on_block: str = DEFAULT_ON_BLOCK,
) -> GateStatus:
    """Raise unless Stage 5 scoring is permitted for this pathogen.

    This is the call every Stage 5/6 real-data entry point must make
    first. It raises rather than returning a boolean so that a caller
    cannot ignore the result by accident.

    ``on_block="annotate"`` returns the closed status instead of raising,
    leaving the caller responsible for carrying the reason into its
    output. Callers that do this MUST mark the result non-authoritative;
    :func:`g4watch.pipeline.stage5_driver.run_stage5` does.
    """
    status = evaluate_gate(ledger_path, pathogen, operational_mode=operational_mode)
    if not status.permitted and on_block != ON_BLOCK_ANNOTATE:
        raise ScoringNotPermittedError(status.explain())
    return status
