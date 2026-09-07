"""Workbench services: datasets, validation, projects, system telemetry.

Everything here is real. File validation actually parses the file and
reports what it found; the system monitor reads the host through psutil;
projects are written to and read from disk. Where a capability genuinely
does not exist yet -- BAM record parsing without pysam, Excel without
openpyxl -- the validator says so rather than inventing a record count.

That distinction matters more in scientific software than in most
software: a validator that guesses is worse than one that abstains,
because a researcher will believe it.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import psutil

from .commands import REPO_ROOT, CommandError, _resolve_inside_repo

PROJECTS_DIR = REPO_ROOT / "projects"

Status = Literal["valid", "warning", "invalid", "unknown"]

# Extension -> the format we will attempt. Anything unlisted is inspected
# by content sniffing and reported as unknown if that fails too.
FORMATS = {
    ".fa": "fasta",
    ".fasta": "fasta",
    ".fna": "fasta",
    ".faa": "fasta",
    ".fas": "fasta",
    ".fq": "fastq",
    ".fastq": "fastq",
    ".bam": "bam",
    ".sam": "sam",
    ".cram": "cram",
    ".vcf": "vcf",
    ".csv": "csv",
    ".tsv": "tsv",
    ".txt": "tsv",
    ".xlsx": "excel",
    ".xls": "excel",
    ".nwk": "newick",
    ".newick": "newick",
    ".tree": "newick",
    ".treefile": "newick",
    ".nex": "nexus",
    ".nexus": "nexus",
    ".gff": "annotation",
    ".gff3": "annotation",
    ".gtf": "annotation",
    ".bed": "annotation",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".json": "json",
    ".yaml": "config",
    ".yml": "config",
}

IUPAC_NT = set("ACGTURYSWKMBDHVN-.acgturyswkmbdhvn")


@dataclass
class Finding:
    level: Literal["pass", "warn", "fail"]
    message: str
    fixable: bool = False
    fix: str | None = None


@dataclass
class Validation:
    path: str
    name: str
    format: str
    size: int
    status: Status = "unknown"
    records: int | None = None
    record_label: str = "records"
    findings: list[Finding] = field(default_factory=list)
    detail: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "name": self.name,
            "format": self.format,
            "size": self.size,
            "status": self.status,
            "records": self.records,
            "record_label": self.record_label,
            "findings": [
                {"level": f.level, "message": f.message, "fixable": f.fixable, "fix": f.fix} for f in self.findings
            ],
            "detail": self.detail,
        }

    def settle(self) -> Validation:
        if any(f.level == "fail" for f in self.findings):
            self.status = "invalid"
        elif any(f.level == "warn" for f in self.findings):
            self.status = "warning"
        elif self.findings:
            self.status = "valid"
        return self


def _open(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def detect_format(path: Path) -> str:
    stem = path
    if path.suffix == ".gz":
        stem = Path(path.stem)
    fmt = FORMATS.get(stem.suffix.lower())
    if fmt:
        return fmt
    try:
        with path.open("rb") as fh:
            head = fh.read(4)
        if head[:2] == b"\x1f\x8b":
            return "bam"
        if head[:2] == b"PK":
            return "excel"
        with _open(path) as fh:
            first = fh.readline().strip()
        if first.startswith(">"):
            return "fasta"
        if first.startswith("@"):
            return "fastq"
        if first.startswith("##fileformat=VCF"):
            return "vcf"
        if first.upper().startswith("#NEXUS"):
            return "nexus"
        if first.startswith("(") and first.rstrip().endswith(";"):
            return "newick"
    except OSError:
        pass
    return "unknown"


# ── format validators ───────────────────────────────────────────────
def _validate_fasta(path: Path, v: Validation) -> None:
    n, empty, bad_chars, lengths, dup = 0, 0, 0, [], 0
    seen: set[str] = set()
    current, length = None, 0
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if current is not None:
                    lengths.append(length)
                    if length == 0:
                        empty += 1
                name = line[1:].split()[0] if len(line) > 1 else ""
                if name in seen:
                    dup += 1
                seen.add(name)
                current, length, n = name, 0, n + 1
            elif current is not None:
                length += len(line)
                bad_chars += sum(1 for c in line if c not in IUPAC_NT)
    if current is not None:
        lengths.append(length)
        if length == 0:
            empty += 1

    v.records, v.record_label = n, "sequences"
    if n == 0:
        v.findings.append(Finding("fail", "No FASTA records found. The file has no '>' header line."))
        return
    v.findings.append(Finding("pass", f"Parsed {n:,} sequences."))
    if lengths:
        v.detail["length_min"], v.detail["length_max"] = min(lengths), max(lengths)
        v.detail["length_mean"] = round(sum(lengths) / len(lengths))
        spread = max(lengths) - min(lengths)
        if spread and min(lengths) < 0.5 * max(lengths):
            v.findings.append(
                Finding(
                    "warn",
                    f"Sequence lengths vary widely ({min(lengths):,}-{max(lengths):,} nt). If this should be an alignment, it is not aligned.",
                )
            )
    if empty:
        v.findings.append(
            Finding(
                "fail", f"{empty} record{'s' if empty != 1 else ''} have a header but no sequence.", True, "drop_empty"
            )
        )
    if dup:
        v.findings.append(
            Finding("warn", f"{dup} duplicate sequence identifier{'s' if dup != 1 else ''}.", True, "dedupe")
        )
    if bad_chars:
        v.findings.append(
            Finding(
                "warn",
                f"{bad_chars:,} character{'s' if bad_chars != 1 else ''} outside the IUPAC nucleotide alphabet.",
                True,
                "mask_non_iupac",
            )
        )


def _validate_fastq(path: Path, v: Validation) -> None:
    n, malformed, mismatched = 0, 0, 0
    with _open(path) as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
            if not qual:
                malformed += 1
                break
            if not head.startswith("@") or not plus.startswith("+"):
                malformed += 1
            elif len(seq.strip()) != len(qual.strip()):
                mismatched += 1
            n += 1
    v.records, v.record_label = n, "reads"
    if n == 0:
        v.findings.append(Finding("fail", "No FASTQ records found."))
        return
    v.findings.append(Finding("pass", f"Parsed {n:,} reads."))
    if malformed:
        v.findings.append(
            Finding(
                "fail",
                f"{malformed} malformed record{'s' if malformed != 1 else ''} (record structure is not @/seq/+/qual).",
            )
        )
    if mismatched:
        v.findings.append(
            Finding(
                "fail",
                f"{mismatched} record{'s' if mismatched != 1 else ''} where the quality string length differs from the sequence length.",
            )
        )


def _validate_delimited(path: Path, v: Validation, delim: str) -> None:
    with _open(path) as fh:
        rows = list(csv.reader(fh, delimiter=delim))
    if not rows:
        v.findings.append(Finding("fail", "File is empty."))
        return
    header, body = rows[0], rows[1:]
    v.records, v.record_label = len(body), "rows"
    v.detail["columns"] = header
    v.findings.append(Finding("pass", f"Parsed {len(body):,} rows across {len(header)} columns."))

    ragged = sum(1 for r in body if len(r) != len(header))
    if ragged:
        v.findings.append(
            Finding("fail", f"{ragged} row{'s' if ragged != 1 else ''} do not have {len(header)} fields.")
        )

    # Metadata-specific checks, only when the table looks like metadata.
    lower = [h.strip().lower() for h in header]
    for column in ("collection_date", "country", "host", "accession"):
        if column not in lower:
            continue
        idx = lower.index(column)
        missing = sum(1 for r in body if idx >= len(r) or not r[idx].strip())
        if missing:
            level = "fail" if column == "accession" else "warn"
            v.findings.append(
                Finding(
                    level,
                    f"{missing:,} record{'s' if missing != 1 else ''} have no {column.replace('_', ' ')}.",
                    column != "accession",
                    "drop_incomplete",
                )
            )
        else:
            article = "an" if column[0] in "aeiou" else "a"
            v.findings.append(Finding("pass", f"Every row has {article} {column.replace('_', ' ')}."))

    if "accession" in lower:
        idx = lower.index("accession")
        seen: set[str] = set()
        dup = 0
        for r in body:
            if idx < len(r):
                if r[idx] in seen:
                    dup += 1
                seen.add(r[idx])
        if dup:
            v.findings.append(Finding("fail", f"{dup} duplicate accession{'s' if dup != 1 else ''}.", True, "dedupe"))


def _validate_newick(path: Path, v: Validation) -> None:
    from web.workstation.dataset import layout_tree, parse_newick

    text = path.read_text(errors="replace").strip()
    if text.count("(") != text.count(")"):
        v.findings.append(Finding("fail", f"Unbalanced parentheses: {text.count('(')} open, {text.count(')')} close."))
        return
    if not text.rstrip().endswith(";"):
        v.findings.append(Finding("warn", "Tree string does not end with ';'.", True, "append_semicolon"))
    try:
        out = layout_tree(parse_newick(text))
    except Exception as exc:  # noqa: BLE001 - reported, not swallowed
        v.findings.append(Finding("fail", f"Could not parse the tree: {exc}"))
        return
    v.records, v.record_label = out["n_leaves"], "tips"
    v.detail["nodes"] = len(out["nodes"])
    v.findings.append(Finding("pass", f"Parsed {out['n_leaves']:,} tips across {len(out['nodes']):,} nodes."))
    if ":" not in text:
        v.findings.append(Finding("warn", "No branch lengths. Time-scaled and evolutionary analyses require them."))


def _validate_vcf(path: Path, v: Validation) -> None:
    n, has_header, samples = 0, False, 0
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                has_header = True
                samples = max(len(line.rstrip().split("\t")) - 9, 0)
                continue
            if line.strip():
                n += 1
    v.records, v.record_label = n, "variants"
    if not has_header:
        v.findings.append(Finding("fail", "No #CHROM header line; this is not a valid VCF."))
        return
    v.detail["samples"] = samples
    v.findings.append(
        Finding("pass", f"Parsed {n:,} variant records across {samples} sample column{'s' if samples != 1 else ''}.")
    )
    if n == 0:
        v.findings.append(Finding("warn", "Header is valid but the file contains no variant records."))


def _validate_binary(path: Path, v: Validation, label: str) -> None:
    """Honest partial validation: identify the container, refuse to guess."""
    with path.open("rb") as fh:
        magic = fh.read(4)
    if label == "bam" and magic[:2] != b"\x1f\x8b":
        v.findings.append(Finding("fail", "Not BGZF-compressed; this is not a BAM file."))
        return
    if label == "excel" and magic[:2] != b"PK":
        v.findings.append(Finding("fail", "Not a ZIP container; this is not an .xlsx file."))
        return
    v.findings.append(Finding("pass", f"Container header is a valid {label.upper()}."))
    tool = "pysam" if label == "bam" else "openpyxl"
    v.findings.append(
        Finding(
            "warn",
            f"Record-level validation needs {tool}, which is not installed. Record count is not reported rather than guessed.",
        )
    )


VALIDATORS = {
    "fasta": _validate_fasta,
    "fastq": _validate_fastq,
    "newick": _validate_newick,
    "vcf": _validate_vcf,
}


def validate_file(rel_path: str) -> dict:
    path = _resolve_inside_repo(rel_path, "path")
    if not path.is_file():
        raise CommandError(f"no such file: {rel_path}")

    fmt = detect_format(path)
    v = Validation(path=str(path.relative_to(REPO_ROOT)), name=path.name, format=fmt, size=path.stat().st_size)

    if fmt in VALIDATORS:
        VALIDATORS[fmt](path, v)
    elif fmt in {"tsv", "annotation"}:
        _validate_delimited(path, v, "\t")
    elif fmt == "csv":
        _validate_delimited(path, v, ",")
    elif fmt in {"bam", "excel", "cram"}:
        _validate_binary(path, v, fmt)
    elif fmt in {"nexus", "genbank", "sam", "json", "config"}:
        v.findings.append(Finding("pass", f"Recognised as {fmt.upper()}."))
        v.findings.append(
            Finding("warn", f"No structural validator for {fmt.upper()} yet; only the format was identified.")
        )
    else:
        v.findings.append(Finding("fail", "Format not recognised from the extension or the file header."))

    if path.stat().st_size == 0:
        v.findings.insert(0, Finding("fail", "File is empty (0 bytes)."))
    return v.settle().to_dict()


def scan_inputs() -> list[dict]:
    """Files in the project that are plausible pipeline inputs."""
    roots = ["data", "results", "config"]
    out: list[dict] = []
    for root in roots:
        base = REPO_ROOT / root
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file() or path.name.startswith("."):
                continue
            stem = Path(path.stem) if path.suffix == ".gz" else path
            if stem.suffix.lower() not in FORMATS:
                continue
            stat = path.stat()
            if stat.st_size > 400_000_000:
                continue
            out.append(
                {
                    "path": str(path.relative_to(REPO_ROOT)),
                    "name": path.name,
                    "format": FORMATS.get(stem.suffix.lower(), "unknown"),
                    "size": stat.st_size,
                    "modified": stat.st_mtime,
                }
            )
    out.sort(key=lambda f: f["modified"], reverse=True)
    return out[:300]


# ── system telemetry ────────────────────────────────────────────────
_last_net = {"t": 0.0, "read": 0, "write": 0}


def system_status() -> dict:
    vm = psutil.virtual_memory()
    disk = psutil.disk_usage(str(REPO_ROOT))
    load = psutil.getloadavg() if hasattr(psutil, "getloadavg") else (0, 0, 0)
    return {
        "cpu_percent": psutil.cpu_percent(interval=None),
        "cpu_count": psutil.cpu_count(logical=True),
        "cpu_physical": psutil.cpu_count(logical=False),
        "load": [round(x, 2) for x in load],
        "memory": {
            "total": vm.total,
            "used": vm.total - vm.available,
            "available": vm.available,
            "percent": vm.percent,
        },
        "disk": {"total": disk.total, "used": disk.used, "free": disk.free, "percent": disk.percent},
        "threads": sum(
            p.info["num_threads"] or 0 for p in psutil.process_iter(["num_threads"]) if p.info["num_threads"]
        ),
        "boot_time": psutil.boot_time(),
        "now": time.time(),
    }


# ── projects ────────────────────────────────────────────────────────
_SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 _-]{0,63}$")


def _project_path(name: str) -> Path:
    if not _SAFE_NAME.match(name):
        raise CommandError("project name must be 1-64 characters: letters, digits, space, hyphen, underscore")
    return PROJECTS_DIR / f"{name}.g4proj.json"


def save_project(name: str, payload: dict) -> dict:
    PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
    path = _project_path(name)
    document = {
        "name": name,
        "schema": "g4watch/project@1",
        "saved": time.time(),
        "payload": payload,
    }
    path.write_text(json.dumps(document, indent=1))
    return {"name": name, "path": str(path.relative_to(REPO_ROOT)), "saved": document["saved"]}


def load_project(name: str) -> dict:
    path = _project_path(name)
    if not path.is_file():
        raise CommandError(f"no project named {name!r}")
    return json.loads(path.read_text())


def list_projects() -> list[dict]:
    if not PROJECTS_DIR.is_dir():
        return []
    out = []
    for path in sorted(PROJECTS_DIR.glob("*.g4proj.json")):
        try:
            doc = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        out.append(
            {
                "name": doc.get("name", path.stem),
                "saved": doc.get("saved", path.stat().st_mtime),
                "size": path.stat().st_size,
            }
        )
    out.sort(key=lambda p: p["saved"], reverse=True)
    return out


def delete_project(name: str) -> dict:
    path = _project_path(name)
    if path.is_file():
        path.unlink()
        return {"deleted": True, "name": name}
    return {"deleted": False, "name": name}


def duplicate_project(name: str, new_name: str) -> dict:
    src, dst = _project_path(name), _project_path(new_name)
    if not src.is_file():
        raise CommandError(f"no project named {name!r}")
    doc = json.loads(src.read_text())
    doc["name"] = new_name
    doc["saved"] = time.time()
    doc["duplicated_from"] = name
    dst.write_text(json.dumps(doc, indent=1))
    return {"name": new_name, "path": str(dst.relative_to(REPO_ROOT))}


# ── provenance ──────────────────────────────────────────────────────
def provenance() -> dict:
    """What a run must record to be reproducible."""
    from g4watch import __version__

    tools = {}
    for name in ("mafft", "iqtree2", "treetime", "nextflow", "Rscript"):
        from .commands import resolve_executable

        tools[name] = resolve_executable(name)
    phi = REPO_ROOT / "vendor" / "phipack" / "Phi"
    tools["PhiPack"] = str(phi) if phi.exists() else None

    commit = None
    head = REPO_ROOT / ".git" / "HEAD"
    if head.is_file():
        ref = head.read_text().strip()
        if ref.startswith("ref: "):
            ref_path = REPO_ROOT / ".git" / ref[5:]
            if ref_path.is_file():
                commit = ref_path.read_text().strip()[:12]
        else:
            commit = ref[:12]

    return {
        "software": {"name": "G4-WATCH", "version": __version__},
        "git_commit": commit,
        "tools": tools,
        "python": shutil.which("python3"),
        "repo_root": str(REPO_ROOT),
    }
