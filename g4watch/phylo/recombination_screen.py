"""Recombination screening via the PHI test (Bruen, Philippe & Bryant 2006),
mandatory before ancestral-state reconstruction (architecture Section 11),
tiered by pathogen: `high_priority` for LSDV (documented recombinant
vaccine/wild-type mixing), `standard` for the four RNA viruses.

Wraps the real `Phi` binary (vendored at vendor/phipack/, built from
source — see vendor/README.md for provenance) as a subprocess; does not
reimplement the PHI test itself (unlike G4Hunter, this is not a simple
enough algorithm to safely reimplement from a paper description — it is a
nontrivial pairwise-incompatibility statistic with its own normal-approximation
and permutation-test machinery, exactly the kind of tool this project's own
house rule says to orchestrate, not reimplement).
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

_VENDOR_PHI_BINARY = Path(__file__).resolve().parents[2] / "vendor" / "phipack" / "Phi"

_FOUND_SEQ_RE = re.compile(r"Found (\d+) sequences of length (\d+)")
_INFORMATIVE_SITES_RE = re.compile(r"Found (\d+) informative sites")
_PHI_NORMAL_RE = re.compile(r"PHI \(Normal\):\s*(--|[\d.eE+-]+)")
_TOO_FEW_SITES_RE = re.compile(r"Too few informative sites")


class RecombinationTier(str, Enum):
    HIGH_PRIORITY = "high_priority"
    STANDARD = "standard"


@dataclass(frozen=True)
class PhiTestResult:
    n_sequences: int
    alignment_length: int
    n_informative_sites: int
    phi_p_value: float | None  # None when PhiPack reports "--" (too few informative sites)
    undefined_reason: str | None
    raw_stdout: str


@dataclass(frozen=True)
class RecombinationScreenResult:
    phi: PhiTestResult
    tier: RecombinationTier
    significant: bool
    route_to_non_tree_based_estimate: bool


class PhiExecutionError(RuntimeError):
    """Raised when the Phi binary itself fails or is missing — never
    silently swallowed, per this project's no-silent-fallback rule."""


def parse_phi_output(stdout: str) -> PhiTestResult:
    """Pure parsing logic, independently unit-testable against captured
    real output text without invoking the binary."""
    seq_match = _FOUND_SEQ_RE.search(stdout)
    if seq_match is None:
        raise PhiExecutionError(f"Could not parse sequence count from Phi output:\n{stdout}")
    n_sequences, alignment_length = int(seq_match.group(1)), int(seq_match.group(2))

    info_match = _INFORMATIVE_SITES_RE.search(stdout)
    n_informative = int(info_match.group(1)) if info_match else 0

    phi_match = _PHI_NORMAL_RE.search(stdout)
    if phi_match is None:
        raise PhiExecutionError(f"Could not find a 'PHI (Normal):' line in Phi output:\n{stdout}")

    raw_value = phi_match.group(1)
    if raw_value == "--":
        undefined_reason = (
            "too few informative sites for the normal approximation"
            if _TOO_FEW_SITES_RE.search(stdout)
            else "PHI value undefined"
        )
        return PhiTestResult(
            n_sequences=n_sequences,
            alignment_length=alignment_length,
            n_informative_sites=n_informative,
            phi_p_value=None,
            undefined_reason=undefined_reason,
            raw_stdout=stdout,
        )

    return PhiTestResult(
        n_sequences=n_sequences,
        alignment_length=alignment_length,
        n_informative_sites=n_informative,
        phi_p_value=float(raw_value),
        undefined_reason=None,
        raw_stdout=stdout,
    )


def run_phi_test(fasta_path: str | Path, phi_binary: str | Path = _VENDOR_PHI_BINARY) -> PhiTestResult:
    binary = Path(phi_binary)
    if not binary.exists():
        raise PhiExecutionError(f"Phi binary not found at {binary} — build it via `make` in vendor/phipack/ first.")

    resolved_fasta = Path(fasta_path).resolve()
    result = subprocess.run(
        [str(binary), "-f", str(resolved_fasta)],
        capture_output=True,
        text=True,
        cwd=resolved_fasta.parent,  # Phi writes side-effect files (Phi.inf.sites etc.) into cwd
    )
    if result.returncode != 0:
        raise PhiExecutionError(
            f"Phi exited with code {result.returncode} on {fasta_path}:\n{result.stdout}\n{result.stderr}"
        )
    return parse_phi_output(result.stdout)


def screen_recombination(
    fasta_path: str | Path,
    tier: RecombinationTier,
    alpha: float = 0.05,
    phi_binary: str | Path = _VENDOR_PHI_BINARY,
) -> RecombinationScreenResult:
    """Runs the PHI test unconditionally (the test always runs — `tier`
    only controls how a positive/undefined result is handled downstream,
    per architecture Section 11). A significant result under
    `HIGH_PRIORITY` routes affected loci to a non-tree-based, descriptive
    conservation estimate (recombination corrupts ancestral-state
    reconstruction); under `STANDARD` it is logged for the dashboard's
    data-quality section and flagged for domain-reviewer sign-off, but does
    not automatically exclude the locus."""
    phi_result = run_phi_test(fasta_path, phi_binary)
    significant = phi_result.phi_p_value is not None and phi_result.phi_p_value < alpha

    return RecombinationScreenResult(
        phi=phi_result,
        tier=tier,
        significant=significant,
        route_to_non_tree_based_estimate=(significant and tier == RecombinationTier.HIGH_PRIORITY),
    )
