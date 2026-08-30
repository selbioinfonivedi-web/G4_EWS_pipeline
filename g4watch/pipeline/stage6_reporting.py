"""Stage 6 — reporting. Gate status is always reportable; scores are not.

Two different things live behind Stage 6, and conflating them is exactly
what Section 17 warns against:

* The **surveillance report** (scores, CUSUM/EWMA alerts, warning levels)
  is Phase 3/4 work and is blocked by the same gate as Stage 5.
* The **gate-status report** is never blocked. A pathogen whose D.H1
  verdict is NOT_SUPPORTED, SIGNAL_EXPLAINED_BY_GC or INSUFFICIENT_DATA
  must show that status prominently — Section 17 requires the dashboard
  to display it rather than leave the scoring section blank, and a
  negative result is a reportable finding (Concept Paper v2 Section 6.6),
  not an absence of one.

So :func:`render_gate_status_report` always works, and
:func:`run_stage6_surveillance_report` fails closed.
"""

from __future__ import annotations

from ..config import PathogenConfig
from ..gating import GateStatus, evaluate_gate
from .stage5_scoring import check_gate


def gate_status(config: PathogenConfig) -> GateStatus:
    return evaluate_gate(
        config.ledger_path,
        config.pathogen,
        operational_mode=config.operational_mode,
    )


def render_gate_status_report(config: PathogenConfig) -> str:
    """A plain-text status report. Never blocked, never blank.

    This is the text the dashboard's always-visible gate panel renders,
    and the text a run prints when it stops at the gate.
    """
    status = gate_status(config)
    banner = "SCORING PERMITTED" if status.permitted else "SCORING BLOCKED"
    lines = [
        f"{'=' * 72}",
        f"G4-WATCH — {config.display_name} ({config.pathogen}) — D.H1 GATE: {banner}",
        f"{'=' * 72}",
        f"  permission        : {status.permission.value}",
        f"  operational_mode  : {status.operational_mode}",
        f"  ledger rows       : {status.n_ledger_rows}",
        f"  supported loci    : {status.n_supported_loci}"
        + (f" ({', '.join(status.supported_loci)})" if status.supported_loci else ""),
        f"  latest run        : {status.latest_timestamp or 'never'}",
    ]
    if status.failing_checks:
        lines.append(f"  failing checks    : {', '.join(status.failing_checks)}")
    lines += ["", status.explain(), ""]
    if not status.permitted:
        lines += [
            "This is a reportable result, not a missing one. Per Concept Paper v2 Section 6.6 a",
            "negative or insufficient-data outcome is published as such; the surveillance score",
            "section stays unwired rather than being filled with numbers the evidence does not",
            "support.",
            "",
        ]
    return "\n".join(lines)


def run_stage6_surveillance_report(config: PathogenConfig) -> GateStatus:
    """Entry point for the scored surveillance report. Blocked by the gate."""
    status = check_gate(config)
    raise NotImplementedError(
        f"{status.explain()}\n\n"
        "Stage 6's scored surveillance report depends on Stage 5, which is not implemented "
        "(see g4watch/pipeline/stage5_scoring.py). The gate-status report is available now via "
        "`g4watch gate-status`."
    )
