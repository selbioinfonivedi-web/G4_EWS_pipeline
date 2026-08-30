"""Stage 5 — scoring. Fail-closed behind the D.H1 gate.

Build Architecture Section 16 is explicit that no Phase 3 code is written
speculatively in parallel with Phase 2, and Section 6 makes Stage 5
contingent on a passed D.H1 test per pathogen. This module is therefore
deliberately *not* an implementation of the score hierarchy. It is the
gate check plus the point at which that implementation will attach.

What already exists and is not blocked:

* ``g4watch/scoring/g4_ews_core.py`` (M3) and
  ``g4watch/scoring/integrated_score.py`` (M4) are implemented as Phase 1
  infrastructure and validated by the non-circular ground-truth suite
  against synthetic data. Importing and unit-testing them is fine.
* What is blocked is *wiring them to a real corpus* — which is what this
  entry point would do, and what :func:`run_stage5_scoring` refuses until
  the ledger says the gate passed.

When D.H1 returns SUPPORTED for a pathogen and its config sets
``operational_mode: true``, the work to do here is enumerated in
:data:`PHASE_3_SCOPE`; nothing below should be implemented before then.
"""

from __future__ import annotations

from ..config import PathogenConfig
from ..gating import GateStatus, assert_scoring_permitted

#: Section 10/13 deliverables that attach here once the gate opens. Kept
#: as data rather than as stub functions so that no half-written scoring
#: code can be called by accident.
PHASE_3_SCOPE = (
    "g4watch/scoring/weight_fitting.py — regularized by default (Section 10.3)",
    "g4watch/scoring/cusum.py — autocorrelation-aware calibration (Section 13.5)",
    "g4watch/scoring/ewma.py",
    "g4watch/validation/power_analysis.py — underpowered_analysis flag (Section 13.4)",
    "g4watch/validation/model_comparison.py — M1-M4 distinct term sets (Section 10)",
    "g4watch/warning/classifier.py",
    "mandatory Bayesian prior-sensitivity report (Section 10.3)",
)


def check_gate(config: PathogenConfig) -> GateStatus:
    """Raise :class:`~g4watch.gating.ScoringNotPermittedError` unless permitted."""
    return assert_scoring_permitted(
        config.ledger_path,
        config.pathogen,
        operational_mode=config.operational_mode,
    )


def run_stage5_scoring(config: PathogenConfig) -> GateStatus:
    """Entry point for Stage 5 on real data.

    Checks the gate first and, on the current FMDV evidence, stops there.
    The ``NotImplementedError`` below is only reachable once a real
    SUPPORTED verdict exists — at which point implementing
    :data:`PHASE_3_SCOPE` becomes in-scope work rather than speculative.
    """
    status = check_gate(config)
    raise NotImplementedError(
        f"{status.explain()}\n\n"
        "The D.H1 gate has opened, so Stage 5 is now in scope — but it is not implemented yet. "
        "Phase 3 work to do, in order:\n  - " + "\n  - ".join(PHASE_3_SCOPE)
    )
