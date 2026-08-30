"""Appendix C minimum-data floor (architecture Section 14): an absolute,
unconditional data-sufficiency check that runs BEFORE `dh1_gate`, not
instead of it and not replaced by the power analysis (Section 13.4) --
passing this floor never implies a pathogen/window is adequately powered,
and failing power analysis never excuses failing this floor. The two are
always reported together, per Section 14.

Pulled forward from its originally-scheduled Sprint (13) because Sprint 8's
own full-pipeline integration test needs it to distinguish `INSUFFICIENT_DATA`
(this gate fires) from `NOT_SUPPORTED` (dh1_gate ran and found no effect) --
these are different outcomes the dashboard must never conflate (Section 15).

The architecture text states the floor's three headline numbers (>=30
sequences/window, >=20/lineage, >=3 timepoints, v1 Part I.5 / Appendix C)
and says "9 checks" without enumerating all of them. Rather than invent
numbers for the unstated remainder, the 9 checks below are the headline
three plus six more pipeline-completeness preconditions this project's own
downstream stages already depend on implicitly (a matched control region,
minimally-informative clade counts for both gc_confound_gate and dh1_gate,
corpus alignment/QC quality, metadata completeness for temporal binning,
and recombination screening having actually been run) -- each individually
hard-fails to INSUFFICIENT_DATA, exactly like the stated three. This is a
documented judgment call, not a literature-derived or architecture-specified
set, consistent with how this codebase already handles underspecified
points elsewhere (see gc_confound_gate.py's and dh1_gate.py's own
design-history comments)."""

from __future__ import annotations

from dataclasses import dataclass

MIN_SEQUENCES_PER_WINDOW = 30
MIN_SEQUENCES_PER_LINEAGE = 20
MIN_TIMEPOINTS = 3
MIN_METADATA_COMPLETENESS_FRACTION = 0.90
MIN_INFORMATIVE_CLADES_PER_GROUP = 3  # matches gc_confound_gate.DEFAULT_MIN_CLADES_PER_GROUP
MIN_ALIGNMENT_QC_PASS_FRACTION = 0.50


@dataclass(frozen=True)
class MinimumDataInput:
    n_sequences_in_window: int
    min_sequences_per_lineage: int  # smallest per-lineage count across all lineages in the window
    n_timepoints: int
    metadata_completeness_fraction: float  # fraction with complete date/host/country
    n_locus_informative_clades: int
    n_control_informative_clades: int
    control_region_found: bool
    alignment_qc_pass_fraction: float  # fraction of the raw corpus surviving QC + alignment
    recombination_screen_completed: bool


@dataclass(frozen=True)
class MinimumDataResult:
    passed_minimum_floor: bool
    failing_checks: tuple[str, ...]
    checks: dict[str, bool]


def minimum_data_gate(data: MinimumDataInput) -> MinimumDataResult:
    """Runs all 9 checks unconditionally (never short-circuits on the first
    failure) so a caller sees every reason a pathogen/window is data-
    insufficient at once, not just the first one hit."""
    checks: dict[str, bool] = {
        "n_sequences_in_window": data.n_sequences_in_window >= MIN_SEQUENCES_PER_WINDOW,
        "min_sequences_per_lineage": data.min_sequences_per_lineage >= MIN_SEQUENCES_PER_LINEAGE,
        "n_timepoints": data.n_timepoints >= MIN_TIMEPOINTS,
        "metadata_completeness": data.metadata_completeness_fraction >= MIN_METADATA_COMPLETENESS_FRACTION,
        "n_locus_informative_clades": data.n_locus_informative_clades >= MIN_INFORMATIVE_CLADES_PER_GROUP,
        "n_control_informative_clades": data.n_control_informative_clades >= MIN_INFORMATIVE_CLADES_PER_GROUP,
        "control_region_found": data.control_region_found,
        "alignment_qc_pass_fraction": data.alignment_qc_pass_fraction >= MIN_ALIGNMENT_QC_PASS_FRACTION,
        "recombination_screen_completed": data.recombination_screen_completed,
    }
    failing = tuple(name for name, ok in checks.items() if not ok)
    return MinimumDataResult(passed_minimum_floor=len(failing) == 0, failing_checks=failing, checks=checks)
