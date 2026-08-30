"""Warning-level classification from a calibrated alert stream.

Translates control-chart output into the small set of levels an operator
can actually act on. Three properties are enforced here rather than left
to the caller, because each is a way a surveillance dashboard can mislead:

**A warning level is never issued without its evidence.** Every
:class:`WarningAssessment` carries the gate status, the underpowered
flag, and the confidence axes of the locus it concerns. A level shown
without them invites the reader to supply their own, more confident,
interpretation.

**An underpowered analysis cannot produce a high warning level.** If the
comparison could not have detected the effect, an alarm from it is not
evidence of the effect. Such an assessment is capped at
``INSUFFICIENT_EVIDENCE`` no matter what the chart did.

**Structural confidence caps the level.** A locus that is a weak
computational candidate (``WC``) with a single concordant tool cannot
raise a HIGH warning, however cleanly its chart alarms — the underlying
claim that there is a G4 there at all is not strong enough to carry it.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from ..atlas.schema import StructuralConfidence


class WarningLevel(str, Enum):
    """Ordered from least to most actionable."""

    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"
    NONE = "NONE"
    WATCH = "WATCH"
    ELEVATED = "ELEVATED"
    HIGH = "HIGH"


#: Rank for cap comparisons.
_ORDER = {
    WarningLevel.INSUFFICIENT_EVIDENCE: 0,
    WarningLevel.NONE: 1,
    WarningLevel.WATCH: 2,
    WarningLevel.ELEVATED: 3,
    WarningLevel.HIGH: 4,
}

#: The highest level each structural-confidence class may reach. EC/BC are
#: experimentally or biophysically confirmed; SC/MC are computational; WC
#: is weak and AA is an algorithm artefact that must never warn at all.
CONFIDENCE_CAP = {
    StructuralConfidence.EC: WarningLevel.HIGH,
    StructuralConfidence.BC: WarningLevel.HIGH,
    StructuralConfidence.SC: WarningLevel.ELEVATED,
    StructuralConfidence.MC: WarningLevel.WATCH,
    StructuralConfidence.WC: WarningLevel.WATCH,
    StructuralConfidence.AA: WarningLevel.NONE,
}


@dataclass(frozen=True)
class WarningAssessment:
    """A warning level together with everything needed to read it honestly."""

    atlas_id: str
    level: WarningLevel
    uncapped_level: WarningLevel
    n_consecutive_alarms: int
    structural_confidence: StructuralConfidence
    underpowered: bool
    gate_permitted: bool
    capped_by: str | None

    @property
    def actionable(self) -> bool:
        return _ORDER[self.level] >= _ORDER[WarningLevel.ELEVATED]

    def summary(self) -> str:
        lines = [f"{self.atlas_id}: {self.level.value}"]
        if self.capped_by:
            lines.append(f"  (chart alone indicated {self.uncapped_level.value}; capped by {self.capped_by})")
        lines.append(
            f"  evidence: {self.n_consecutive_alarms} consecutive alarm window(s), "
            f"structural confidence {self.structural_confidence.name}, "
            f"underpowered={self.underpowered}, D.H1 gate permitted={self.gate_permitted}"
        )
        return "\n".join(lines)


def _level_from_alarms(n_consecutive_alarms: int) -> WarningLevel:
    """Consecutive alarm windows -> raw level.

    Consecutive windows rather than a single crossing: one alarm in an
    autocorrelated series is weak evidence, and requiring persistence is
    the cheapest available guard against acting on a transient.
    """
    if n_consecutive_alarms <= 0:
        return WarningLevel.NONE
    if n_consecutive_alarms == 1:
        return WarningLevel.WATCH
    if n_consecutive_alarms == 2:
        return WarningLevel.ELEVATED
    return WarningLevel.HIGH


def classify_warning(
    atlas_id: str,
    *,
    n_consecutive_alarms: int,
    structural_confidence: StructuralConfidence,
    underpowered: bool,
    gate_permitted: bool,
) -> WarningAssessment:
    """Classify one locus's warning level, applying every cap.

    Caps are applied in order of severity so that ``capped_by`` names the
    binding constraint rather than the last one checked.
    """
    uncapped = _level_from_alarms(n_consecutive_alarms)
    level = uncapped
    capped_by: str | None = None

    if not gate_permitted:
        # No scored warning is meaningful before the D.H1 gate passes.
        level = WarningLevel.INSUFFICIENT_EVIDENCE
        capped_by = "the D.H1 gate (scoring is not permitted for this pathogen)"
    elif underpowered:
        level = WarningLevel.INSUFFICIENT_EVIDENCE
        capped_by = "an underpowered analysis (an alarm from a test that could not have detected the effect is not evidence of it)"
    else:
        cap = CONFIDENCE_CAP[structural_confidence]
        if _ORDER[level] > _ORDER[cap]:
            level = cap
            capped_by = f"structural confidence {structural_confidence.name} (cap: {cap.value})"

    return WarningAssessment(
        atlas_id=atlas_id,
        level=level,
        uncapped_level=uncapped,
        n_consecutive_alarms=n_consecutive_alarms,
        structural_confidence=structural_confidence,
        underpowered=underpowered,
        gate_permitted=gate_permitted,
        capped_by=capped_by,
    )


def count_trailing_alarms(alarm_indices: tuple[int, ...], n_windows: int) -> int:
    """Consecutive alarm windows ending at the most recent window.

    Counts backwards from the end: a burst of alarms two years ago is not
    a current warning, and only the run that reaches the present window
    counts towards the level.
    """
    if not alarm_indices or n_windows <= 0:
        return 0
    alarms = set(alarm_indices)
    count = 0
    for index in range(n_windows - 1, -1, -1):
        if index in alarms:
            count += 1
        else:
            break
    return count
