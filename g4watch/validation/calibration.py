"""Calibrate the structural-confidence operating point against known G4s.

WHAT THIS IS FOR. The SC rule -- two concordant tools, |G4Hunter| >= 1.5,
conservation >= 85% -- was carried forward from Revision 1 and never
re-derived for cross-species use. Open question Q1.1 recorded that and
concluded calibration was impossible because no FMDV G4 has been
biophysically confirmed.

That reasoning was right about FMDV and wrong about the framework.
Confirmed G4s exist in other viruses and are a legitimate target for a
tool that claims to work across species. Measured against them on
2026-09-09, the current rule classified every confirmed G4 tested as WC:
the one tier that can never be scoring-eligible.

    HIV-1 LTR-III/IV (5' U3)   |1.107|  1 tool  -> WC
    HIV-1 nef                  |0.944|  1 tool  -> WC
    HIV-1 LTR-III/IV (3' U3)   |1.154|  1 tool  -> WC

WHAT THIS MODULE DOES NOT DO. It does not change a threshold. It reports
what the current one achieves and what alternatives would achieve, and
stops. Moving an operating point is a scientific decision that needs a
person to make it and a revision-log entry to record it -- and moving one
while looking at the loci it would admit is precisely the failure this
whole framework is built to prevent.

WHY POSITIVES AND NEGATIVES ARE MATCHED THE SAME WAY. Negatives come from
``validation/control_regions.py`` -- the same GC- and length-matching the
pipeline already uses to pick controls for D.H1. Matching them any other
way would calibrate against a contrast the pipeline never actually draws.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from ..g4prediction import g4hunter, pattern_motif

#: A row whose coordinates were located computationally rather than read
#: from the publication. Calibrating on these is mildly circular: the span
#: sits where the predictor already looked. Reported separately.
PROVENANCE_DERIVED = "derived"
PROVENANCE_STATED = "stated"

#: Below this many positives, an operating point derived here is not
#: usable. Stated rather than assumed so the refusal is visible.
MIN_POSITIVES_FOR_CALIBRATION = 30

#: And from fewer than this many distinct viruses, a threshold has learned
#: the viruses rather than the biology.
MIN_VIRUSES_FOR_CALIBRATION = 4


@dataclass(frozen=True)
class ScoredLocus:
    locus_id: str
    virus: str
    is_positive: bool
    g4hunter_score: float | None
    n_tools: int
    provenance: str = ""

    @property
    def magnitude(self) -> float:
        """|G4Hunter|. The sign is the strand, not the quality."""
        return abs(self.g4hunter_score) if self.g4hunter_score is not None else 0.0


@dataclass(frozen=True)
class CalibrationReport:
    positives: list[ScoredLocus]
    negatives: list[ScoredLocus]
    curve: list[dict]
    usable: bool
    blocking_reasons: tuple[str, ...]

    @property
    def n_viruses(self) -> int:
        return len({locus.virus for locus in self.positives})

    def sensitivity_at(self, threshold: float, min_tools: int = 2) -> float:
        """Share of confirmed G4s the rule would admit."""
        if not self.positives:
            return 0.0
        found = sum(
            1 for p in self.positives if p.magnitude >= threshold and p.n_tools >= min_tools
        )
        return found / len(self.positives)

    def explain(self) -> str:
        if not self.usable:
            return (
                "This calibration set cannot support an operating point: "
                + "; ".join(self.blocking_reasons)
                + ". The numbers below describe the set, not a recommended threshold."
            )
        return (
            f"{len(self.positives)} confirmed G4s across {self.n_viruses} viruses, "
            f"{len(self.negatives)} matched negatives."
        )


def load_confirmed_set(path: str | Path) -> list[dict]:
    """Read the curated ground-truth TSV."""
    with Path(path).open(newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row.get("locus_id")]


def score_region(sequence: str, start: int, end: int, *, window: int = 25,
                 threshold: float = 1.2) -> tuple[float | None, int]:
    """Best |G4Hunter| in a 1-based inclusive span, and the tool count.

    Returns ``(None, 0)`` when the predictor finds nothing there at all --
    a miss, which is different from a weak hit and must not be scored as
    one.
    """
    span = sequence[start - 1 : end]
    hits = g4hunter.predict(span, window, threshold)
    if not hits:
        return None, 0
    best = max(hits, key=lambda h: abs(h.score))
    motifs = pattern_motif.predict(span)
    overlapping = [
        m for m in motifs if m.start < best.end and best.start < m.end
    ]
    return best.score, 2 if overlapping else 1


def roc_curve(positives: list[ScoredLocus], negatives: list[ScoredLocus],
              min_tools: int = 1) -> list[dict]:
    """Sensitivity and specificity across candidate |G4Hunter| cut-offs.

    A locus the predictor missed entirely counts against sensitivity at
    every threshold, which is correct: no cut-off recovers a hit that was
    never made.
    """
    cuts = sorted({round(x * 0.05, 2) for x in range(0, 61)})
    out = []
    for cut in cuts:
        tp = sum(1 for p in positives if p.magnitude >= cut and p.n_tools >= min_tools)
        fp = sum(1 for n in negatives if n.magnitude >= cut and n.n_tools >= min_tools)
        sens = tp / len(positives) if positives else 0.0
        spec = 1 - (fp / len(negatives)) if negatives else 0.0
        out.append({
            "threshold": cut,
            "min_tools": min_tools,
            "sensitivity": sens,
            "specificity": spec,
            "youden_j": sens + spec - 1,
            "n_true_positive": tp,
            "n_false_positive": fp,
        })
    return out


def build_report(positives: list[ScoredLocus], negatives: list[ScoredLocus],
                 min_tools: int = 1) -> CalibrationReport:
    reasons: list[str] = []
    if len(positives) < MIN_POSITIVES_FOR_CALIBRATION:
        reasons.append(
            f"{len(positives)} confirmed positives, {MIN_POSITIVES_FOR_CALIBRATION} needed"
        )
    viruses = len({p.virus for p in positives})
    if viruses < MIN_VIRUSES_FOR_CALIBRATION:
        reasons.append(f"positives from {viruses} virus(es), {MIN_VIRUSES_FOR_CALIBRATION} needed")
    derived = sum(1 for p in positives if p.provenance == PROVENANCE_DERIVED)
    if derived:
        reasons.append(
            f"{derived} of {len(positives)} positives have DERIVED coordinates, located by the "
            "predictor being calibrated -- mildly circular until upgraded to stated"
        )
    if not negatives:
        reasons.append("no matched negatives were supplied, so specificity is unmeasured")
    return CalibrationReport(
        positives=positives,
        negatives=negatives,
        curve=roc_curve(positives, negatives, min_tools=min_tools),
        usable=not reasons,
        blocking_reasons=tuple(reasons),
    )
