"""The pathogen-specific report card.

One document per pathogen carrying every section the framework asks for:
identity, dataset, quality, genomic metrics, the four hypothesis verdicts,
early-warning indicators, the score, confidence, validation status,
limitations and a final interpretation.

Three rules govern what may appear here, and they are enforced in code
rather than left to whoever assembles a report:

1. **A section that has no data says so.** Every section carries a status
   of ``reported``, ``blocked`` or ``unavailable`` with a reason. A blank
   section and an absent section mean different things and are never
   conflated.

2. **Scored sections are gated.** The surveillance score, warning level
   and early-warning indicators appear only when the D.H1 gate permits
   scoring for this pathogen. Otherwise they are marked ``blocked`` with
   the gate's own explanation, so a reader learns why rather than seeing
   an empty box.

3. **Synthetic input is stamped, not silently carried.** A card built
   from a demonstration corpus is marked non-authoritative at the top
   level, and every consumer can see it without inspecting the data.

The card is a plain dictionary. Rendering it -- to the dashboard, to HTML,
to a Word document -- is somebody else's job, which keeps the content
decisions in one place.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

REPORTED = "reported"
BLOCKED = "blocked"
UNAVAILABLE = "unavailable"


@dataclass
class Section:
    key: str
    title: str
    status: str
    reason: str = ""
    data: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"key": self.key, "title": self.title, "status": self.status, "reason": self.reason, "data": self.data}


def _section(key: str, title: str, data: dict | None, reason: str = "") -> Section:
    if data:
        return Section(key, title, REPORTED, reason, data)
    return Section(key, title, UNAVAILABLE, reason or "no data available for this section")


def build_report_card(
    dataset: dict,
    *,
    stage5: dict | None = None,
    dh3: dict | None = None,
    tracks: dict | None = None,
    spectrum: dict | None = None,
) -> dict:
    """Assemble the card from whatever the pipeline has produced.

    Every argument except ``dataset`` is optional. A missing input yields
    an ``unavailable`` section naming what would have produced it, which
    is more useful to a reader than its silent absence.
    """
    # D.H3's result is already part of a Stage 5 payload. Requiring it to
    # be passed a second time meant every command-line card reported
    # "clade growth classification has not been run" directly after D.H3
    # had run, with its verdict sitting unread in ``stage5["dh3"]``. An
    # explicit argument still wins, so callers that compute it separately
    # (the web app does) are unaffected.
    if dh3 is None and stage5:
        dh3 = stage5.get("dh3") or None

    identity = dataset.get("identity", {})
    gate = dataset.get("gate", {})
    floor = dataset.get("floor", {})
    loci = dataset.get("loci", [])
    synthetic = bool(identity.get("synthetic"))
    permitted = bool(gate.get("permitted"))

    sections: list[Section] = []

    # ── 1 identity ──────────────────────────────────────────────
    sections.append(
        _section(
            "identity",
            "Pathogen identity",
            {
                "pathogen": identity.get("pathogen"),
                "display_name": identity.get("display_name"),
                "genome_type": identity.get("genome_type"),
                "reference": identity.get("reference"),
                "genome_length": identity.get("genome_length"),
                "atlas_version": identity.get("atlas_version"),
                "lineage_field": identity.get("lineage_field"),
            },
        )
    )

    # ── 2 dataset summary ───────────────────────────────────────
    sections.append(
        _section(
            "dataset",
            "Dataset summary",
            {
                "genomes_analysed": identity.get("n_samples"),
                "genomes_acquired": identity.get("n_raw"),
                "lineages": identity.get("n_lineages"),
                "countries": identity.get("n_countries"),
                "period": identity.get("period"),
                "lineage_counts": dataset.get("lineages", {}),
                "top_countries": dict(list(dataset.get("countries", {}).items())[:10]),
                "genomes_per_year": dataset.get("years", {}),
            },
        )
    )

    # ── 3 data quality ──────────────────────────────────────────
    checks = []
    for name, row in (floor or {}).items():
        if row.get("value") is None:
            continue
        checks.append(
            {
                "check": name,
                "value": row["value"],
                "floor": row["floor"],
                "unit": row.get("unit"),
                "passed": row["value"] >= row["floor"],
                "which": row.get("which"),
            }
        )
    failed = [c["check"] for c in checks if not c["passed"]]
    sections.append(
        Section(
            "quality",
            "Data quality and minimum-data floor",
            REPORTED if checks else UNAVAILABLE,
            "" if checks else "the Appendix C floor has not been evaluated",
            {
                "checks": checks,
                "n_failing": len(failed),
                "failing": failed,
                "note": (
                    "Clearing the floor permits a test to run. It does not imply "
                    "the corpus is adequately powered, which is assessed separately."
                ),
            },
        )
    )

    # ── 4 genomic metrics ───────────────────────────────────────
    genomic: dict[str, Any] = {
        "n_loci": len(loci),
        "loci": [
            {
                k: locus.get(k)
                for k in (
                    "id",
                    "start",
                    "end",
                    "strand",
                    "feature",
                    "g4hunter",
                    "tools",
                    "conservation",
                    "confidence",
                    "context",
                    "gc",
                )
            }
            for locus in loci
        ],
        "composition": dataset.get("composition", {}),
    }
    if tracks:
        genomic["locus_vs_background"] = tracks.get("loci", [])
        genomic["window"] = tracks.get("window")
    if spectrum:
        genomic["mutation_spectrum"] = {
            "ti_tv": spectrum.get("ti_tv"),
            "transitions": spectrum.get("transitions"),
            "transversions": spectrum.get("transversions"),
        }
    sections.append(_section("genomic", "Key genomic metrics", genomic))

    # ── 5 hypothesis verdicts ───────────────────────────────────
    hyps = {
        "D.H1": {
            "claim": "G4 regions differ in conservation from matched non-G4 controls",
            "verdict": gate.get("permission"),
            "status": REPORTED if gate.get("permission") else UNAVAILABLE,
            "detail": gate.get("explanation"),
            "supported_loci": gate.get("supported_loci", []),
            "ledger_rows": gate.get("n_ledger_rows"),
        },
    }
    mc = (stage5 or {}).get("model_comparison") or {}
    for tag, claim, evidence in (
        (
            "D.H2",
            "G4 changes associate with subsequent lineage expansion",
            "the nested M4-vs-M2 likelihood ratio test",
        ),
        (
            "D.H4",
            "G4 features add predictive value over conventional metrics",
            "the M4-vs-M2 LRT together with held-out AUC",
        ),
    ):
        if mc.get("status") in {"failed", "not_attempted"} or not mc:
            hyps[tag] = {
                "claim": claim,
                "verdict": None,
                "status": BLOCKED,
                "detail": mc.get("detail", "the M1-M4 comparison has not been run"),
            }
            continue

        # "EVALUATED" reported only that a test had happened. The finding
        # itself -- ComparisonResult.g4_adds_value(), which the comparison
        # module documents as the D.H2/D.H4 question -- was computed and
        # then left out of the card, and both hypotheses were given the
        # same generic method description regardless of the answer.
        adds = mc.get("g4_adds_value")
        if adds is None:
            verdict, detail = "EVALUATED", mc.get("note")
        else:
            verdict = "SUPPORTED" if adds else "NOT_SUPPORTED"
            detail = f"{mc.get('verdict')}, on {evidence}."
        hyps[tag] = {
            "claim": claim,
            "verdict": verdict,
            "status": REPORTED,
            "detail": detail,
            "method": mc.get("note"),
            "best_by_aic": mc.get("best_by_aic"),
            "best_by_holdout_auc": mc.get("best_by_holdout_auc"),
            "likelihood_ratio_tests": mc.get("likelihood_ratio_tests"),
            "n_train": mc.get("n_train"),
            "n_holdout": mc.get("n_holdout"),
        }
    if dh3:
        hyps["D.H3"] = {
            "claim": "G4 transitions cluster in expanding clades",
            "verdict": dh3.get("verdict"),
            "status": REPORTED if dh3.get("verdict") != "INSUFFICIENT_DATA" else BLOCKED,
            "detail": dh3.get("explanation"),
            "p_value": dh3.get("p_value"),
            "underpowered": dh3.get("underpowered"),
        }
    else:
        hyps["D.H3"] = {
            "claim": "G4 transitions cluster in expanding clades",
            "verdict": None,
            "status": BLOCKED,
            "detail": "clade growth classification has not been run for this pathogen",
        }
    sections.append(Section("hypotheses", "Research hypothesis verdicts", REPORTED, "", {"hypotheses": hyps}))

    # ── 6 early warning + score (gated) ─────────────────────────
    if not permitted:
        sections.append(
            Section(
                "early_warning",
                "Early-warning indicators",
                BLOCKED,
                gate.get("explanation", "the D.H1 gate does not permit scoring"),
                {"gate": gate.get("permission")},
            )
        )
        sections.append(
            Section(
                "score",
                "Surveillance score",
                BLOCKED,
                "No score, alert level or warning is computed while the gate is "
                "closed. This is a reported outcome, not a missing one.",
                {"gate": gate.get("permission")},
            )
        )
    elif stage5:
        series = stage5.get("score_series", [])
        detection = stage5.get("detection", {})
        sections.append(
            _section(
                "early_warning",
                "Early-warning indicators",
                {
                    "cusum": detection.get("cusum"),
                    "ewma": detection.get("ewma"),
                    "warning": stage5.get("warning"),
                },
                "detection charts calibrated on the first half of the series",
            )
        )
        sections.append(
            _section(
                "score",
                "Surveillance score",
                {
                    "n_windows_scored": len(series),
                    "latest": series[-1] if series else None,
                    "series": series,
                    "weights": stage5.get("weights"),
                },
            )
        )
    else:
        sections.append(
            Section(
                "early_warning",
                "Early-warning indicators",
                UNAVAILABLE,
                "the gate permits scoring but Stage 5 has not been run",
            )
        )
        sections.append(
            Section("score", "Surveillance score", UNAVAILABLE, "the gate permits scoring but Stage 5 has not been run")
        )

    # ── 7 emerging signals ──────────────────────────────────────
    observations = dataset.get("observations", [])
    sections.append(
        _section(
            "signals",
            "Emerging signals and observations",
            {
                "n_blocking": sum(1 for o in observations if o.get("severity") == "block"),
                "n_caution": sum(1 for o in observations if o.get("severity") == "caution"),
                "observations": observations,
            },
            "each observation states the rule that produced it",
        )
    )

    # ── 8 confidence ────────────────────────────────────────────
    tiers: dict[str, int] = {}
    for locus in loci:
        tiers[locus.get("confidence", "?")] = tiers.get(locus.get("confidence", "?"), 0) + 1
    single_tool = sum(1 for locus in loci if (locus.get("tools") or 0) <= 1)
    sections.append(
        _section(
            "confidence",
            "Confidence assessment",
            {
                "structural_confidence": tiers,
                "single_algorithm_loci": single_tool,
                "evidence_class": "computational only" if single_tool == len(loci) and loci else "mixed",
                "note": (
                    "Published inter-algorithm discordance for G4 prediction runs 30-60%. A locus "
                    "supported by one tool is a candidate motif, not a demonstrated structure."
                ),
            },
        )
    )

    # ── 9 validation status ─────────────────────────────────────
    validation = {
        "minimum_data_floor": "evaluated" if checks else "not evaluated",
        "gc_confound_adjustment": "applied at D.H1",
        "multiple_testing": "Benjamini-Hochberg across loci",
        "recombination_screen": "mandatory; run before D.H1",
        "holdout_partition": (stage5 or {}).get("partition", {}).get("method", "not run"),
        "external_validation": "not performed",
        "calibration": "not performed",
    }
    sections.append(_section("validation", "Validation status", validation))

    # ── 10 limitations ──────────────────────────────────────────
    limitations = []
    if failed:
        limitations.append(f"Minimum-data floor not met: {', '.join(failed)}.")
    if single_tool and single_tool == len(loci):
        limitations.append("Every Atlas locus rests on a single prediction algorithm.")
    if not permitted:
        limitations.append("Scoring is not authorised, so no surveillance claim can be made.")
    if synthetic:
        limitations.append("This card was built from a synthetic demonstration corpus and describes no real organism.")
    if dh3 and dh3.get("underpowered"):
        limitations.append("The D.H3 comparison is underpowered; see its group sizes.")
    limitations.append(
        "Sequencing in this field is outbreak-driven rather than systematic, so "
        "temporal and geographic comparisons are not sampling-balanced."
    )
    sections.append(_section("limitations", "Limitations", {"items": limitations}))

    # ── 11 interpretation ───────────────────────────────────────
    if synthetic:
        headline = (
            "Synthetic demonstration corpus. Nothing in this card is a finding about any "
            "organism; it exists so the open-gate interface can be reviewed."
        )
    elif not permitted:
        headline = f"No surveillance claim is made for {identity.get('pathogen')}. {gate.get('explanation', '')}"
    else:
        headline = (
            f"{identity.get('pathogen')} cleared the minimum-data floor and the D.H1 gate. "
            "Scored outputs in this card are provisional until externally validated."
        )
    sections.append(
        _section(
            "interpretation",
            "Final scientific interpretation",
            {
                "headline": headline,
                "research_use_only": True,
                "statement": (
                    "G4-WATCH is a research framework, not a validated diagnostic or an "
                    "operational alert system. No output here should drive a control decision "
                    "on its own."
                ),
            },
        )
    )

    return {
        "schema": "g4watch/report-card@1",
        "pathogen": identity.get("pathogen"),
        "display_name": identity.get("display_name"),
        "authoritative": not synthetic,
        "synthetic": synthetic,
        "scoring_permitted": permitted,
        "generated_from": {
            "dataset": True,
            "stage5": stage5 is not None,
            "dh3": dh3 is not None,
            "tracks": tracks is not None,
            "spectrum": spectrum is not None,
        },
        "sections": [s.as_dict() for s in sections],
        "section_status": {s.key: s.status for s in sections},
    }
