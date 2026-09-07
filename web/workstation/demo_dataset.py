"""SYNTHETIC DEMONSTRATION DATASET -- NOT REAL SCIENCE.

Why this exists: several states of the interface cannot be reached with
the real corpus. FMDV's D.H1 gate is BLOCKED_INSUFFICIENT_DATA and,
correctly, will not open -- so the scored model view, the completed
pipeline, the SUPPORTED verdict and the significant-finding styling can
never be seen, reviewed, or shown to a stakeholder. This module supplies
a fabricated corpus so those screens can be exercised.

Every value is invented here and nothing is read from ``data/``:

* ``DEMO-A`` ... ``DEMO-G`` are not serotypes.
* ``DEMO-REF-000`` is not an accession.
* The genome is a seeded pseudo-random sequence, not a virus.
* Countries are placeholders ("Region 1" ...), never real nations.

Two guarantees hold, and there are tests for both:

1. ``identity.synthetic`` is ``True`` and ``display_name`` carries the
   word DEMO, so every surface can mark it. The UI refuses to render it
   without a standing watermark.
2. Nothing here can reach the real testing ledger. This module has no
   write path at all -- it returns a dict.

The real FMDV verdict is unaffected and still stands: INSUFFICIENT_DATA,
because serotype C has 8 genomes against a floor of 20.
"""

from __future__ import annotations

import math
import random

from .dataset import layout_tree, parse_newick

SEED = 20260901
GENOME_LEN = 9_450
CDS = [1_120, 9_180]
LINEAGES = ["DEMO-A", "DEMO-B", "DEMO-C", "DEMO-D", "DEMO-E", "DEMO-F", "DEMO-G"]
REGIONS = [f"Region {i}" for i in range(1, 15)]
HOSTS = ["host-alpha", "host-beta", "host-gamma"]


def _newick(rng: random.Random, labels: list[str]) -> str:
    """A random rooted binary topology with branch lengths.

    Built as a string and handed to the same parser the real tree uses,
    so the demo exercises the production layout code rather than a
    parallel implementation that could drift from it.
    """
    nodes = [f"{name}:{rng.uniform(0.002, 0.05):.6f}" for name in labels]
    rng.shuffle(nodes)
    while len(nodes) > 1:
        i = rng.randrange(len(nodes) - 1)
        a = nodes.pop(i)
        b = nodes.pop(rng.randrange(len(nodes)))
        nodes.append(f"({a},{b}):{rng.uniform(0.001, 0.03):.6f}")
    return nodes[0].rsplit(":", 1)[0] + ";"


def build_demo_dataset(n_samples: int = 420) -> dict:
    rng = random.Random(SEED)

    # ── corpus ───────────────────────────────────────────────────
    weights = [0.26, 0.20, 0.16, 0.13, 0.11, 0.08, 0.06]
    samples = []
    for i in range(n_samples):
        lineage = rng.choices(LINEAGES, weights=weights)[0]
        samples.append(
            {
                "a": f"DEMO{i:04d}.1",
                "l": lineage,
                "c": rng.choices(REGIONS, weights=[14, 12, 11, 9, 8, 7, 6, 6, 5, 5, 4, 4, 3, 3])[0],
                "y": rng.choices(
                    range(2009, 2027), weights=[1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 11, 9, 6]
                )[0],
                "h": rng.choice(HOSTS),
            }
        )

    lineages: dict[str, int] = {}
    countries: dict[str, int] = {}
    years: dict[int, int] = {}
    for s in samples:
        lineages[s["l"]] = lineages.get(s["l"], 0) + 1
        countries[s["c"]] = countries.get(s["c"], 0) + 1
        years[s["y"]] = years.get(s["y"], 0) + 1

    # ── tree ─────────────────────────────────────────────────────
    tree = layout_tree(parse_newick(_newick(rng, [s["a"] for s in samples])))
    index = {s["a"]: i for i, s in enumerate(samples)}
    for node in tree["nodes"]:
        node["s"] = index.get(node["n"], -1) if node["n"] else -1

    # ── genome + loci ────────────────────────────────────────────
    genome = "".join(rng.choices("ACGT", weights=[24, 27, 26, 23], k=GENOME_LEN))
    counts = {b: genome.count(b) for b in "ACGT"}
    composition = {
        "total": GENOME_LEN,
        "counts": counts,
        "percent": {b: round(counts[b] / GENOME_LEN * 100, 2) for b in "ACGT"},
        "gc": round((counts["G"] + counts["C"]) / GENOME_LEN * 100, 2),
        "other": 0,
    }

    spec = [
        ("DEMO-G4-001", 402, 449, "-", "5' UTR", 3, -1.81, 78.4, "EC", "known_functional_region"),
        ("DEMO-G4-002", 1_042, 1_088, "+", "5' UTR", 3, 1.94, 71.2, "EC", "known_functional_region"),
        ("DEMO-G4-003", 2_910, 2_948, "+", "CDS (polyprotein)", 2, 1.42, 55.8, "BC", "unannotated"),
        ("DEMO-G4-004", 4_688, 4_722, "-", "CDS (polyprotein)", 2, -1.36, 61.0, "BC", "unannotated"),
        ("DEMO-G4-005", 6_401, 6_433, "+", "CDS (polyprotein)", 1, 0.88, 29.7, "WC", "unannotated"),
        ("DEMO-G4-006", 9_221, 9_262, "-", "3' UTR", 3, -1.66, 66.5, "EC", "known_functional_region"),
    ]
    loci = []
    for atlas_id, start, end, strand, feature, tools, score, cons, conf, ctx in spec:
        seq = genome[start - 1 : end]
        n = len(seq)
        loci.append(
            {
                "id": atlas_id,
                "start": start,
                "end": end,
                "strand": strand,
                "feature": feature,
                "seq": seq,
                "g4hunter": score,
                "tools": tools,
                "topology": "parallel",
                "g4_type": "RNA_G4",
                "gc_flank": round(rng.uniform(52, 71), 2),
                "conservation": cons,
                "confidence": conf,
                "context": ctx,
                "composition": {b: round(seq.count(b) / n * 100, 1) for b in "ACGT"},
                "gc": round((seq.count("G") + seq.count("C")) / n * 100, 1),
            }
        )

    supported = [locus["id"] for locus in loci if locus["confidence"] in {"EC", "BC"}]

    # ── floor: everything clears, which is why the gate can open ──
    smallest = min(lineages.values())
    floor = {
        "n_sequences_in_window": {"value": len(samples), "floor": 30, "unit": "seq"},
        "min_sequences_per_lineage": {
            "value": smallest,
            "floor": 20,
            "unit": "seq",
            "which": min(lineages, key=lineages.get),
        },
        "n_timepoints": {"value": len(years), "floor": 3, "unit": "yr"},
        "metadata_completeness": {"value": 0.9881, "floor": 0.90, "unit": "frac"},
        "alignment_qc_pass_fraction": {"value": 0.9143, "floor": 0.50, "unit": "frac"},
    }

    real_countries = len(countries)
    entropy = -sum((v / len(samples)) * math.log(v / len(samples)) for v in countries.values())

    return {
        "identity": {
            "pathogen": "DEMO",
            "display_name": "DEMO — Synthetic Pathogen (fabricated data)",
            "synthetic": True,
            "warning": "SYNTHETIC DEMONSTRATION DATA — NOT A RESULT. Every value is fabricated. "
            "Nothing here describes any real organism and none of it may be reported or cited.",
            "genome_type": "ssRNA_positive",
            "reference": "DEMO-REF-000",
            "genome_length": GENOME_LEN,
            "cds": CDS,
            "atlas_version": "0.0-demo",
            "lineage_field": "lineage",
            "n_samples": len(samples),
            "n_raw": int(len(samples) / 0.9143),
            "n_lineages": len(lineages),
            "n_countries": real_countries,
            "period": [min(years), max(years)],
            "operational_mode": True,
        },
        "samples": samples,
        "lineages": dict(sorted(lineages.items(), key=lambda kv: -kv[1])),
        "countries": dict(sorted(countries.items(), key=lambda kv: -kv[1])),
        "years": dict(sorted(years.items())),
        "tree": tree,
        "ancestral": {},
        "loci": loci,
        "composition": composition,
        "gate": {
            "permission": "PERMITTED",
            "permitted": True,
            "operational_mode": True,
            "failing_checks": [],
            "supported_loci": supported,
            "latest_run": "2026-09-01T09:00:00+00:00",
            "n_ledger_rows": len(loci),
            "explanation": "DEMO: the synthetic corpus clears every Appendix C check, so the D.H1 gate "
            "returns SUPPORTED and the scored views become reachable. This is a fabricated "
            "demonstration of the open-gate interface, not a finding about any organism.",
        },
        "observations": [
            {
                "severity": "block",
                "domain": "provenance",
                "title": "This dataset is fabricated",
                "detail": "Every genome, lineage, coordinate and score on this screen was generated by "
                "web/workstation/demo_dataset.py. It exists so the open-gate interface can be "
                "reviewed. It must never be exported as a result, screenshotted into a report, "
                "or compared against the real FMDV Atlas.",
                "rule": "identity.synthetic == True",
                "focus": {"type": "gate", "value": None},
            },
            {
                "severity": "caution",
                "domain": "structure",
                "title": f"{sum(1 for locus in loci if locus['tools'] <= 1)} locus rests on a single algorithm",
                "detail": "Even in the demo, single-tool support is reported as weak. Inter-algorithm "
                "discordance for G4 prediction runs 30–60%.",
                "rule": "concordant_tool_count <= 1",
                "focus": {"type": "locus", "value": "DEMO-G4-005"},
            },
            {
                "severity": "note",
                "domain": "evolution",
                "title": "Conservation spans 29.7–78.4% across loci",
                "detail": "DEMO-G4-001 is the most conserved candidate; DEMO-G4-005 the least.",
                "rule": "range(conservation_pct_phylo)",
                "focus": {"type": "locus", "value": "DEMO-G4-001"},
            },
            {
                "severity": "note",
                "domain": "geography",
                "title": f"{real_countries} regions, Shannon H = {entropy:.2f}",
                "detail": "Placeholder regions, not real nations.",
                "rule": "Shannon entropy over region counts",
                "focus": {"type": "country", "value": max(countries, key=countries.get)},
            },
        ],
        "floor": floor,
    }


# ── synthetic equivalents of the derived views ──────────────────────
def demo_tracks(window: int = 40) -> dict:
    """Fabricated genome tracks matching the real endpoint's shape."""
    rng = random.Random(SEED + 1)
    d = build_demo_dataset()
    n_win = GENOME_LEN // window + 1
    starts = [i * window + 1 for i in range(n_win)]

    diversity, gc, gaps = [], [], []
    for start in starts:
        base = 0.06 + 0.05 * math.sin(start / 900.0)
        inside = any(locus["start"] - 60 <= start <= locus["end"] + 60 for locus in d["loci"])
        # Loci are drawn less variable, which is the pattern the real
        # hypothesis predicts. In the demo it is asserted, not measured.
        diversity.append(round(max(0.0, base * (0.35 if inside else 1.0) + rng.gauss(0, 0.012)), 5))
        gc.append(round(min(0.95, max(0.20, 0.52 + (0.16 if inside else 0.0) + rng.gauss(0, 0.03))), 5))
        gaps.append(round(max(0.0, rng.gauss(0.02, 0.015)), 5))

    background = sum(diversity) / len(diversity)
    loci = []
    for locus in d["loci"]:
        lo, hi = (locus["start"] - 1) // window, (locus["end"] - 1) // window + 1
        seg = diversity[lo:hi] or [background]
        value = sum(seg) / len(seg)
        loci.append(
            {
                "id": locus["id"],
                "diversity": round(value, 5),
                "background": round(background, 5),
                "ratio": round(value / background, 3),
            }
        )

    return {
        "window": window,
        "n_windows": n_win,
        "genome_length": GENOME_LEN,
        "n_sequences": d["identity"]["n_samples"],
        "starts": starts,
        "diversity": diversity,
        "gc": gc,
        "gaps": gaps,
        "loci": loci,
        "synthetic": True,
    }


def demo_ordination() -> dict:
    """Fabricated ordination: lineage clusters with scatter."""
    rng = random.Random(SEED + 2)
    d = build_demo_dataset()
    centres = {
        L: (math.cos(i * 2 * math.pi / len(LINEAGES)) * 0.9, math.sin(i * 2 * math.pi / len(LINEAGES)) * 0.9)
        for i, L in enumerate(LINEAGES)
    }
    xs, ys, ids = [], [], []
    for s in d["samples"]:
        cx, cy = centres[s["l"]]
        ids.append(s["a"])
        xs.append(round(cx + rng.gauss(0, 0.18), 5))
        ys.append(round(cy + rng.gauss(0, 0.18), 5))
    return {"ids": ids, "x": xs, "y": ys, "explained": [38.4, 22.7], "n": len(ids), "synthetic": True}


# ── a corpus with per-genome G4 states, for the Stage 5 driver ──────
def demo_samples(n: int = 2600):
    """Synthetic ``surveillance_metrics.Sample`` objects with tip states.

    The generator plants a real signal so the downstream chain has
    something to find: from 2019 a rising fraction of one lineage's
    genomes carry a disrupted G4, and that lineage subsequently gains
    frequency. Nothing here is a claim about biology -- it is a test
    fixture whose only job is to exercise every branch of the driver on
    data whose right answer is known by construction.
    """
    from g4watch.metrics.surveillance_metrics import Sample

    rng = random.Random(SEED + 7)
    d = build_demo_dataset()
    loci = [locus["id"] for locus in d["loci"]]
    # A 46-year span: CUSUM and EWMA need >= 20 baseline windows to
    # calibrate a control limit, and calibrating on fewer is refused.
    years = list(range(1981, 2027))
    focal = "DEMO-B"

    out = []
    for i in range(n):
        # the focal lineage grows from ~12% to ~34% across the span
        year = rng.choices(years, weights=[4 + i * 0.6 for i in range(len(years))])[0]
        t = (year - years[0]) / (years[-1] - years[0])
        weights = [0.24 - 0.06 * t, 0.12 + 0.22 * t, 0.15, 0.13, 0.13, 0.13 - 0.08 * t, 0.10 - 0.08 * t]
        lineage = rng.choices(LINEAGES, weights=[max(w, 0.01) for w in weights])[0]

        # disruption rises with time in the focal lineage only
        p_disrupt = (0.05 + 0.45 * max(0.0, t - 0.45) / 0.55) if lineage == focal else 0.06
        p_gain = 0.03
        states = {}
        for locus in loci:
            r = rng.random()
            states[locus] = "disrupted" if r < p_disrupt else "gained" if r < p_disrupt + p_gain else "present"

        total_mut = max(1, int(rng.gauss(60 + 26 * t, 12)))
        g4_mut = max(0, int(total_mut * rng.gauss(0.09 + 0.05 * p_disrupt, 0.02)))
        out.append(Sample(
            accession=f"DEMO{i:04d}.1", lineage=lineage,
            country=rng.choices(REGIONS, weights=[14, 12, 11, 9, 8, 7, 6, 6, 5, 5, 4, 4, 3, 3])[0],
            year=year, states=states, g4_mutations=g4_mut, total_mutations=total_mut,
        ))
    return out
