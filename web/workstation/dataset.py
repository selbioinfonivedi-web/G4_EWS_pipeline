"""Assemble the linked-view payload the workstation canvas renders.

Every field here is read from a real pipeline artifact -- the corpus
metadata table, the IQ-TREE rooted Newick, the ancestral-state
reconstruction, the G4 Reference Atlas and the testing ledger. Nothing is
synthesised. Where an artifact is absent the corresponding view is
reported as unavailable rather than filled with a plausible substitute,
because a genomics workstation that invents a tree is worse than one that
admits it has none.

The payload is deliberately one document: every view in the workstation
is a projection of the same sample list, and cross-highlighting only
works because they share `accession` as the join key.
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from g4watch.config import PathogenConfig, load_config
from g4watch.gating import evaluate_gate, read_ledger
from g4watch.io.fasta import read_fasta
from g4watch.pipeline.stage45_dh1 import _extract_year

REPO_ROOT = Path(__file__).resolve().parents[2]

# Appendix C floors, mirrored here for display only. The gate itself is
# always computed by g4watch.validation, never by this module.
FLOORS = {
    "n_sequences_in_window": 30,
    "min_sequences_per_lineage": 20,
    "n_timepoints": 3,
    "metadata_completeness": 0.90,
    "alignment_qc_pass_fraction": 0.50,
}


# ── newick ──────────────────────────────────────────────────────────
@dataclass
class TreeNode:
    idx: int
    parent: int | None
    name: str
    branch: float
    children: list[int]


def parse_newick(text: str) -> list[TreeNode]:
    """A small Newick reader: internal labels and branch lengths, no comments.

    Written here rather than pulled in as a dependency because the tree
    only has to be laid out for display; nothing downstream of this
    module does phylogenetics with the result.
    """
    text = text.strip().rstrip(";")
    nodes: list[TreeNode] = []
    stack: list[int] = []
    token = ""

    def new_node(parent: int | None) -> int:
        nodes.append(TreeNode(idx=len(nodes), parent=parent, name="", branch=0.0, children=[]))
        if parent is not None:
            nodes[parent].children.append(len(nodes) - 1)
        return len(nodes) - 1

    def finish(idx: int, raw: str) -> None:
        raw = raw.strip()
        if not raw:
            return
        if ":" in raw:
            label, _, length = raw.rpartition(":")
            nodes[idx].name = label
            try:
                nodes[idx].branch = float(length)
            except ValueError:
                nodes[idx].branch = 0.0
        else:
            nodes[idx].name = raw

    root = new_node(None)
    stack.append(root)
    current: int | None = None

    for char in text:
        if char == "(":
            parent = stack[-1] if current is None else current
            child = new_node(parent)
            stack.append(child)
            current = None
            token = ""
        elif char == ",":
            if current is not None:
                finish(current, token)
                current = None
            else:
                current = new_node(stack[-1])
                finish(current, token)
                current = None
            token = ""
        elif char == ")":
            if current is not None:
                finish(current, token)
                current = None
            else:
                leaf = new_node(stack[-1])
                finish(leaf, token)
            token = ""
            current = stack.pop()
        else:
            token += char

    if current is not None:
        finish(current, token)
    return nodes


def layout_tree(nodes: list[TreeNode]) -> dict:
    """Rectangular cladogram layout: x = divergence from root, y = leaf order."""
    root = 0
    depth: dict[int, float] = {}
    order: list[int] = []

    # Iterative post-order so a deep tree cannot blow the recursion limit.
    stack: list[tuple[int, bool]] = [(root, False)]
    while stack:
        idx, expanded = stack.pop()
        if expanded:
            order.append(idx)
            continue
        parent = nodes[idx].parent
        depth[idx] = (depth[parent] if parent is not None else 0.0) + nodes[idx].branch
        stack.append((idx, True))
        for child in reversed(nodes[idx].children):
            stack.append((child, False))

    leaf_y: dict[int, float] = {}
    counter = 0
    for idx in order:
        if not nodes[idx].children:
            leaf_y[idx] = float(counter)
            counter += 1

    y: dict[int, float] = dict(leaf_y)
    for idx in order:
        kids = nodes[idx].children
        if kids:
            y[idx] = sum(y[k] for k in kids) / len(kids)

    max_depth = max(depth.values()) or 1.0
    return {
        "n_leaves": counter,
        "max_depth": max_depth,
        "nodes": [
            {
                "i": n.idx,
                "p": n.parent,
                "x": round(depth[n.idx] / max_depth, 6),
                "y": round(y[n.idx] / max(counter - 1, 1), 6),
                "n": n.name if not n.children else "",
            }
            for n in nodes
        ],
    }


# ── helpers ─────────────────────────────────────────────────────────
def _year(raw: str) -> int | None:
    """Year of collection, using the pipeline's own parser.

    Deliberately delegates to ``stage45_dh1._extract_year`` rather than
    re-implementing it. GenBank dates in this corpus look like
    ``18-Jul-2024/21-Jul-2024``, and a second, subtly different parser
    here would make the workstation's timepoint count disagree with the
    one the Appendix C floor is actually evaluated on.
    """
    text = _extract_year(raw or "")
    if not text:
        return None
    value = int(text)
    return value if 1900 <= value <= 2100 else None


def _resolve(config: PathogenConfig, *parts: str) -> Path:
    return REPO_ROOT.joinpath(*parts)


# ── payload ─────────────────────────────────────────────────────────
@lru_cache(maxsize=4)
def build_dataset(pathogen: str) -> dict:
    config = load_config(pathogen)

    aligned_path = _resolve(
        config,
        "data/reference_genomes",
        pathogen,
        "corpus/aligned",
        f"{pathogen}_qc_passed_aligned_to_ref.fasta",
    )
    tree_path = _resolve(
        config, "data/reference_genomes", pathogen, "corpus/phylogenetics", f"{pathogen}_iqtree_rooted.nwk"
    )
    ancestral_path = _resolve(
        config, "data/reference_genomes", pathogen, "corpus/phylogenetics", "ancestral_serotype_reconstruction.tsv"
    )

    aligned_ids: set[str] = set(read_fasta(aligned_path)) if aligned_path.is_file() else set()

    # ── samples ──────────────────────────────────────────────────
    samples: list[dict] = []
    metadata_path = Path(config.corpus_metadata_tsv) if config.corpus_metadata_tsv else None
    if metadata_path and metadata_path.is_file():
        from g4watch.qc.metadata_normalization import LineageVocabulary

        vocabulary = LineageVocabulary.from_config(config)
        with open(metadata_path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                accession = row["accession"]
                if aligned_ids and accession not in aligned_ids:
                    continue
                fields = [row.get(config.lineage_field, "")]
                fields += [row.get(name, "") for name in config.lineage_fallback_fields]
                lineage = vocabulary.resolve(*fields) or (row.get(config.lineage_field, "") or "").strip()
                samples.append(
                    {
                        "a": accession,
                        "l": lineage or "—",
                        "c": (row.get("country", "") or "").split(":")[0].strip() or "—",
                        "y": _year(row.get("collection_date", "")),
                        "h": (row.get("host", "") or "").strip() or "—",
                    }
                )

    index = {s["a"]: i for i, s in enumerate(samples)}

    # ── lineage / geography / time ───────────────────────────────
    lineages: dict[str, int] = {}
    countries: dict[str, int] = {}
    years: dict[int, int] = {}
    for sample in samples:
        lineages[sample["l"]] = lineages.get(sample["l"], 0) + 1
        countries[sample["c"]] = countries.get(sample["c"], 0) + 1
        if sample["y"]:
            years[sample["y"]] = years.get(sample["y"], 0) + 1

    named = {k: v for k, v in lineages.items() if k != "—"}

    # ── tree ─────────────────────────────────────────────────────
    tree: dict | None = None
    if tree_path.is_file():
        tree = layout_tree(parse_newick(tree_path.read_text()))
        for node in tree["nodes"]:
            node["s"] = index.get(node["n"], -1) if node["n"] else -1

    ancestral: dict = {}
    if ancestral_path.is_file():
        with open(ancestral_path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                ancestral[row["node_id"]] = row.get("most_likely_state", "")

    # ── atlas ────────────────────────────────────────────────────
    loci: list[dict] = []
    atlas_path = Path(config.atlas_path) if config.atlas_path else None
    if atlas_path and atlas_path.is_file():
        with open(atlas_path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                loci.append(
                    {
                        "id": row["atlas_id"],
                        "start": int(row["genome_start"]),
                        "end": int(row["genome_end"]),
                        "strand": row["strand"],
                        "feature": row["gene_feature"],
                        "seq": row.get("sequence", ""),
                        "g4hunter": float(row["g4hunter_score"]) if row.get("g4hunter_score") else None,
                        "tools": int(row["concordant_tool_count"] or 0),
                        "topology": row.get("predicted_topology", ""),
                        "g4_type": row.get("g4_type", ""),
                        "gc_flank": float(row["gc_content_flanking"]) if row.get("gc_content_flanking") else None,
                        "conservation": float(row["conservation_pct_phylo"])
                        if row.get("conservation_pct_phylo")
                        else None,
                        "confidence": row.get("structural_confidence", ""),
                        "context": row.get("functional_context", ""),
                    }
                )

    # ── composition ──────────────────────────────────────────────
    # Real base counts, so the inspector shows the genome rather than a
    # decorative bar chart.
    composition = {}
    reference_fasta = Path(config.reference_fasta) if config.reference_fasta else None
    if reference_fasta and reference_fasta.is_file():
        seq = "".join(read_fasta(reference_fasta).values()).upper()
        total = len(seq) or 1
        counts = {b: seq.count(b) for b in "ACGT"}
        composition = {
            "total": len(seq),
            "counts": counts,
            "percent": {b: round(counts[b] / total * 100, 2) for b in "ACGT"},
            "gc": round((counts["G"] + counts["C"]) / total * 100, 2),
            "other": len(seq) - sum(counts.values()),
        }
    for locus in loci:
        sq = (locus.get("seq") or "").upper()
        if sq:
            t = len(sq)
            locus["composition"] = {b: round(sq.count(b) / t * 100, 1) for b in "ACGT"}
            locus["gc"] = round((sq.count("G") + sq.count("C")) / t * 100, 1)

    # ── gate ─────────────────────────────────────────────────────
    gate = evaluate_gate(config.ledger_path, config.pathogen, operational_mode=config.operational_mode)

    # ── floor ────────────────────────────────────────────────────
    n_complete = sum(1 for s in samples if s["y"] and s["c"] != "—" and s["h"] != "—")
    raw_total = 0
    corpus_fasta = Path(config.corpus_sequences_fasta) if config.corpus_sequences_fasta else None
    if corpus_fasta and corpus_fasta.is_file():
        raw_total = sum(1 for line in corpus_fasta.open() if line.startswith(">"))

    floor = {
        "n_sequences_in_window": {"value": len(samples), "floor": FLOORS["n_sequences_in_window"], "unit": "seq"},
        "min_sequences_per_lineage": {
            "value": min(named.values()) if named else 0,
            "floor": FLOORS["min_sequences_per_lineage"],
            "unit": "seq",
            "which": min(named, key=named.get) if named else None,
        },
        "n_timepoints": {"value": len(years), "floor": FLOORS["n_timepoints"], "unit": "yr"},
        "metadata_completeness": {
            "value": round(n_complete / len(samples), 4) if samples else 0,
            "floor": FLOORS["metadata_completeness"],
            "unit": "frac",
        },
        "alignment_qc_pass_fraction": {
            "value": round(len(samples) / raw_total, 4) if raw_total else None,
            "floor": FLOORS["alignment_qc_pass_fraction"],
            "unit": "frac",
        },
    }

    return {
        "identity": {
            "pathogen": config.pathogen,
            "display_name": config.display_name,
            "genome_type": config.genome_type,
            "reference": config.reference_accession,
            "genome_length": _genome_length(config),
            "cds": _cds_bounds(config),
            "atlas_version": _atlas_version(loci, atlas_path),
            "lineage_field": config.lineage_field,
            "n_samples": len(samples),
            "n_raw": raw_total,
            "n_lineages": len(named),
            "n_countries": len([c for c in countries if c != "—"]),
            "period": [min(years), max(years)] if years else None,
            "operational_mode": config.operational_mode,
        },
        "samples": samples,
        "lineages": dict(sorted(lineages.items(), key=lambda kv: -kv[1])),
        "countries": dict(sorted(countries.items(), key=lambda kv: -kv[1])),
        "years": dict(sorted(years.items())),
        "tree": tree,
        "ancestral": ancestral,
        "loci": loci,
        "composition": composition,
        "gate": {
            "permission": gate.permission.value,
            "permitted": gate.permitted,
            "operational_mode": gate.operational_mode,
            "failing_checks": list(gate.failing_checks),
            "supported_loci": list(gate.supported_loci),
            "latest_run": gate.latest_timestamp,
            "n_ledger_rows": gate.n_ledger_rows,
            "explanation": gate.explain(),
        },
        "floor": floor,
        # Stage 1.5 and Stage 3, surfaced beside the data they qualify.
        # The recombination verdict in particular belongs next to the tree:
        # it is the check that decides whether a single tree describes this
        # corpus at all, and a reader should not have to know the log exists.
        "recombination": _recombination_or_none(config.pathogen),
        "variants": _variants_or_none(config.pathogen),
        # A weak clock undermines anything reading calendar time, and the
        # number otherwise lives in a file nobody opens.
        "molecular_clock": _clock_or_none(config.pathogen),
        # Where this pathogen's artifacts actually live, repo-relative.
        # The Nextflow panel supplies them as --atlas/--alignment/
        # --rooted_tree so the orchestrated run SKIPS the stages that would
        # rebuild them; without the paths it would silently supply nothing
        # and quietly rebuild a published alignment instead.
        "paths": _artifact_paths(config),
        # The per-locus D.H1 evidence behind the gate's one-word verdict.
        # `gate` says SUPPORTED or BLOCKED; this says which loci, on what
        # p-values, against how many control clades. A verdict resting on
        # one locus with a thin control arm and a verdict resting on
        # twenty look identical until the rows are shown.
        "dh1": _dh1_or_none(config),
        "observations": observe(samples, named, years, countries, loci, floor, gate),
    }


def _artifact_paths(config) -> dict:
    """Repo-relative paths to the artifacts a run may supply rather than rebuild.

    Only paths that actually exist are returned. A key present but absent
    from disk would be passed to Nextflow as ``--alignment <missing>`` and
    fail the run's input check, which is a worse outcome than rebuilding.
    """
    from .tracks import _paths

    try:
        resolved = _paths(config.pathogen.lower())
    except Exception:  # noqa: BLE001 - a pathogen with no corpus has no artifacts
        return {}

    # _paths has no rooted-tree entry, so look beside the tree artifacts
    # it does know about. The divergence-rooted tree is preferred: it is
    # the one D.H1 reads (revision log R-16).
    phylo_dir = Path(resolved["mldist"]).parent if resolved.get("mldist") else None
    stem = config.pathogen.lower()
    rooted = None
    if phylo_dir:
        for name in (f"{stem}_rooted_div.nwk", f"{stem}_rooted.nwk", f"{stem}_iqtree_rooted.nwk"):
            if (phylo_dir / name).is_file():
                rooted = phylo_dir / name
                break

    out: dict[str, str] = {}
    candidates = {
        "alignment": resolved.get("alignment"),
        "rooted_tree": rooted,
        "atlas": getattr(config, "atlas_path", None),
    }
    for name, path in candidates.items():
        if path and Path(path).is_file():
            try:
                out[name] = str(Path(path).resolve().relative_to(REPO_ROOT))
            except ValueError:
                out[name] = str(path)
    return out


def _dh1_or_none(config) -> dict | None:
    """Per-locus D.H1 rows from the latest run recorded for this pathogen.

    Latest run only, matching :func:`g4watch.gating.evaluate_gate` — the
    ledger is append-only so earlier runs stay on the record, but showing
    a superseded run beside the current verdict would invite reading the
    two as one result.

    Nothing here is recomputed. These are the rows the gate itself read,
    so the panel cannot disagree with the verdict above it.
    """
    try:
        rows = read_ledger(config.ledger_path, config.pathogen)
    except Exception:  # noqa: BLE001 - a missing or malformed ledger must not break the payload
        return None
    if not rows:
        return None

    latest_timestamp = max(row.get("timestamp", "") for row in rows)
    latest = [row for row in rows if row.get("timestamp", "") == latest_timestamp]

    def _num(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _int_or_none(value):
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    loci = [
        {
            "atlas_id": row.get("atlas_id", ""),
            "verdict": row.get("verdict", ""),
            # Empty, not zero. A locus halted at the minimum-data floor
            # has no p-value, and 0.0 is a p-value.
            "raw_p": _num(row.get("raw_p_value")),
            "gc_adjusted_p_fdr": _num(row.get("gc_adjusted_p_value_fdr")),
            "locus_rate": _num(row.get("locus_disruption_rate")),
            "control_rate": _num(row.get("control_disruption_rate")),
            "tested": bool((row.get("raw_p_value") or "").strip()),
            "minimum_data_passed": (row.get("minimum_data_passed") or "").strip().lower() == "true",
            "failing_checks": [c for c in (row.get("minimum_data_failing_checks") or "").split(";") if c],
            "underpowered": (row.get("underpowered") or "").strip().lower() == "true",
            # Control provenance. Empty for runs predating R-20, which is
            # the honest value — it was not recorded then — and the UI says
            # "not recorded" rather than implying zero controls.
            "n_controls": _int_or_none(row.get("n_controls")),
            "control_regions": [r for r in (row.get("control_regions") or "").split(";") if r],
            "n_controls_same_compartment": _int_or_none(row.get("n_controls_same_compartment")),
        }
        for row in latest
    ]
    verdicts: dict[str, int] = {}
    for locus in loci:
        verdicts[locus["verdict"]] = verdicts.get(locus["verdict"], 0) + 1

    supported = [locus for locus in loci if locus["verdict"] == "SUPPORTED"]
    return {
        "timestamp": latest_timestamp or None,
        "decision_rule": (config.raw.get("dh1_gate") or {}).get("decision_rule", "conjunction"),
        "alpha": (config.raw.get("dh1_gate") or {}).get("alpha"),
        "min_carriers": ((config.raw.get("dh1_gate") or {}).get("locus_selection") or {}).get("min_carriers"),
        "n_loci": len(loci),
        "n_tested": sum(1 for locus in loci if locus["tested"]),
        "verdict_counts": verdicts,
        "loci": sorted(loci, key=lambda r: (r["gc_adjusted_p_fdr"] is None, r["gc_adjusted_p_fdr"] or 0)),
        # Named so the UI can say it rather than leaving a reader to count
        # the rows: a verdict carried by one locus is a different claim
        # from the same verdict carried by twenty.
        "rests_on": [locus["atlas_id"] for locus in supported],
        "ledger_path": str(config.ledger_path),
    }


def _recombination_or_none(pathogen: str):
    from .tracks import recombination

    try:
        return recombination(pathogen)
    except Exception:  # noqa: BLE001 - a missing or malformed log must not break the payload
        return None


def _clock_or_none(pathogen: str):
    from .tracks import molecular_clock

    try:
        return molecular_clock(pathogen)
    except Exception:  # noqa: BLE001
        return None


def _variants_or_none(pathogen: str):
    from .tracks import variant_summary

    try:
        return variant_summary(pathogen)
    except Exception:  # noqa: BLE001
        return None


def _genome_length(config: PathogenConfig) -> int | None:
    raw = getattr(config, "raw", None)
    if isinstance(raw, dict):
        return (raw.get("reference") or {}).get("genome_length")
    return getattr(config, "genome_length", None)


def _cds_bounds(config: PathogenConfig) -> list[int] | None:
    raw = getattr(config, "raw", None)
    if isinstance(raw, dict):
        reference = raw.get("reference") or {}
        if reference.get("cds_start") and reference.get("cds_end"):
            return [reference["cds_start"], reference["cds_end"]]
    return None


def _atlas_version(loci: list[dict], atlas_path: Path | None) -> str:
    """The Atlas version, preferring what the records themselves record.

    This used to read only the filename, so an Atlas stored as
    ``atlas.tsv`` reported no version at all even though every record in
    it carries an ``atlas_version`` field. The filename is kept as a
    fallback for older files written before that field existed.
    """
    for locus in loci:
        version = str(locus.get("atlas_version") or "").strip()
        if version and version != "unversioned":
            return version
    if atlas_path is None:
        return ""
    match = re.search(r"v(\d+(?:\.\d+)*)", atlas_path.name)
    return match.group(1) if match else ""


# ── observation engine ──────────────────────────────────────────────
def observe(samples, named, years, countries, loci, floor, gate) -> list[dict]:
    """Derive the interpretation layer from the data, by explicit rule.

    Each observation states the rule that fired and the numbers behind
    it, so a reader can check the claim rather than trust it. Severity is
    structural: ``block`` means a gate condition is unmet, ``caution``
    means a result is admissible but weakly supported, ``note`` is
    descriptive.
    """
    out: list[dict] = []

    smallest = floor["min_sequences_per_lineage"]
    if smallest["value"] < smallest["floor"]:
        out.append(
            {
                "severity": "block",
                "domain": "power",
                "title": f"Lineage {smallest['which']} below the per-lineage floor",
                "detail": (
                    f"{smallest['which']} is represented by {smallest['value']} genomes against an Appendix C floor "
                    f"of {smallest['floor']}. Corpus size cannot compensate: the floor is per lineage."
                ),
                "rule": f"min_sequences_per_lineage = {smallest['value']} < {smallest['floor']}",
                "focus": {"type": "lineage", "value": smallest["which"]},
            }
        )

    unresolved = sum(1 for s in samples if s["l"] == "—")
    if unresolved:
        out.append(
            {
                "severity": "caution",
                "domain": "metadata",
                "title": f"{unresolved} genome{'s' if unresolved != 1 else ''} with no resolvable lineage",
                "detail": (
                    "Held as an explicit unknown rather than folded into a named lineage, which would inflate the "
                    "smallest per-lineage count and could convert INSUFFICIENT_DATA into a false pass."
                ),
                "rule": "the pathogen's lineage vocabulary resolved nothing",
                "focus": {"type": "lineage", "value": "—"},
            }
        )

    single_tool = [locus for locus in loci if locus["tools"] <= 1]
    if single_tool:
        out.append(
            {
                "severity": "caution",
                "domain": "structure",
                "title": f"All {len(single_tool)} Atlas loci rest on a single algorithm",
                "detail": (
                    "Published inter-algorithm discordance for G4 prediction runs 30-60%. A locus supported by one "
                    "tool is a candidate motif, not a demonstrated structure, and no experimental G4 validation "
                    "exists for this pathogen."
                ),
                "rule": "concordant_tool_count <= 1 for every locus",
                "focus": {"type": "locus", "value": None},
            }
        )

    conserved = [locus for locus in loci if locus["conservation"] is not None]
    if conserved:
        low = min(conserved, key=lambda locus: locus["conservation"])
        high = max(conserved, key=lambda locus: locus["conservation"])
        out.append(
            {
                "severity": "note",
                "domain": "evolution",
                "title": f"Conservation spans {low['conservation']:.1f}-{high['conservation']:.1f}% across loci",
                "detail": (
                    f"{low['id']} is the least conserved and {high['id']} the most. None approaches the level at "
                    "which a structural constraint would be the obvious explanation."
                ),
                "rule": "range(conservation_pct_phylo)",
                "focus": {"type": "locus", "value": high["id"]},
            }
        )

    utr = [locus for locus in loci if "UTR" in (locus["feature"] or "")]
    if utr and len(utr) < len(loci):
        out.append(
            {
                "severity": "note",
                "domain": "architecture",
                "title": f"{len(utr)} of {len(loci)} candidates fall in untranslated regions",
                "detail": (
                    "UTRs are a small fraction of the genome, so this is a positional skew worth testing rather than "
                    "reporting. It is not evidence of function on its own."
                ),
                "rule": "gene_feature contains 'UTR'",
                "focus": {"type": "region", "value": "UTR"},
            }
        )

    if years:
        recent = sorted(years.items())[-5:]
        share = sum(n for _, n in recent) / max(sum(years.values()), 1)
        if share > 0.5:
            out.append(
                {
                    "severity": "caution",
                    "domain": "sampling",
                    "title": f"{share * 100:.0f}% of genomes come from the last five sampled years",
                    "detail": (
                        "Outbreak-driven sequencing concentrates recent time points. Temporal comparisons across the "
                        "full span are not sampling-balanced."
                    ),
                    "rule": f"share of {recent[0][0]}-{recent[-1][0]} = {share:.2f} > 0.50",
                    "focus": {"type": "period", "value": [recent[0][0], recent[-1][0]]},
                }
            )

    real_countries = {k: v for k, v in countries.items() if k != "—"}
    if real_countries:
        total = sum(real_countries.values())
        top, top_n = max(real_countries.items(), key=lambda kv: kv[1])
        entropy = -sum((n / total) * math.log(n / total) for n in real_countries.values())
        out.append(
            {
                "severity": "note",
                "domain": "geography",
                "title": f"{len(real_countries)} countries, Shannon H = {entropy:.2f}",
                "detail": f"{top} contributes the largest share at {top_n} genomes ({top_n / total * 100:.0f}%).",
                "rule": "Shannon entropy over country counts",
                "focus": {"type": "country", "value": top},
            }
        )

    if not gate.permitted:
        out.append(
            {
                "severity": "block",
                "domain": "gate",
                "title": "Stage 5 scoring is not authorised",
                "detail": gate.explain(),
                "rule": f"gate = {gate.permission.value}",
                "focus": {"type": "gate", "value": None},
            }
        )

    order = {"block": 0, "caution": 1, "note": 2}
    out.sort(key=lambda o: order[o["severity"]])
    return out
