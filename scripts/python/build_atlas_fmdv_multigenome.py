#!/usr/bin/env python3
"""Sprint 0/2 scale-up: rebuilds the FMDV G4 Reference Atlas from the
curated 5-10 genome reference set (data/reference_genomes/fmdv/
curated_reference_set.tsv), not just the single AY593823 bootstrap genome
`build_atlas_fmdv.py` used. Then re-populates conservation_pct_phylo on the
resulting candidate set exactly as `populate_conservation_fmdv.py` did,
producing the final `G4_Reference_Atlas_v1.0.fmdv.tsv`.

Design, and why it's structured this way:

1. Each curated genome is scanned independently, in its own native
   coordinates, via the existing, unmodified `atlas.stage0.scan_genome_stage0`
   -- no prediction logic is duplicated or reimplemented here.
2. Every curated genome is already a member of the real, QC-passed corpus
   alignment (`corpus/aligned/fmdv_qc_passed_aligned_to_ref.fasta`), so each
   non-reference genome's candidate positions are projected onto AY593823's
   coordinate system via that EXISTING alignment (`atlas.multi_genome`) --
   no new alignment run is needed.
3. Projected candidates are merged by reference-coordinate overlap. Where
   AY593823 itself already has a native candidate at that locus, the
   AY593823 record stays canonical (same position/sequence/score as before)
   and gains a cross-genome support tally in its evidence_note -- this is
   deliberately NOT folded into `concordant_tool_count`, which means
   *algorithm* concordance (G4Hunter + pattern-motif on ONE sequence), not
   *genome* concordance; conflating the two would repeat exactly the kind
   of construct-conflation this codebase's own confidence.py module exists
   to prevent.
4. Where a projected candidate has NO overlapping AY593823-native candidate
   (a PQS-forming sequence in another strain that AY593823 itself doesn't
   have), AY593823's own sequence at the projected span (+-100nt flank) is
   independently re-scanned with the SAME `scan_genome_stage0` pipeline. If
   AY593823's own sequence reproduces a qualifying hit there, it becomes a
   new, freshly-scored Atlas record (never inheriting the other genome's
   score). If it does not, the locus is reported in a side file, NOT
   silently added to the reference-coordinate Atlas -- claiming "AY593823
   has a PQS here" would be false if AY593823's own sequence doesn't
   actually form one.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from Bio import Phylo  # noqa: E402

from g4watch.atlas.confidence import structural_confidence  # noqa: E402
from g4watch.atlas.io import write_atlas_tsv  # noqa: E402
from g4watch.atlas.multi_genome import build_coordinate_map, project_reference_span  # noqa: E402
from g4watch.atlas.schema import AtlasCandidate, AtlasRecord  # noqa: E402
from g4watch.atlas.stage0 import GenomeAnnotation, scan_genome_stage0  # noqa: E402
from g4watch.io.fasta import read_fasta  # noqa: E402
from g4watch.metrics.conservation import g4c_phylo, g4d_phylo, naive_tip_proportion  # noqa: E402
from g4watch.metrics.tip_state_classifier import TipState, classify_tip_state  # noqa: E402
from g4watch.phylo.ancestral_states import reconstruct_ancestral_states  # noqa: E402
from g4watch.phylo.clade_collapse import collapse_to_maximal_clades  # noqa: E402

FMDV_DIR = REPO_ROOT / "data" / "reference_genomes" / "fmdv"
CURATED_SET_TSV = FMDV_DIR / "curated_reference_set.tsv"
ALIGNED_FASTA = FMDV_DIR / "corpus" / "aligned" / "fmdv_qc_passed_aligned_to_ref.fasta"
ROOTED_TREE = FMDV_DIR / "corpus" / "phylogenetics" / "fmdv_iqtree_rooted.nwk"
ATLAS_OUT = REPO_ROOT / "data" / "atlases" / "G4_Reference_Atlas_v1.0.fmdv.tsv"
UNCONFIRMED_OUT = REPO_ROOT / "data" / "atlases" / "fmdv_non_reference_transferable_loci.tsv"

REFERENCE_ACCESSION = "AY593823.1"
DISCOVERY_FLANK = 100
# G4Hunter reports the merged span of qualifying windows, so a real,
# homologous PQS's reported boundaries shift by up to roughly half a
# window's width between strains that differ by even a single G-run-length
# indel nearby -- verified against this real data (e.g. AY593823's own
# locus at nt 2692-2724 projects, across the 7 other curated genomes, to a
# smear of hits spanning roughly nt 2696-2768, not one fixed span). Two
# spans are treated as the "same" locus if they overlap at all, OR the gap
# between them is small relative to a typical motif length -- a documented
# judgment call, not a literature-derived constant.
MERGE_GAP_TOLERANCE = 15


def _load_curated_set() -> list[dict]:
    with open(CURATED_SET_TSV) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _load_native_sequence(accession_no_version: str) -> str:
    lines = (FMDV_DIR / f"{accession_no_version}.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def _overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """Same convention as g4prediction/concordance.py: overlap length over
    the shorter region's length, 1-based inclusive spans."""
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    shorter = min(a_end - a_start + 1, b_end - b_start + 1)
    return overlap / shorter if shorter > 0 else 0.0


def _same_locus(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """True if two reference-coordinate spans overlap at all, or are
    separated by no more than MERGE_GAP_TOLERANCE -- see that constant's
    docstring for why a fixed overlap-fraction threshold isn't used here."""
    if _overlap_fraction(a_start, a_end, b_start, b_end) > 0:
        return True
    gap = max(a_start, b_start) - min(a_end, b_end) - 1
    return gap <= MERGE_GAP_TOLERANCE


def scan_curated_genomes(curated_rows: list[dict]) -> dict[str, list[AtlasRecord]]:
    """Runs Stage 0 independently on every curated genome, in each genome's
    own native coordinates. Returns {accession: [AtlasRecord, ...]}."""
    per_genome: dict[str, list[AtlasRecord]] = {}
    for row in curated_rows:
        accession = row["accession"]
        accession_no_version = accession.split(".")[0]
        sequence = _load_native_sequence(accession_no_version)
        annotation = GenomeAnnotation(cds_start=int(row["cds_start"]), cds_end=int(row["cds_end"]))
        records = scan_genome_stage0(
            sequence,
            virus="FMDV",
            reference_accession=accession,
            atlas_version="v0.1-preconservation-multigenome-scan",
            annotation=annotation,
        )
        per_genome[accession] = records
        print(f"{accession} ({row['serotype']}, {row['country']}): {len(records)} native-coordinate candidates")
    return per_genome


class _Cluster:
    def __init__(self, ref_start: int, ref_end: int, is_reference_native: bool = False) -> None:
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.is_reference_native = is_reference_native
        self.members: list[tuple[str, AtlasRecord]] = []  # (accession, native record)

    def matches(self, ref_start: int, ref_end: int) -> bool:
        return _same_locus(self.ref_start, self.ref_end, ref_start, ref_end)

    def absorb(self, accession: str, record: AtlasRecord, ref_start: int, ref_end: int) -> None:
        self.members.append((accession, record))
        # Reference-native clusters keep AY593823's own native span fixed
        # (that IS the locus) -- only "novel" clusters (no reference member)
        # grow their envelope as more genomes' hits are added, since there
        # is no single authoritative span for a locus AY593823 itself
        # doesn't have a native candidate at.
        if not self.is_reference_native:
            self.ref_start = min(self.ref_start, ref_start)
            self.ref_end = max(self.ref_end, ref_end)


def project_and_cluster(
    per_genome: dict[str, list[AtlasRecord]],
    reference_native_seq: str,
    aligned: dict[str, str],
) -> tuple[list[_Cluster], list[dict]]:
    """Projects every curated genome's native-coordinate candidates onto
    AY593823's coordinate system, in two passes:

    1. Every projected candidate is first checked against AY593823's own
       FIXED native loci (Sprint 2's original 4) -- these anchor clusters
       are matched against, never grown, so a chain of slightly-shifted
       cross-strain hits can't drift a cluster away from its true locus
       (found the hard way: an earlier version grew reference clusters too,
       which let one genome's shifted hit widen a cluster enough for an
       unrelated nearby genome's hit to also match it, corrupting the
       merge).
    2. Anything left over (no AY593823-native locus nearby) is clustered
       among itself -- these are the "novel, discovered-via-another-genome"
       candidates handled by `build_evidence_and_records`.

    Returns (clusters, unmapped_report_rows)."""
    aligned_reference = aligned[REFERENCE_ACCESSION]
    reference_clusters = [
        _Cluster(record.genome_start, record.genome_end, is_reference_native=True)
        for record in per_genome[REFERENCE_ACCESSION]
    ]
    for cluster, record in zip(reference_clusters, per_genome[REFERENCE_ACCESSION]):
        cluster.absorb(REFERENCE_ACCESSION, record, record.genome_start, record.genome_end)

    novel_clusters: list[_Cluster] = []
    unmapped_rows: list[dict] = []

    for accession, records in per_genome.items():
        if accession == REFERENCE_ACCESSION:
            continue
        aligned_query = aligned[accession]
        coord_map = build_coordinate_map(aligned_query, aligned_reference)

        for record in records:
            projected = project_reference_span(record.genome_start, record.genome_end, coord_map)
            if projected is None:
                unmapped_rows.append(
                    {
                        "origin_accession": accession,
                        "native_start": record.genome_start,
                        "native_end": record.genome_end,
                        "reason": "span not mappable onto AY593823 (falls in a region unique to this genome)",
                    }
                )
                continue
            ref_start, ref_end = projected

            ref_match = next((c for c in reference_clusters if c.matches(ref_start, ref_end)), None)
            if ref_match is not None:
                ref_match.absorb(accession, record, ref_start, ref_end)
                continue

            novel_match = next((c for c in novel_clusters if c.matches(ref_start, ref_end)), None)
            if novel_match is None:
                novel_match = _Cluster(ref_start, ref_end, is_reference_native=False)
                novel_clusters.append(novel_match)
            novel_match.absorb(accession, record, ref_start, ref_end)

    return reference_clusters + novel_clusters, unmapped_rows


def build_evidence_and_records(
    clusters: list[_Cluster],
    reference_native_seq: str,
    n_curated_genomes: int,
) -> tuple[list[AtlasRecord], list[dict]]:
    """For each cluster: if AY593823 is a member, keep its native record as
    canonical and annotate cross-genome support. If not, re-scan AY593823's
    own sequence at the projected span (+-flank) fresh; only becomes an
    Atlas record if AY593823's own sequence independently reproduces a
    qualifying hit there AND that hit isn't just a re-detection of one of
    AY593823's own already-accepted native loci (the discovery flank can
    overlap a neighboring true locus -- guarded against explicitly)."""
    final_records: list[AtlasRecord] = []
    not_reproduced_rows: list[dict] = []
    reference_native_spans = [
        (c.ref_start, c.ref_end) for c in clusters if c.is_reference_native
    ]

    for cluster in clusters:
        ref_member = next((rec for acc, rec in cluster.members if acc == REFERENCE_ACCESSION), None)
        n_genomes_with_hit = len({acc for acc, _ in cluster.members})
        n_genomes_concordant = len({acc for acc, rec in cluster.members if rec.concordant_tool_count >= 2})
        support_note = (
            f"Cross-genome check across {n_curated_genomes} curated reference genomes "
            f"(O/A/Asia1; see data/reference_genomes/fmdv/curated_reference_set.tsv): "
            f"G4Hunter-predicted at the homologous position in {n_genomes_with_hit}/{n_curated_genomes} genomes; "
            f"pattern-motif-concordant (2-tool) in {n_genomes_concordant}/{n_curated_genomes} genomes."
        )

        if ref_member is not None:
            final_records.append(
                replace(ref_member, evidence_note=ref_member.evidence_note + " | " + support_note)
            )
            continue

        # Discovery case: no AY593823-native candidate at this locus.
        # Re-scan AY593823's own sequence at the projected span, fresh.
        n = len(reference_native_seq)
        flank_start = max(0, cluster.ref_start - 1 - DISCOVERY_FLANK)
        flank_end = min(n, cluster.ref_end + DISCOVERY_FLANK)
        local_seq = reference_native_seq[flank_start:flank_end]

        local_records = scan_genome_stage0(
            local_seq,
            virus="FMDV",
            reference_accession=REFERENCE_ACCESSION,
            atlas_version="v0.1-preconservation-multigenome-scan",
        )
        # Shift local (fragment-relative) coordinates back to full-genome
        # AY593823 coordinates and keep only a hit that actually overlaps
        # the span the other genome(s) implicated.
        reproduced = None
        best_overlap = 0.0
        for local_record in local_records:
            shifted_start = local_record.genome_start + flank_start
            shifted_end = local_record.genome_end + flank_start
            # Skip anything that's really just a re-detection of one of
            # AY593823's own already-accepted native loci, caught in the
            # flank of this different, nearby cluster.
            if any(
                _same_locus(shifted_start, shifted_end, native_start, native_end)
                for native_start, native_end in reference_native_spans
            ):
                continue
            overlap = _overlap_fraction(shifted_start, shifted_end, cluster.ref_start, cluster.ref_end)
            if overlap > 0 and overlap > best_overlap:
                best_overlap = overlap
                reproduced = replace(
                    local_record,
                    genome_start=shifted_start,
                    genome_end=shifted_end,
                    evidence_note=(
                        local_record.evidence_note
                        + f" | Locus discovered via curated-set genome(s) "
                        f"{sorted({acc for acc, _ in cluster.members})}, not natively present as a candidate "
                        f"in AY593823's own Stage-0 scan; independently re-scanned and reproduced on "
                        f"AY593823's own sequence at this projected position. | " + support_note
                    ),
                )

        if reproduced is not None:
            final_records.append(reproduced)
        else:
            not_reproduced_rows.append(
                {
                    "reference_projected_start": cluster.ref_start,
                    "reference_projected_end": cluster.ref_end,
                    "origin_genomes": sorted({acc for acc, _ in cluster.members}),
                    "n_genomes_with_hit": n_genomes_with_hit,
                    "reason": "AY593823's own sequence at this position does not independently reproduce a "
                    "qualifying G4Hunter hit -- not added to the reference-coordinate Atlas.",
                }
            )

    final_records.sort(key=lambda r: r.genome_start)
    renumbered = [replace(r, atlas_id=f"FMDV-G4-{i:03d}") for i, r in enumerate(final_records, start=1)]
    return renumbered, not_reproduced_rows


def populate_conservation(records: list[AtlasRecord]) -> list[AtlasRecord]:
    """Same real pipeline as populate_conservation_fmdv.py: tip-state
    classification -> ancestral reconstruction -> clade collapse ->
    g4c_phylo/g4d_phylo -> recomputed structural_confidence."""
    aligned = read_fasta(ALIGNED_FASTA)
    reference_seq = aligned[REFERENCE_ACCESSION]
    tree = Phylo.read(str(ROOTED_TREE), "newick")

    updated: list[AtlasRecord] = []
    for locus in records:
        print(f"--- {locus.atlas_id} (nt {locus.genome_start}-{locus.genome_end}) ---")
        tip_states: dict[str, str] = {}
        for accession, seq in aligned.items():
            state = classify_tip_state(reference_seq, seq, locus.genome_start, locus.genome_end)
            tip_states[accession] = state.value

        if len(set(tip_states.values())) < 2:
            print("  Only one distinct state across all tips -- skipping reconstruction.")
            updated.append(locus)
            continue

        naive = naive_tip_proportion(tip_states, TipState.CONSERVED.value)
        ancestral_run = reconstruct_ancestral_states(ROOTED_TREE, tip_states)
        clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)
        conservation_result = g4c_phylo(clades)
        disruption_result = g4d_phylo(clades)

        print(f"  naive={naive:.4f} g4c_phylo={conservation_result.status}"
              f"({conservation_result.value}) g4d_phylo={disruption_result.status}({disruption_result.value})")

        new_conservation_pct = round(conservation_result.value * 100, 2) if conservation_result.status == "OK" else None

        candidate = AtlasCandidate(
            concordant_tool_count=locus.concordant_tool_count,
            g4hunter_score=locus.g4hunter_score,
            g4rna_screener_score=locus.g4rna_screener_score,
            pqsfinder_score=locus.pqsfinder_score,
            conservation_pct_phylo=new_conservation_pct,
            overlaps_annotated_functional_region=(locus.functional_context.name == "KNOWN_FUNCTIONAL"),
        )
        recomputed_confidence = structural_confidence(candidate)
        if recomputed_confidence != locus.structural_confidence:
            print(f"  structural_confidence: {locus.structural_confidence.name} -> {recomputed_confidence.name}")

        updated.append(
            replace(
                locus,
                conservation_pct_phylo=new_conservation_pct,
                structural_confidence=recomputed_confidence,
                atlas_version="v1.0",
                evidence_note=locus.evidence_note
                + f" | conservation_pct_phylo is phylogenetically-corrected "
                f"(n_informative_clades={conservation_result.n_informative_clades}), naive tip-proportion "
                f"would have read {round(naive * 100, 2)}%.",
            )
        )
    return updated


def _write_report_tsv(rows: list[dict], path: Path) -> None:
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    curated_rows = _load_curated_set()
    print(f"Curated reference set: {len(curated_rows)} genomes\n")

    per_genome = scan_curated_genomes(curated_rows)
    print()

    aligned = read_fasta(ALIGNED_FASTA)
    reference_native_seq = _load_native_sequence(REFERENCE_ACCESSION.split(".")[0])

    clusters, unmapped_rows = project_and_cluster(per_genome, reference_native_seq, aligned)
    print(f"\n{len(clusters)} candidate clusters after reference-coordinate merge "
          f"({len(unmapped_rows)} candidates unmappable onto AY593823).\n")

    merged_records, not_reproduced_rows = build_evidence_and_records(
        clusters, reference_native_seq, n_curated_genomes=len(curated_rows)
    )
    print(f"{len(merged_records)} records kept in the reference-coordinate Atlas "
          f"({len(not_reproduced_rows)} projected loci NOT reproduced on AY593823's own sequence).\n")

    final_records = populate_conservation(merged_records)

    write_atlas_tsv(final_records, ATLAS_OUT)
    print(f"\nWrote {ATLAS_OUT} ({len(final_records)} loci)")

    _write_report_tsv(unmapped_rows + not_reproduced_rows, UNCONFIRMED_OUT)
    print(f"Wrote {UNCONFIRMED_OUT} ({len(unmapped_rows) + len(not_reproduced_rows)} excluded/unmapped loci, "
          f"documented not fabricated)")

    n_sc_plus = sum(1 for r in final_records if r.is_scoring_eligible())
    print(f"\nScoring-eligible (SC/BC/EC) loci: {n_sc_plus} / {len(final_records)}")


if __name__ == "__main__":
    main()
