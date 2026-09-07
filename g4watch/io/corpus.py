"""Load a pathogen's corpus as analysis-ready samples.

This lives in the library rather than in the web app because the library
must not depend on the app. The CLI needs the same corpus the workstation
shows, and reaching into ``web.workstation`` to get it made ``g4watch``
unimportable outside the repository root -- the installed console script
failed with ``ModuleNotFoundError: No module named 'web'``.

The dependency now runs one way: ``web`` imports ``g4watch``, never the
reverse.

Only the metadata needed for the downstream chain is read here. Tree
layout and the derived genome tracks stay where they are; this is the
corpus, not the whole payload.

Metadata alone is not a corpus the surveillance layer can measure.
:func:`load_samples` returns genomes with empty ``states`` and no mutation
counts, which leaves four of the seven G.2 terms permanently ``None``.
:func:`load_annotated_samples` is the entry point that also reads the
alignment and the Atlas and fills those in -- see
``variants/corpus_states.py`` for what it derives and why.
"""

from __future__ import annotations

import csv
from pathlib import Path

from ..atlas.io import read_atlas_tsv
from ..config import PathogenConfig
from ..metrics.surveillance_metrics import Sample
from ..pipeline.stage45_dh1 import _extract_year
from ..qc.metadata_normalization import LineageVocabulary
from .fasta import read_fasta


def _year(raw: str) -> int | None:
    text = _extract_year(raw or "")
    if not text:
        return None
    value = int(text)
    return value if 1900 <= value <= 2100 else None


def aligned_path(config: PathogenConfig) -> Path | None:
    """The conventional alignment location for a pathogen, if it exists."""
    name = config.pathogen.lower()
    candidate = Path(config.corpus_metadata_tsv).parent / "aligned" / f"{name}_qc_passed_aligned_to_ref.fasta"
    return candidate if candidate.is_file() else None


def load_samples(
    config: PathogenConfig,
    *,
    restrict_to_aligned: bool = True,
    exclude_lineages: tuple[str, ...] = (),
) -> list[Sample]:
    """Read the corpus metadata into :class:`Sample` objects.

    ``restrict_to_aligned`` keeps only genomes present in the alignment,
    which is what every downstream stage actually operates on. Counting
    genomes the analysis never sees would overstate the corpus.

    ``exclude_lineages`` is added to whatever ``corpus.exclude_lineages``
    the config declares -- it never replaces it, so an exploratory command
    line cannot quietly re-admit a lineage the config excluded. It exists
    for a real case: a lineage that no longer circulates cannot be the
    subject of early warning, and carrying it only to fail the per-lineage
    floor blocks the analysis of everything that does circulate. Any
    exclusion is recorded alongside the result.
    """
    metadata = config.corpus_metadata_tsv
    if not metadata or not Path(metadata).is_file():
        return []

    keep: set[str] | None = None
    if restrict_to_aligned:
        path = aligned_path(config)
        if path is not None:
            keep = set(read_fasta(path))

    excluded = {e.strip().upper() for e in (*config.exclude_lineages, *exclude_lineages) if e.strip()}
    vocabulary = LineageVocabulary.from_config(config)
    field = config.lineage_field
    out: list[Sample] = []

    with open(metadata, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = (row.get("accession") or "").strip()
            if not accession or (keep is not None and accession not in keep):
                continue

            fields = [row.get(field, "")]
            fields += [row.get(name, "") for name in config.lineage_fallback_fields]
            lineage = vocabulary.resolve(*fields) or (row.get(field, "") or "").strip()

            if lineage.upper() in excluded:
                continue

            out.append(
                Sample(
                    accession=accession,
                    lineage=lineage or "—",
                    country=(row.get("country", "") or "").split(":")[0].strip() or "—",
                    year=_year(row.get("collection_date", "")),
                    states={},
                    g4_mutations=None,
                    total_mutations=None,
                )
            )
    return out


def lineage_counts(samples: list[Sample]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for sample in samples:
        counts[sample.lineage] = counts.get(sample.lineage, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))


def load_annotated_samples(
    config: PathogenConfig,
    *,
    exclude_lineages: tuple[str, ...] = (),
    include_ineligible_loci: bool = False,
    detect_gains: bool = True,
):
    """The corpus with G4 tip states and mutation counts filled in.

    This is what the surveillance layer needs and what :func:`load_samples`
    alone does not provide. Returns ``(samples, report)``; the report is
    ``None`` only when annotation could not be attempted at all, and says
    why in that case via :func:`annotation_blocked_reason`.

    Nothing is silently degraded: a missing alignment, a missing Atlas or
    an Atlas with no scoring-eligible locus each produce unannotated
    samples and a stated reason, never plausible-looking zeros.
    """
    from ..variants.corpus_states import annotate_samples

    samples = load_samples(config, exclude_lineages=exclude_lineages)
    reason = annotation_blocked_reason(config)
    if reason is not None:
        return samples, None

    g4 = config.raw.get("g4_prediction", {}).get("g4hunter", {}) or {}
    alignment = read_fasta(aligned_path(config))
    return annotate_samples(
        samples,
        alignment=alignment,
        reference_id=str(config.reference_accession),
        atlas=read_atlas_tsv(config.atlas_path),
        include_ineligible=include_ineligible_loci,
        detect_gains=detect_gains,
        g4hunter_window=int(g4.get("window", 25)),
        g4hunter_threshold=float(g4.get("threshold", 1.2)),
    )


def annotation_blocked_reason(config: PathogenConfig) -> str | None:
    """Why tip states cannot be assigned, or ``None`` if they can be.

    Checked before the work starts so the caller gets a sentence naming
    the missing artifact rather than an empty result to interpret.
    """
    path = aligned_path(config)
    if path is None:
        return (
            "no reference-pinned alignment was found for this pathogen; tip states are defined "
            "against alignment columns, so Stage 1 must run first"
        )
    if not config.atlas_path or not Path(config.atlas_path).is_file():
        return f"no G4 Reference Atlas at {config.atlas_path}; Stage 0 must run first"
    if not config.reference_accession:
        return "config sets no reference.accession, so there is no coordinate system to pin states to"
    if config.reference_accession not in read_fasta(path):
        return (
            f"reference {config.reference_accession} is not present in the alignment at {path}; "
            "the alignment was not built against the Atlas's reference"
        )
    return None
