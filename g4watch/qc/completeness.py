"""Is this sequence a whole genome, or a piece of one?

WHY THIS IS A SEPARATE QUESTION FROM QC. ``sequence_qc`` already computes
a completeness fraction, but only to pass or fail a sequence against a
threshold — the number is used as a gate and then discarded. Nothing
recorded whether a sequence that PASSED was a complete genome or a long
fragment, and nothing downstream could tell.

That distinction is not cosmetic here. Every Atlas coordinate is
reference-relative, and the surveillance metrics read tip states at those
coordinates. A partial sequence aligned to the reference is mostly gap,
and ``classify_tip_state`` correctly returns UNKNOWN across the regions
it does not cover — so a corpus of fragments does not produce false
disruption, it produces silent absence of evidence. The failure mode is
an analysis that looks fully powered and is not, because half its
genomes never spanned the loci being tested.

So completeness is classified, recorded, and shown. A user is told
whether they are analysing complete genomes or genomic regions, and the
answer travels with the analysis rather than living in a threshold.

WHAT THE CATEGORIES MEAN

    complete   >= COMPLETE_MIN of the reference length, contiguous
    partial    a real sequence covering materially less than the reference
    fragment   very short relative to the reference — a gene or a segment

A segmented virus is deliberately NOT handled by guessing. Influenza's
segments are each "partial" against a whole-genome reference and that is
the honest answer for a per-segment reference-relative pipeline; the
right fix is a per-segment reference, declared in config, not a
heuristic here that decides a 2.3 kb sequence "is probably segment 4".
"""

from __future__ import annotations

from dataclasses import dataclass

#: At or above this fraction of the reference, a sequence is treated as a
#: complete genome. 0.95 rather than 1.0 because real complete genomes are
#: routinely a little short: terminal repeats are hard to assemble and
#: submitters trim primer regions.
COMPLETE_MIN = 0.95

#: Below this, a sequence is a fragment rather than a partial genome — a
#: gene or a segment, not a genome missing some pieces.
FRAGMENT_MAX = 0.50

COMPLETE = "complete"
PARTIAL = "partial"
FRAGMENT = "fragment"

#: What a whole corpus is, when its sequences disagree.
MIXED = "mixed"


@dataclass(frozen=True)
class Completeness:
    """One sequence's completeness relative to the reference."""

    category: str
    fraction: float
    length: int
    reference_length: int
    #: Callable bases as a fraction of the sequence's own length. A
    #: sequence can be reference-length and still mostly N — the EBV
    #: corpus has whole submission batches like that — so length alone
    #: does not establish coverage.
    callable_fraction: float | None = None

    @property
    def is_complete(self) -> bool:
        return self.category == COMPLETE

    def describe(self) -> str:
        pct = f"{self.fraction:.0%}"
        if self.category == COMPLETE:
            return f"Complete genome ({pct} of reference, {self.length:,} nt)"
        if self.category == PARTIAL:
            return f"Partial genome ({pct} of reference, {self.length:,} nt)"
        return f"Genomic region / fragment ({pct} of reference, {self.length:,} nt)"


def classify_length(
    seq_length: int,
    reference_length: int,
    *,
    callable_fraction: float | None = None,
) -> Completeness:
    """Classify one sequence. Raises on a non-positive reference length."""
    if reference_length <= 0:
        raise ValueError("reference_length must be positive")
    fraction = seq_length / reference_length
    if fraction >= COMPLETE_MIN:
        category = COMPLETE
    elif fraction >= FRAGMENT_MAX:
        category = PARTIAL
    else:
        category = FRAGMENT
    return Completeness(
        category=category,
        fraction=round(fraction, 4),
        length=seq_length,
        reference_length=reference_length,
        callable_fraction=callable_fraction,
    )


def classify_sequence(sequence: str, reference_length: int) -> Completeness:
    """Classify from the sequence itself, measuring callable bases too.

    Length is taken WITHOUT gap characters: an aligned sequence padded to
    reference length is not thereby a complete genome, and counting the
    padding would make every fragment in an alignment look complete.
    """
    ungapped = sequence.replace("-", "").replace(".", "")
    callable_bases = sum(1 for b in ungapped.upper() if b in "ACGTU")
    return classify_length(
        len(ungapped),
        reference_length,
        callable_fraction=round(callable_bases / len(ungapped), 4) if ungapped else 0.0,
    )


@dataclass(frozen=True)
class CorpusCompleteness:
    """What a whole corpus is made of."""

    category: str
    counts: dict[str, int]
    n_sequences: int
    median_fraction: float

    def describe(self) -> str:
        if self.category == MIXED:
            parts = ", ".join(f"{n} {k}" for k, n in self.counts.items() if n)
            return f"Mixed corpus ({parts})"
        label = {COMPLETE: "Complete genomes", PARTIAL: "Partial genomes",
                 FRAGMENT: "Genomic regions / fragments"}[self.category]
        return f"{label} ({self.n_sequences:,} sequences)"

    @property
    def analysis_caveat(self) -> str | None:
        """What a reader must know before trusting a result from this corpus."""
        if self.category == COMPLETE:
            return None
        if self.category == MIXED:
            return (
                "This corpus mixes complete genomes with partial ones. Atlas coordinates are "
                "reference-relative, so a partial sequence contributes UNKNOWN — not evidence of "
                "loss — at every locus it does not span. Power therefore varies by locus in a way "
                "the sequence count does not show."
            )
        return (
            "These are not complete genomes. Every Atlas locus outside the region they cover is "
            "unassessable, so an analysis over them is narrower than the sequence count suggests."
        )


def classify_corpus(
    completenesses: list[Completeness],
    *,
    dominant_share: float = 0.9,
) -> CorpusCompleteness:
    """Summarise a corpus. ``mixed`` unless one category clearly dominates.

    A corpus is only called by a single category when at least
    ``dominant_share`` of its sequences fall in it. Calling a 60/40 split
    "complete" would be exactly the silent treatment of partial sequences
    as whole genomes this module exists to prevent.
    """
    counts = {COMPLETE: 0, PARTIAL: 0, FRAGMENT: 0}
    for item in completenesses:
        counts[item.category] += 1
    total = len(completenesses)
    if total == 0:
        return CorpusCompleteness(MIXED, counts, 0, 0.0)

    fractions = sorted(c.fraction for c in completenesses)
    median = fractions[total // 2]
    category = MIXED
    for name, n in counts.items():
        if n / total >= dominant_share:
            category = name
            break
    return CorpusCompleteness(category, counts, total, round(median, 4))
