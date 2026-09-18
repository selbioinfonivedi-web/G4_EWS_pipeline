"""FMDV polyprotein cleavage-product coordinates, mapped onto AY593823.1.

Built to test one specific external claim (Sindhu, ICAR-NIVEDI concept
deck, 2026): that G4 motifs in the 5' UTR and the 2B/2C non-structural
region ("sanctuary zones") stay conserved while motifs in the VP4-VP2-
VP3-VP1 capsid region are lost during outbreak years, as the surface
antigen escapes immunity. That is a genuine, testable directional claim
about WHERE in the genome disruption differs, and this module exists
only to answer "where in the genome is each Atlas locus" so D.H1 can be
run separately on each compartment -- with its own matched controls and
GC-adjustment, not the raw motif counts the original claim used.

THE PROBLEM THIS SOLVES. AY593823.1 (the FMDV2026 reference) carries a
single undivided ``CDS (polyprotein)`` feature -- no mat_peptide records,
so ``gene_feature`` in the Atlas cannot distinguish VP1 from 2C. Most
FMDV GenBank submissions have the same gap; only some carry the full
12-product cleavage map (L, VP4, VP2, VP3, VP1, 2A, 2B, 2C, 3A, 3B, 3C,
3D) as explicit mat_peptide features.

THE METHOD. 140 genomes in the FMDV corpus (``fmdv_complete_genomes.gb``)
carry all 12 mat_peptide features. Each was aligned onto AY593823.1 with
the same ``mafft --keeplength --addfragments`` invocation Stage 1 already
uses for reference-anchored alignment, which guarantees the output has
exactly len(AY593823.1) columns and column i IS AY593823.1 position i.
Walking each genome's own (ungapped) cleavage-site positions against that
alignment gives that genome's opinion of where each product falls in
AY593823.1's coordinates. 132 of 140 resolved cleanly; the coordinates
below are the MEDIAN across those 132, not any single genome's opinion.

Two genomes were checked by hand before trusting the batch: MF372126.1
and PX864607.1, aligned individually first. Their liftovers of the L
start disagreed by 72 nt -- not plausible as strain variation, since
every well-annotated FMDV genome has L begin exactly at the CDS start (a
biological fact confirmed across all 140: gap is 0 in every case), and
AY593823.1's own CDS start (1099, config/fmdv2026.yaml) sits 14 nt from
PX864607's liftover and 86 nt from MF372126's. MF372126's 5' UTR is
independently annotated as containing a poly-C tract -- FMDV's
notoriously variable-length, low-complexity repeat -- and a generic
aligner misplacing that one repeat would explain a roughly CONSTANT
offset propagating through every downstream boundary, which is what was
observed. Median-of-132 is the fix: an isolated misalignment like that is
exactly what a median should reject, and did -- the consensus L start
(1101) needed only a 2 nt reconciliation against AY593823.1's own known
CDS start of 1099, not 86.

WHAT IS AND ISN'T PRE-SPECIFIED HERE. The compartment split itself
(sanctuary = 5' UTR + 2B + 2C, capsid = VP4-VP2-VP3-VP1) is fixed by
Sindhu's independently-stated hypothesis, not chosen after looking at
G4-WATCH's own D.H1 p-values -- that is what makes running it a fair
test rather than a fishing expedition. But it IS a second, non-nested
test alongside the pathogen-wide D.H1 run already in the ledger, so a
result here should be reported as its own row (see cli.py's --lineage
convention: a stratified run is recorded as PATHOGEN:LABEL, and cannot by
itself open the pathogen's main gate) rather than silently replacing or
averaging into it.

WHERE THE NUMBERS LIVE. ``CONSENSUS_CLEAVAGE_SITES`` is the per-product
(start, end) table, 1-based inclusive, in AY593823.1 coordinates. Every
boundary is derived from the median START only -- product N's end is
fixed as (product N+1's median start) - 1, rather than also taking an
independent median of product N's own observed ends. The two do not
agree exactly (start-medians and end-medians are separate empirical
statistics, even for boundaries that are truly adjacent in every single
genome), and an early version of this table used both, which left the
3B/3C junction with a 1 nt gap (6050 -> 6052) that belonged to neither
product -- caught by
``test_the_twelve_products_are_contiguous_and_non_overlapping``.
Deriving every end from the next start makes exact adjacency a property
of the construction, not something to re-verify by eye.

The IQR alongside these medians in the derivation notebook (not shipped
here) showed the VP1/2A and 3D/3'UTR boundaries carrying the widest
spread (up to ~48 nt) -- both junctions with known biological length
variability (2A is ~16 residues and its exact cleavage point is
strain-sensitive; the 3' end abuts the 3'UTR, whose own length varies).
Loci near those two boundaries specifically should be read as
compartment-ambiguous even though the classifier below assigns them a
definite side.
"""

from __future__ import annotations

#: 1-based inclusive (start, end) in AY593823.1 coordinates. Median across
#: 132 independently mat_peptide-annotated FMDV genomes; see module
#: docstring for the alignment method and the cross-check that caught one
#: bad liftover before trusting the batch.
CONSENSUS_CLEAVAGE_SITES: dict[str, tuple[int, int]] = {
    "5UTR": (1, 1100),
    "L": (1101, 1703),
    "VP4": (1704, 1958),
    "VP2": (1959, 2615),
    "VP3": (2616, 3275),
    "VP1": (3276, 3915),
    "2A": (3916, 3964),
    "2B": (3965, 4426),
    "2C": (4427, 5380),
    "3A": (5381, 5838),
    "3B": (5839, 6051),
    "3C": (6052, 6690),
    "3D": (6691, 8100),
}

#: The two compartments Sindhu's hypothesis names. "Other" (L, 2A, the
#: 3-series, 3' UTR) is neither -- her claim makes no prediction there,
#: so a locus landing in it is excluded from both tests rather than
#: forced into the nearer one.
SANCTUARY_PRODUCTS = ("5UTR", "2B", "2C")
CAPSID_PRODUCTS = ("VP4", "VP2", "VP3", "VP1")


def _region_of(product: str) -> tuple[int, int]:
    return CONSENSUS_CLEAVAGE_SITES[product]


def classify_compartment(genome_start: int, genome_end: int) -> str:
    """'sanctuary', 'capsid', or 'other' for an AY593823.1-coordinate span.

    A locus overlapping both a sanctuary and a capsid region (possible
    only right at a junction, given the compartments are non-adjacent
    except through the intervening 'other' products) is reported as
    'ambiguous' rather than assigned to either -- silently picking a side
    at a boundary is exactly the kind of small decision that should not
    be invisible.
    """
    sanctuary_regions = [_region_of(p) for p in SANCTUARY_PRODUCTS]
    capsid_regions = [_region_of(p) for p in CAPSID_PRODUCTS]

    def overlaps(regions: list[tuple[int, int]]) -> bool:
        return any(genome_start <= r_end and genome_end >= r_start for r_start, r_end in regions)

    in_sanctuary = overlaps(sanctuary_regions)
    in_capsid = overlaps(capsid_regions)
    if in_sanctuary and in_capsid:
        return "ambiguous"
    if in_sanctuary:
        return "sanctuary"
    if in_capsid:
        return "capsid"
    return "other"
