"""Builds the G4-WATCH overview deck: how it works, and what it found.

Deliberately separate from ``ppt_assets/build_deck.py``, which builds the
long technical briefing. This one is short and answers two questions for
a reader who has not seen the code.

EVERY number on these slides is read from a file in this repository --
the Atlas TSVs, the append-only testing ledger, the pathogen configs --
by ``make_figures.py`` or by the constants block below, which cites its
source for each value. The technical briefing carries some deliberately
illustrative charts, labelled as such on their slides; none of them are
reused here, because a deck about what the pipeline found may only show
what it found.

The headline result is a negative one. That is the result.
"""
import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
DECK = os.path.join(OUT, "G4-WATCH_overview.pptx")

NAVY = RGBColor(0x1B, 0x3A, 0x5C)
TEAL = RGBColor(0x1F, 0x8A, 0x70)
AMBER = RGBColor(0xD9, 0x8C, 0x1A)
RED = RGBColor(0xC0, 0x39, 0x2B)
GREY = RGBColor(0x5A, 0x66, 0x70)
LGREY = RGBColor(0xEE, 0xF1, 0xF3)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
DARK = RGBColor(0x22, 0x2A, 0x30)
PALE = RGBColor(0x9F, 0xC7, 0xD8)

# ── facts, with their sources ───────────────────────────────────────
# data/atlases/G4_Reference_Atlas_v*.tsv and config/*.yaml, counted by
# make_figures.py; re-run it and these must still agree.
N_PATHOGENS = 11
N_FETCHED = 7289
N_ALIGNED = 6485
# data/atlases/testing_ledger.tsv, latest run per pathogen
FMDV_LOCI_TESTED = 37
FMDV_PAST_FLOOR = 18
FMDV_NOT_SUPPORTED = 15
FMDV_GC = 2
FMDV_OPPOSITE = 1
# the standout locus, from the Atlas TSV and the ledger row that tested it
G4_004 = {
    "id": "FMDV2026-G4-004",
    "span": "4313–4337",
    "feature": "CDS (polyprotein)",
    "type": "RNA_G4",
    "g4hunter": "-1.28",
    "conservation": "94.86%",
    "gc_flank": "55.56%",
    "locus_rate": "0.775",
    "control_rate": "0.288",
    "q": "1.08e-09",
}
N_TESTS = 1199

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]
SW, SH = prs.slide_width, prs.slide_height


def add_slide():
    return prs.slides.add_slide(BLANK)


def set_bg(slide, color=WHITE):
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = color


def add_rect(slide, x, y, w, h, color, line=None):
    shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, x, y, w, h)
    shp.fill.solid()
    shp.fill.fore_color.rgb = color
    if line is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line
        shp.line.width = Pt(1)
    shp.shadow.inherit = False
    return shp


def add_text(slide, x, y, w, h, text, size=18, color=DARK, bold=False, italic=False,
             align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, font="Calibri", line_spacing=1.0):
    box = slide.shapes.add_textbox(x, y, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = line_spacing
        run = p.add_run()
        run.text = line
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.italic = italic
        run.font.name = font
        run.font.color.rgb = color
    return box


def bullets(slide, x, y, w, h, items, size=16, color=DARK, space_after=9,
            font="Calibri", line_spacing=1.08):
    """items: (level, text) or (level, text, {overrides})."""
    box = slide.shapes.add_textbox(x, y, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    first = True
    for item in items:
        level, text = item[0], item[1]
        over = item[2] if len(item) > 2 else {}
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.level = level
        p.space_after = Pt(over.get("space_after", space_after))
        p.line_spacing = over.get("line_spacing", line_spacing)
        marker = "" if over.get("no_marker") else {0: "▪  ", 1: "–  ", 2: "·  "}.get(level, "")
        run = p.add_run()
        run.text = marker + text
        run.font.size = Pt(over.get("size", size - level * 1.5))
        run.font.bold = over.get("bold", False)
        run.font.italic = over.get("italic", False)
        run.font.name = over.get("font", font)
        run.font.color.rgb = over.get("color", color)
    return box


def header(slide, kicker, title, rule=TEAL):
    add_rect(slide, 0, 0, SW, Inches(1.02), NAVY)
    add_text(slide, Inches(0.55), Inches(0.11), Inches(11), Inches(0.3),
             kicker.upper(), size=12.5, color=PALE, bold=True)
    add_text(slide, Inches(0.55), Inches(0.37), Inches(12.2), Inches(0.6),
             title, size=26, color=WHITE, bold=True)
    add_rect(slide, 0, Inches(1.02), SW, Pt(3), rule)


def footer(slide, note, page):
    add_text(slide, Inches(0.55), Inches(7.12), Inches(11.2), Inches(0.3),
             note, size=9.5, color=GREY, italic=True)
    add_text(slide, Inches(12.2), Inches(7.12), Inches(0.6), Inches(0.3),
             str(page), size=9.5, color=GREY, align=PP_ALIGN.RIGHT)


def stat(slide, x, y, w, value, label, color=NAVY):
    add_text(slide, x, y, w, Inches(0.62), value, size=38, color=color, bold=True)
    add_text(slide, x, y + Inches(0.62), w, Inches(0.5), label, size=11.5, color=GREY)


def chip(slide, x, y, w, h, title, sub="", fill=NAVY, fg=WHITE, ts=12.5, ss=10):
    add_rect(slide, x, y, w, h, fill)
    if sub:
        add_text(slide, x + Inches(0.1), y + Inches(0.13), w - Inches(0.2), Inches(0.3),
                 title, size=ts, color=fg, bold=True, align=PP_ALIGN.CENTER)
        add_text(slide, x + Inches(0.08), y + Inches(0.44), w - Inches(0.16), h - Inches(0.5),
                 sub, size=ss, color=fg, align=PP_ALIGN.CENTER, line_spacing=0.95)
    else:
        add_text(slide, x + Inches(0.1), y, w - Inches(0.2), h, title, size=ts, color=fg,
                 bold=True, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)


ARROW_W = Inches(0.26)


def arrow(slide, x, y, w=ARROW_W):
    """A thin rule between stages.

    This was a ▶ glyph in a textbox, which LibreOffice renders as a filled
    amber block -- it applies the theme accent to a shape whose fill
    python-pptx left unset. A rectangle with an explicit fill has no such
    ambiguity, and the sequence is already numbered.
    """
    add_rect(slide, x + Inches(0.05), y + Inches(0.14), w - Inches(0.10), Pt(2.5),
             RGBColor(0xC3, 0xCC, 0xD2))


def pic(slide, name, x, y, w):
    return slide.shapes.add_picture(os.path.join(HERE, name), x, y, width=w)


# ══════════════════════════════════════════════════════ 1 · title ═══
s = add_slide(); set_bg(s, NAVY)
add_rect(s, 0, Inches(3.30), Inches(1.5), Pt(4), TEAL)
add_text(s, Inches(0.9), Inches(1.55), Inches(11.5), Inches(0.4),
         "ICAR–NIVEDI  ·  GENOMIC SURVEILLANCE", size=13, color=PALE, bold=True)
add_text(s, Inches(0.9), Inches(2.05), Inches(11.5), Inches(1.1),
         "G4-WATCH", size=54, color=WHITE, bold=True)
add_text(s, Inches(0.9), Inches(3.65), Inches(10.6), Inches(1.5),
         "Do G-quadruplex sites in livestock virus genomes behave differently\n"
         "from the rest of the genome — and can that be used as an early warning?",
         size=21, color=WHITE, line_spacing=1.28)
add_text(s, Inches(0.9), Inches(5.35), Inches(10.6), Inches(0.9),
         "How the pipeline works, and what it found across "
         f"{N_PATHOGENS} pathogens and {N_ALIGNED:,} genomes.",
         size=15, color=PALE, line_spacing=1.2)
add_text(s, Inches(0.9), Inches(6.62), Inches(11), Inches(0.4),
         "Research framework — not a validated diagnostic. No output here should drive a "
         "control decision on its own.", size=11, color=RGBColor(0x7E, 0xA3, 0xB8), italic=True)

# ══════════════════════════════════════════════ 2 · the idea ═══════
s = add_slide(); set_bg(s)
header(s, "the idea", "A G-quadruplex is a structure, not just a sequence")
bullets(s, Inches(0.55), Inches(1.42), Inches(6.5), Inches(4.7), [
    (0, "Four runs of guanine, separated by short loops, can fold the strand back on itself "
        "into a stack of G-tetrads held together by Hoogsteen bonds — a G-quadruplex, or G4."),
    (0, "Where one forms it obstructs polymerase and ribosome. A virus carrying one is under "
        "selection either to keep it or to lose it, so the site is a place worth watching."),
    (0, "Detection here is G4Hunter: a sliding window (w = 25) scoring G-skew from +4 to −4, "
        "with |score| ≥ 1.2 called, plus an independent regex motif. Two methods must overlap "
        "by ≥ 80% to count as concordant.", {"font": "Calibri"}),
    (0, "Negative scores are C-rich windows — a G4 on the opposite strand. They are kept, and "
        "strand is reported as a covariate rather than settled by exclusion."),
    (0, "The framework rests on one testable claim about these sites, stated on the next slide. "
        "It is stated as a claim precisely so that it can fail.",
        {"color": NAVY, "bold": True}),
], size=14.5, space_after=13)

add_rect(s, Inches(7.4), Inches(1.55), Inches(5.4), Inches(4.35), LGREY)
add_text(s, Inches(7.7), Inches(1.80), Inches(4.8), Inches(0.35),
         "WHAT THE PIPELINE STORES PER SITE", size=11.5, color=NAVY, bold=True)
bullets(s, Inches(7.7), Inches(2.25), Inches(4.8), Inches(3.4), [
    (0, "Position on the reference genome"),
    (0, "G4Hunter score and predicted topology"),
    (0, "Gene feature and strand"),
    (0, "GC content of the flanking sequence"),
    (0, "Conservation across the phylogeny"),
    (0, "A structural-confidence tier"),
], size=13, color=GREY, space_after=11)
add_text(s, Inches(7.7), Inches(4.95), Inches(4.8), Inches(0.8),
         "Together these are the G4 Reference Atlas. It is fixed before any\n"
         "statistical test runs, so the loci cannot be chosen after seeing the answer.",
         size=11.5, color=NAVY, italic=True, line_spacing=1.2)
footer(s, "Atlas schema: data/atlases/G4_Reference_Atlas_v*.tsv", 2)

# ══════════════════════════════════════════════ 3 · hypothesis ═════
s = add_slide(); set_bg(s)
header(s, "the claim under test", "D.H1 — a directional hypothesis", rule=AMBER)
add_rect(s, Inches(0.55), Inches(1.45), Inches(12.25), Inches(1.5), LGREY,
         line=RGBColor(0xCE, 0xD8, 0xDE))
add_text(s, Inches(0.95), Inches(1.78), Inches(11.5), Inches(1.0),
         "G4 loci are disrupted LESS often than matched control regions of the same genome —\n"
         "and that difference survives adjustment for GC content.",
         size=19.5, color=NAVY, bold=True, line_spacing=1.3)

add_text(s, Inches(0.55), Inches(3.25), Inches(12.2), Inches(0.35),
         "Three things have to be true at once, and each is a separate way to fail:",
         size=14, color=GREY)
for i, (title, body, colour) in enumerate([
    ("1 · A real difference",
     "Fisher exact test on the locus's clade-level disruption counts against those of its "
     "5 matched controls, then Benjamini–Hochberg across the whole pre-specified set.", NAVY),
    ("2 · In the right direction",
     "Lower, not merely different. Both tests are two-sided, so the sign is checked explicitly: "
     "a locus significantly MORE disrupted returns SIGNAL_OPPOSITE_DIRECTION.", AMBER),
    ("3 · Not just GC",
     "G4 sites are G-rich by definition and GC-rich sequence mutates differently. A pooled "
     "logistic model with flanking GC as a covariate must keep the locus term significant.", TEAL),
]):
    x = Inches(0.55) + i * Inches(4.16)
    add_rect(s, x, Inches(3.75), Inches(3.9), Pt(3.5), colour)
    add_text(s, x, Inches(3.95), Inches(3.9), Inches(0.35), title, size=15, color=colour, bold=True)
    add_text(s, x, Inches(4.38), Inches(3.9), Inches(1.7), body, size=12.5, color=DARK,
             line_spacing=1.22)

add_text(s, Inches(0.55), Inches(6.25), Inches(12.2), Inches(0.75),
         "Each locus therefore returns one of six verdicts, not a yes/no: SUPPORTED, "
         "NOT_SUPPORTED, SIGNAL_OPPOSITE_DIRECTION,\nSIGNAL_EXPLAINED_BY_GC, "
         "INSUFFICIENT_DATA, or UNDERPOWERED. The distinction between the last two and "
         "\"no effect\" is the point.", size=13, color=GREY, line_spacing=1.3)
footer(s, "Implemented in g4watch/validation/dh1_gate.py; every verdict is written to an append-only ledger", 3)

# ══════════════════════════════════════════════ 4 · pipeline ═══════
s = add_slide(); set_bg(s)
header(s, "how it works", "Seven stages, in a fixed order")
STAGES = [
    ("0", "Atlas", "predict G4 loci\non the reference"),
    ("1", "QC + align", "drop incomplete,\nundated genomes"),
    ("1.5", "Recombination", "PHI test —\nmandatory"),
    ("2", "Phylogeny", "ML tree,\ntime-scaled root"),
    ("3", "Variants", "call against\nthe reference"),
    ("4", "Metrics", "per-locus disruption,\nGC control gate"),
]
x = Inches(0.55)
for num, title, sub in STAGES:
    chip(s, x, Inches(1.55), Inches(1.72), Inches(1.15), f"{num} · {title}", sub)
    x += Inches(1.72)
    if num != "4":
        arrow(s, x, Inches(2.0)); x += Inches(0.3)

add_rect(s, Inches(0.55), Inches(3.05), Inches(12.25), Inches(0.95), AMBER)
add_text(s, Inches(0.9), Inches(3.19), Inches(11.5), Inches(0.7),
         "GATE  ·  D.H1 must return SUPPORTED, and the config must set operational_mode: true",
         size=16, color=WHITE, bold=True)
add_text(s, Inches(0.9), Inches(3.56), Inches(11.5), Inches(0.4),
         "Both conditions, checked in code against a persisted ledger — editing the workflow "
         "cannot open the gate.", size=11.5, color=WHITE)

chip(s, Inches(0.55), Inches(4.35), Inches(6.0), Inches(1.15),
     "5 · Scoring", "G4-EWS score, CUSUM / EWMA control limits, detection", fill=TEAL, ts=14, ss=11)
chip(s, Inches(6.82), Inches(4.35), Inches(5.98), Inches(1.15),
     "6 · Reporting", "report card, gate status, provenance — ALWAYS runs", fill=NAVY, ts=14, ss=11)

add_text(s, Inches(0.55), Inches(5.8), Inches(12.2), Inches(1.1),
         "Stage 1.5 has no skip flag. Reconstructing ancestral states across a recombinant "
         "alignment reconstructs a history that never happened, so the screen is\nmandatory for "
         "every pathogen. Stage 6 runs whatever the gate says: a blocked run still produces a "
         "report explaining what blocked it.",
         size=13, color=GREY, line_spacing=1.35)
footer(s, "Orchestrated with Nextflow DSL2. Tools: MAFFT · PhiPack · IQ-TREE 2 · TreeTime", 4)

# ══════════════════════════════════════════════ 5 · the cone ═══════
# The pipeline's defining property is attrition: 936 genomes in, 67
# candidate loci, 18 still testable, none supported. A box-and-arrow
# flowchart hides that; a cone makes it the subject.
s = add_slide(); set_bg(s)
header(s, "the workflow, in one figure", "Where everything goes")
# Full width: the per-stage methods down the right of the figure are the
# point of it, and they are unreadable at anything narrower.
pic(s, "fig_cone.png", Inches(0.42), Inches(1.20), Inches(12.5))
add_text(s, Inches(0.55), Inches(6.90), Inches(12.2), Inches(0.3),
         "Green bands are the statistical gates; red is the outcome. Two bands go UP — the "
         "alignment adds the reference, which is why 936 follows 935.",
         size=11, color=GREY, italic=True)
footer(s, "FMDV2026 throughout — one pathogen end to end, not totals that mix corpora stopped at different stages", 5)

# ══════════════════════════════════════════════ 5 · the floor ══════
s = add_slide(); set_bg(s)
header(s, "before any test runs", "A minimum-data floor, fixed in advance")
add_text(s, Inches(0.55), Inches(1.45), Inches(12.2), Inches(0.9),
         "A statistic computed on four sequences is not a weak result — it is not a result. "
         "Every locus is checked against these\nthresholds first, and a locus that fails them "
         "returns INSUFFICIENT_DATA rather than a number nobody should read.",
         size=14.5, color=DARK, line_spacing=1.3)
FLOORS = [
    ("30", "sequences in the\nanalysis window"),
    ("20", "sequences in every\nlineage"),
    ("3", "distinct sampling\ntimepoints"),
    ("90%", "metadata\ncompleteness"),
    ("3", "informative clades\nper group"),
]
x = Inches(0.55)
for value, label in FLOORS:
    add_rect(s, x, Inches(2.75), Inches(2.3), Inches(1.55), LGREY)
    add_text(s, x, Inches(2.92), Inches(2.3), Inches(0.6), value, size=33, color=NAVY,
             bold=True, align=PP_ALIGN.CENTER)
    add_text(s, x, Inches(3.55), Inches(2.3), Inches(0.7), label, size=11.5, color=GREY,
             align=PP_ALIGN.CENTER, line_spacing=1.1)
    x += Inches(2.47)

add_rect(s, Inches(0.55), Inches(4.75), Pt(4), Inches(1.7), TEAL)
add_text(s, Inches(0.95), Inches(4.78), Inches(11.8), Inches(1.7),
         "Why this is a threshold and not a list of exclusions\n"
         "Writing down which lineages to drop AFTER seeing the corpus is selection on the outcome. "
         "A threshold is a property of sampling,\nknowable before any test runs, and it lives in "
         "the config — so a run stays fully described by (commit, config, accession list).",
         size=13.5, color=DARK, line_spacing=1.4)
footer(s, "g4watch/validation/minimum_data_gate.py — thresholds from Build Architecture Appendix C", 6)

# ══════════════════════════════════════════════ 6 · corpora ════════
s = add_slide(); set_bg(s)
header(s, "what it was run on", "Eleven pathogens, real GenBank corpora")
pic(s, "fig_corpus.png", Inches(1.67), Inches(1.40), Inches(10.0))
add_text(s, Inches(0.55), Inches(5.95), Inches(12.2), Inches(1.0),
         "Ten livestock viruses across four genome types — positive- and negative-sense RNA, "
         "segmented dsRNA, and large dsDNA — plus\nEpstein–Barr virus as a cross-species check "
         "on a GC-rich human herpesvirus. No synthetic sequence appears in any of these counts; "
         "the\nsynthetic datasets used to test the software live separately under data/synthetic/ "
         "and are prefixed SYNTH-.",
         size=12.5, color=GREY, line_spacing=1.35)
footer(s, "Counted from config/*.yaml and the aligned FASTA on disk by overview_assets/make_figures.py", 7)

# ══════════════════════════════════════════ 7 · result · atlas ═════
s = add_slide(); set_bg(s)
header(s, "result 1", "The Atlas: plenty of predictions, little evidence")
pic(s, "fig_atlas.png", Inches(1.67), Inches(1.38), Inches(10.0))
add_text(s, Inches(0.55), Inches(5.95), Inches(12.2), Inches(1.0),
         "Prediction is cheap; evidence is not. A locus reaches SC only when its conservation "
         "across phylogenetically independent lineages\nsupports it, and that conservation is "
         "computed from the alignment — deliberately NOT from the disruption rate, which would "
         "make\neligibility restate the very outcome D.H1 is meant to test.",
         size=12.5, color=GREY, line_spacing=1.35)
footer(s, "Tiers: EC / BC / SC / MC / WC / AA. Scoring eligibility is SC and above", 8)

# ══════════════════════════════════════════ 8 · result · D.H1 ══════
s = add_slide(); set_bg(s)
header(s, "result 2", "The D.H1 test on its best-powered corpus", rule=RED)
pic(s, "fig_dh1.png", Inches(0.75), Inches(1.30), Inches(5.35))
bullets(s, Inches(6.75), Inches(1.60), Inches(6.05), Inches(5.0), [
    (0, f"{FMDV_LOCI_TESTED} loci were tested on FMDV2026 (936 genomes, the deepest corpus here). "
        f"{FMDV_LOCI_TESTED - FMDV_PAST_FLOOR} did not clear the data floor."),
    (0, f"Of the {FMDV_PAST_FLOOR} that did: {FMDV_NOT_SUPPORTED} NOT_SUPPORTED, "
        f"{FMDV_GC} explained by GC, {FMDV_OPPOSITE} significant in the opposite direction."),
    (0, "Zero SUPPORTED.", {"color": RED, "bold": True, "size": 19}),
    (0, "Most loci do sit below the diagonal — the direction D.H1 predicts — but none of those "
        "differences survive GC adjustment and FDR correction. A trend that does not reach "
        "significance is not a finding."),
    (0, "The one locus that IS significant sits above the diagonal: more disrupted than its "
        "controls, not less.", {"color": RED}),
], size=14.5, space_after=16)
footer(s, "data/atlases/testing_ledger.tsv, latest run per pathogen. q = Benjamini-Hochberg adjusted p", 9)

# ══════════════════════════════════════════ 9 · G4-004 ═════════════
s = add_slide(); set_bg(s)
header(s, "result 3", "The one locus that moved — and it moved the wrong way", rule=RED)
add_rect(s, Inches(0.55), Inches(1.40), Inches(12.25), Inches(1.28), LGREY,
         line=RGBColor(0xCE, 0xD8, 0xDE))
add_text(s, Inches(0.9), Inches(1.60), Inches(5.6), Inches(0.5), G4_004["id"],
         size=24, color=RED, bold=True, font="Consolas")
add_text(s, Inches(0.9), Inches(2.10), Inches(6.0), Inches(0.45),
         f"nt {G4_004['span']}  ·  {G4_004['feature']}  ·  {G4_004['type']}",
         size=13.5, color=GREY, font="Consolas")
# The q value is twice as wide as a bare rate, so it gets its own width
# and point size rather than wrapping out of the band.
for i, (v, lab, col) in enumerate([
    (G4_004["locus_rate"], "disruption, locus", RED),
    (G4_004["control_rate"], "disruption, controls", NAVY),
]):
    stat(s, Inches(6.85) + i * Inches(1.95), Inches(1.60), Inches(1.85), v, lab, color=col)
add_text(s, Inches(10.75), Inches(1.62), Inches(1.95), Inches(0.55),
         G4_004["q"], size=27, color=RED, bold=True)
add_text(s, Inches(10.75), Inches(2.18), Inches(2.0), Inches(0.5),
         "q, GC-adjusted + FDR", size=11.5, color=GREY)

bullets(s, Inches(0.55), Inches(3.05), Inches(6.0), Inches(3.6), [
    (0, "Highly conserved: " + G4_004["conservation"] + " across independent lineages — "
        "this is not a poorly aligned region."),
    (0, "Not a GC artefact: flanking GC is " + G4_004["gc_flank"] + ", and the signal "
        "STRENGTHENED when GC-matched controls were used."),
    (0, "Not a control artefact either. When control selection was rebuilt to stop it "
        "collapsing onto the genome's 5' end, p moved from 0.14 to 1.08e-09.",
        {"color": NAVY, "bold": True}),
], size=14, space_after=15)

add_rect(s, Inches(6.9), Inches(3.05), Inches(5.92), Inches(2.35), NAVY)
add_text(s, Inches(7.2), Inches(3.28), Inches(5.3), Inches(0.35),
         "THE CAVEAT THAT COMES WITH IT", size=11.5, color=PALE, bold=True)
add_text(s, Inches(7.2), Inches(3.70), Inches(5.35), Inches(1.6),
         "Disruption at this locus is strongly serotype-structured — Asia1 is nearly "
         "invariant (0.03) while SAT2 is completely diverged (1.00). The pathogen-level "
         "verdict pools lineages that behave very differently. Before this locus is "
         "interpreted biologically, D.H1 must be run per serotype.",
         size=12.5, color=WHITE, line_spacing=1.28)

add_text(s, Inches(0.55), Inches(6.15), Inches(12.2), Inches(0.75),
         "Read plainly: at the one locus in this study with a decisive signal, G4 sites are "
         "MORE labile than matched controls, not less.\nThat contradicts D.H1. It is a real "
         "measurement, and reporting it as anything else would be the only actual failure here.",
         size=13.5, color=DARK, bold=True, line_spacing=1.35)
footer(s, "Serotype breakdown recorded in docs/revision_log.md; the superseded verdict is kept in the ledger, not overwritten", 10)

# ══════════════════════════════════════ 10 · what it outputs ═══════
s = add_slide(); set_bg(s)
header(s, "so what does it output today?", "A closed gate, and the reason for it")
add_text(s, Inches(0.55), Inches(1.42), Inches(12.2), Inches(0.5),
         "Across all eleven pathogens, the testing ledger currently holds no SUPPORTED D.H1 "
         "verdict. So:", size=15, color=DARK)

# Bar heights hug their text; a fixed 1.45in rule ran well past the
# two-line block and read as an unfinished list.
for title, body, colour, y, h in [
    ("Produced", "G4 Reference Atlas per pathogen · QC report · reference-anchored alignment · "
     "recombination screen · ML tree and time-scaled root · variant table · per-locus "
     "disruption rates · GC control gate · D.H1 verdict · report card · full provenance",
     TEAL, Inches(1.95), Inches(0.62)),
    ("NOT produced", "G4-EWS surveillance score · CUSUM / EWMA control limits · detection "
     "calls · any statement that a pathogen is or is not changing at a G4 site",
     RED, Inches(2.95), Inches(0.62)),
]:
    add_rect(s, Inches(0.55), y, Pt(4), h, colour)
    add_text(s, Inches(0.95), y, Inches(2.1), Inches(0.4), title, size=17, color=colour, bold=True)
    add_text(s, Inches(3.15), y, Inches(9.6), Inches(0.9), body, size=13.5, color=DARK,
             line_spacing=1.3)

add_text(s, Inches(0.55), Inches(4.05), Inches(12.2), Inches(0.35),
         "What the tool says when you ask it to score anyway:", size=13.5, color=NAVY, bold=True)
add_rect(s, Inches(0.55), Inches(4.45), Inches(12.25), Inches(1.32), LGREY,
         line=RGBColor(0xCE, 0xD8, 0xDE))
# Verbatim from `g4watch score --pathogen fmdv2026`, which exits 3.
add_text(s, Inches(0.85), Inches(4.62), Inches(11.7), Inches(1.05),
         "FMDV2026: Stage 5 scoring is BLOCKED — at least one locus shows a real,\n"
         "GC-adjustment-surviving effect that runs AGAINST D.H1: it is MORE disrupted than its\n"
         "matched control, not less. This is evidence against the hypothesis, not an absence of\n"
         "evidence for it, and it is reported separately for that reason.",
         size=12.5, color=NAVY, font="Consolas", line_spacing=1.22)

add_text(s, Inches(0.55), Inches(5.95), Inches(12.2), Inches(0.9),
         "This is the designed behaviour, not a broken run. A score computed from a hypothesis "
         "its own data does not support would be a\nnumber with nothing behind it — and it "
         "would look exactly like a number with something behind it.",
         size=13.5, color=DARK, line_spacing=1.32)
footer(s, "Quoted verbatim from `g4watch score --pathogen fmdv2026`, which exits 3 rather than emitting a score", 11)

# ══════════════════════════════════════ 11 · limits ════════════════
s = add_slide(); set_bg(s)
header(s, "limits", "What this does not yet establish")
bullets(s, Inches(0.55), Inches(1.50), Inches(6.0), Inches(5.2), [
    (0, "Predictions, not structures.", {"bold": True, "color": NAVY, "space_after": 3}),
    (0, "Every locus here is computational. No G4 in this Atlas has been confirmed in vitro "
        "by CD spectroscopy or NMR.", {"no_marker": True, "size": 13.5}),
    (0, "One test, one corpus.", {"bold": True, "color": NAVY, "space_after": 3}),
    (0, "Only FMDV2026 and its serotype-O subset have cleared the floor for enough loci to "
        "test. Eight pathogens have no D.H1 verdict at all.", {"no_marker": True, "size": 13.5}),
    (0, "Lineage pooling.", {"bold": True, "color": NAVY, "space_after": 3}),
    (0, "The pathogen-level verdict treats each corpus as one population. G4-004 shows that "
        "can hide large between-lineage differences.", {"no_marker": True, "size": 13.5}),
], size=15, space_after=13)
bullets(s, Inches(6.9), Inches(1.50), Inches(5.9), Inches(5.2), [
    (0, "Sampling is opportunistic.", {"bold": True, "color": NAVY, "space_after": 3}),
    (0, "GenBank is not a survey. Country and year distributions reflect who sequenced what, "
        "not where the virus was.", {"no_marker": True, "size": 13.5}),
    (0, "The control-limit calibration is short.", {"bold": True, "color": NAVY, "space_after": 3}),
    (0, "Baselines are annual and few, so limits carry explicit uncertainty intervals rather "
        "than being quoted as exact numbers.", {"no_marker": True, "size": 13.5}),
    (0, "No clinical or field validation.", {"bold": True, "color": NAVY, "space_after": 3}),
    (0, "Nothing here has been tested against outbreak outcomes. This is a hypothesis-testing "
        "framework, not a diagnostic.", {"no_marker": True, "size": 13.5}),
], size=15, space_after=13)
# The two columns filled only the top two thirds. Rather than padding them
# out, the space carries the thing a reader should leave with.
add_rect(s, Inches(0.55), Inches(5.35), Inches(12.25), Inches(1.5), LGREY)
add_text(s, Inches(0.95), Inches(5.58), Inches(11.5), Inches(1.1),
         "None of these limits is a reason to distrust the negative result.\n"
         "A framework that could not return a negative result would be the thing worth "
         "distrusting. What they do bound is how far it\ngeneralises: one hypothesis, one "
         "well-powered corpus, computational predictions, and no field validation yet.",
         size=13.5, color=NAVY, line_spacing=1.32)
footer(s, "Recorded in docs/validation_report.md and docs/revision_log.md", 12)

# ══════════════════════════════════════ 12 · methods table ═════════
# A reader who wants to argue with the result needs the parameters, not a
# description of them. Everything here is read from config/fmdv2026.yaml
# and `g4watch doctor`.
s = add_slide(); set_bg(s)
header(s, "methods", "Every parameter that changes the answer")

METHODS = [
    ("Stage 0 · G4 prediction", "G4Hunter", "window 25 · |score| ≥ 1.2 · flank 100 nt"),
    ("", "regex motif", "enabled · concordance at ≥ 80% overlap"),
    ("", "not used", "G4RNA screener (Python 2 only) · pqsfinder (DNA-scoped)"),
    ("Stage 1 · QC", "thresholds", "completeness ≥ 0.90 · N ≤ 0.05 · year-precision date"),
    ("Stage 1 · Alignment", "MAFFT v7.526", "--auto --keeplength --addfragments · gap cap 0.50"),
    ("Stage 1.5 · Recombination", "PhiPack Φ", "α = 0.05 · mandatory · tier: standard"),
    ("Stage 2 · Phylogeny", "IQ-TREE 2.3.6", "GTR+F+I+G4 · 1,000 UFBoot · seed 20250823"),
    ("", "TreeTime 0.11.4", "least-squares reroot · clock filter 3.0 IQD"),
    ("", "ancestral states", "ape::ace on the DIVERGENCE tree, not the timetree"),
    ("Stage 4 · Controls", "per locus", "5 controls · length ±10% · GC ±5% · 50 nt buffer"),
    ("", "matched on", "the locus's own compartment (5'UTR / CDS / 3'UTR)"),
    ("Stage 4.5 · D.H1", "tests", "Fisher exact + pooled GC-adjusted logistic"),
    ("", "correction", "Benjamini–Hochberg · α = 0.05 · decision: GC-adjusted alone"),
    ("", "analysis set", "loci in ≥ 20 genomes — pre-specified, score-blind"),
]
y = Inches(1.40)
row_h = Inches(0.375)
for i, (stage, item, value) in enumerate(METHODS):
    if i % 2 == 0:
        add_rect(s, Inches(0.55), y, Inches(12.25), row_h, LGREY)
    if stage:
        add_text(s, Inches(0.72), y + Inches(0.085), Inches(2.85), Inches(0.3), stage,
                 size=11.5, color=NAVY, bold=True)
    add_text(s, Inches(3.70), y + Inches(0.085), Inches(2.25), Inches(0.3), item,
             size=11.5, color=GREY)
    add_text(s, Inches(6.10), y + Inches(0.085), Inches(6.6), Inches(0.3), value,
             size=11.5, color=DARK, font="Consolas")
    y += row_h

add_text(s, Inches(0.55), Inches(6.85), Inches(12.2), Inches(0.35),
         "Two of these moved after data were seen, and both are recorded as post-hoc: the "
         "SC operating point (R-11) and the GC-adjusted decision rule (R-12).",
         size=11.5, color=AMBER, italic=True)
footer(s, "config/fmdv2026.yaml and `g4watch doctor`. Changing any of these changes the result, so none is buried in nextflow.config", 13)

# ══════════════════════════════════════ 12 · reproducibility ═══════
s = add_slide(); set_bg(s)
header(s, "reproducibility", "A run is described by three things")
for i, (num, title, body) in enumerate([
    ("1", "The commit", "Every stage records the git commit it ran under. Analysis code and "
     "configuration are versioned together."),
    ("2", "The config", "One YAML per pathogen: reference accession, corpus paths, QC "
     "thresholds, model, seed, and every exclusion rule."),
    ("3", "The accession list", "The exact GenBank records fetched. The acquisition step is a "
     "separate entry point so an analysis run can never silently re-fetch and change its own inputs."),
]):
    x = Inches(0.55) + i * Inches(4.16)
    add_rect(s, x, Inches(1.50), Inches(3.9), Inches(2.0), LGREY)
    add_text(s, x + Inches(0.25), Inches(1.68), Inches(0.5), Inches(0.5), num,
             size=30, color=TEAL, bold=True)
    add_text(s, x + Inches(0.85), Inches(1.78), Inches(2.9), Inches(0.4), title,
             size=15.5, color=NAVY, bold=True)
    add_text(s, x + Inches(0.25), Inches(2.32), Inches(3.4), Inches(1.1), body,
             size=12, color=DARK, line_spacing=1.25)

add_text(s, Inches(0.55), Inches(3.85), Inches(12.2), Inches(0.35),
         "Run it:", size=14, color=NAVY, bold=True)
add_rect(s, Inches(0.55), Inches(4.25), Inches(12.25), Inches(1.15), LGREY,
         line=RGBColor(0xCE, 0xD8, 0xDE))
add_text(s, Inches(0.85), Inches(4.45), Inches(11.7), Inches(0.85),
         "nextflow run workflow/main.nf --pathogen fmdv2026 -profile docker\n"
         "g4watch gate-status --pathogen fmdv2026        # never blocked; reports why",
         size=15, color=NAVY, font="Consolas", line_spacing=1.35)

for i, (v, lab) in enumerate([
    (f"{N_TESTS:,}", "automated tests passing"),
    (f"{N_PATHOGENS}", "pathogens provisioned"),
    (f"{N_ALIGNED:,}", "genomes aligned"),
    ("append-only", "testing ledger"),
]):
    stat(s, Inches(0.55) + i * Inches(3.15), Inches(5.65), Inches(3.0), v, lab)
footer(s, "Profiles: docker · singularity · conda_free · slurm · awsbatch", 14)

# ══════════════════════════════════════ 13 · close ═════════════════
s = add_slide(); set_bg(s, NAVY)
add_rect(s, Inches(0.9), Inches(2.05), Inches(1.5), Pt(4), TEAL)
add_text(s, Inches(0.9), Inches(2.45), Inches(11.4), Inches(2.2),
         "The framework works. The hypothesis,\nso far, does not hold.",
         size=36, color=WHITE, bold=True, line_spacing=1.25)
add_text(s, Inches(0.9), Inches(4.35), Inches(11.0), Inches(1.8),
         "Eleven pathogens, 6,485 genomes and an end-to-end pipeline that runs, records what it "
         "did, and refuses to score when its\nown gate is closed. At the one locus with a "
         "decisive signal, G4 sites are more labile than matched controls — the opposite\n"
         "of what D.H1 predicts. The next step is to test that per serotype, and to get more "
         "pathogens past the data floor.",
         size=15.5, color=PALE, line_spacing=1.4)
add_text(s, Inches(0.9), Inches(6.45), Inches(11.0), Inches(0.4),
         "ICAR–NIVEDI  ·  G4-WATCH  ·  research framework, not a validated diagnostic",
         size=11.5, color=RGBColor(0x7E, 0xA3, 0xB8))

prs.save(DECK)
print(f"wrote {DECK} — {len(prs.slides.__iter__.__self__._sldIdLst)} slides")
