# -*- coding: utf-8 -*-
"""Builds the G4-WATCH technical briefing deck (python-pptx)."""
import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.oxml.ns import qn

HERE = os.path.dirname(os.path.abspath(__file__))
# The finished deck belongs beside the other deliverables, not in the
# asset folder next to this script's scratch partials. It used to be
# written to HERE, so `python build_deck.py` left the published copy in
# deliverables/ untouched and the two silently diverged — the deck shipped
# for weeks describing thresholds the code no longer used.
OUT = os.path.dirname(HERE)
DECK = os.path.join(OUT, "G4-WATCH_technical_briefing.pptx")

NAVY = RGBColor(0x1B, 0x3A, 0x5C)
TEAL = RGBColor(0x1F, 0x8A, 0x70)
AMBER = RGBColor(0xD9, 0x8C, 0x1A)
RED = RGBColor(0xC0, 0x39, 0x2B)
GREY = RGBColor(0x5A, 0x66, 0x70)
LGREY = RGBColor(0xEE, 0xF1, 0xF3)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
DARK = RGBColor(0x22, 0x2A, 0x30)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]
SW, SH = prs.slide_width, prs.slide_height


def add_slide():
    return prs.slides.add_slide(BLANK)


def set_bg(slide, color=WHITE):
    bg = slide.background
    bg.fill.solid()
    bg.fill.fore_color.rgb = color


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
             align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, font="Calibri", line_spacing=1.0,
             wrap=True):
    box = slide.shapes.add_textbox(x, y, w, h)
    tf = box.text_frame
    tf.word_wrap = wrap
    tf.vertical_anchor = anchor
    tf.margin_left = 0
    tf.margin_right = 0
    tf.margin_top = 0
    tf.margin_bottom = 0
    lines = text.split("\n")
    for i, line in enumerate(lines):
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


def bullets(slide, x, y, w, h, items, size=16, color=DARK, bold_head=False,
            space_after=8, font="Calibri", anchor=MSO_ANCHOR.TOP, line_spacing=1.05):
    """items: list of (level, text) or (level, text, {overrides})"""
    box = slide.shapes.add_textbox(x, y, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = 0
    tf.margin_right = 0
    tf.margin_top = 0
    tf.margin_bottom = 0
    first = True
    for item in items:
        level, text = item[0], item[1]
        override = item[2] if len(item) > 2 else {}
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.level = level
        p.space_after = Pt(override.get("space_after", space_after))
        p.line_spacing = override.get("line_spacing", line_spacing)
        marker = "" if override.get("no_marker") else {0: "▪  ", 1: "–  ", 2: "·  "}.get(level, "")
        run = p.add_run()
        run.text = marker + text
        run.font.size = Pt(override.get("size", size - level * 1.5))
        run.font.bold = override.get("bold", bold_head and level == 0)
        run.font.italic = override.get("italic", False)
        run.font.name = override.get("font", font)
        run.font.color.rgb = override.get("color", color)
    return box


def header(slide, kicker, title, kicker_color=TEAL):
    add_rect(slide, 0, 0, SW, Inches(1.05), NAVY)
    add_text(slide, Inches(0.55), Inches(0.10), Inches(11), Inches(0.32), kicker.upper(),
              size=13, color=RGBColor(0x9F, 0xC7, 0xD8), bold=True)
    add_text(slide, Inches(0.55), Inches(0.38), Inches(12.2), Inches(0.65), title,
              size=27, color=WHITE, bold=True)
    add_rect(slide, 0, Inches(1.05), SW, Pt(3), kicker_color)


def footer(slide, note, page):
    add_text(slide, Inches(0.55), Inches(7.16), Inches(10.5), Inches(0.3), note,
              size=10, color=GREY, italic=True)
    add_text(slide, Inches(12.5), Inches(7.16), Inches(0.6), Inches(0.3), str(page),
              size=10, color=GREY, align=PP_ALIGN.RIGHT)


def content_area():
    return Inches(0.55), Inches(1.35), Inches(12.25), Inches(5.65)


def mono_box(slide, x, y, w, h, code_text, size=14.5, color=NAVY, bg=LGREY, line_spacing=1.15):
    add_rect(slide, x, y, w, h, bg, line=RGBColor(0xCE, 0xD8, 0xDE))
    add_text(slide, x + Inches(0.18), y + Inches(0.12), w - Inches(0.36), h - Inches(0.24),
              code_text, size=size, color=color, font="Consolas", line_spacing=line_spacing)


def diagram_box(slide, x, y, w, h, title, subtitle="", fill=NAVY, text_color=WHITE, title_size=13, sub_size=10.5):
    shp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, x, y, w, h)
    shp.adjustments[0] = 0.10
    shp.fill.solid()
    shp.fill.fore_color.rgb = fill
    shp.line.color.rgb = WHITE
    shp.line.width = Pt(0.75)
    shp.shadow.inherit = False
    tf = shp.text_frame
    tf.word_wrap = True
    tf.margin_left = Pt(4); tf.margin_right = Pt(4); tf.margin_top = Pt(3); tf.margin_bottom = Pt(3)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = title; r.font.size = Pt(title_size); r.font.bold = True
    r.font.color.rgb = text_color; r.font.name = "Calibri"
    if subtitle:
        p2 = tf.add_paragraph()
        p2.alignment = PP_ALIGN.CENTER
        r2 = p2.add_run(); r2.text = subtitle; r2.font.size = Pt(sub_size)
        r2.font.color.rgb = text_color; r2.font.name = "Calibri"
    return shp


def arrow(slide, x1, y1, x2, y2, color=GREY, width=1.75):
    conn = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, x1, y1, x2, y2)
    conn.line.color.rgb = color
    conn.line.width = Pt(width)
    line = conn.line._get_or_add_ln()
    tail = line.makeelement(qn('a:tailEnd'), {'type': 'triangle', 'w': 'med', 'len': 'med'})
    line.append(tail)
    return conn


# ======================================================================
# SLIDE 1 — Title
# ======================================================================
s = add_slide(); set_bg(s, NAVY)
add_rect(s, 0, Inches(4.55), SW, Pt(2), TEAL)
add_text(s, Inches(0.9), Inches(2.15), Inches(11.5), Inches(0.5), "ICAR – NIVEDI", size=18,
         color=RGBColor(0x9F, 0xC7, 0xD8), bold=True)
add_text(s, Inches(0.9), Inches(2.6), Inches(11.5), Inches(1.3), "G4-WATCH", size=58, color=WHITE, bold=True)
add_text(s, Inches(0.9), Inches(3.65), Inches(11.5), Inches(0.7),
         "G-Quadruplex Genomic Early-Warning Surveillance for Livestock Viruses",
         size=22, color=RGBColor(0xD8, 0xE4, 0xE9))
add_text(s, Inches(0.9), Inches(4.75), Inches(11.5), Inches(1.6),
         "How it works · the formulas and statistical tools · the workflow · the report card ·\n"
         "how baseline, threshold and weights are derived · what the Early-Warning Score means",
         size=15, color=RGBColor(0xB8, 0xC9, 0xD1), italic=True, line_spacing=1.3)
add_text(s, Inches(0.9), Inches(6.85), Inches(8), Inches(0.4),
         "Technical briefing — all formulas and thresholds cited from the implemented codebase",
         size=11.5, color=RGBColor(0x8A, 0x9F, 0xAA))

# ======================================================================
# SLIDE 2 — Why this exists
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Motivation", "Why a G-quadruplex-based early-warning layer")
x, y, w, h = content_area()
add_text(s, x, y, w, Inches(0.9),
         "G-quadruplexes (G4s) are four-stranded DNA/RNA secondary structures formed by G-rich sequence. "
         "In several viruses they sit over functionally important genome regions and their disruption or loss "
         "has been linked to altered replication, translation, or immune evasion.",
         size=16, color=DARK, line_spacing=1.25)
bullets(s, x, y + Inches(1.1), w, Inches(4.0), [
    (0, "The surveillance question:", {"bold": True, "size": 17, "color": NAVY}),
    (1, "Are conserved G4 motifs in a pathogen's genome being disrupted, lost, or newly gained "
        "at a rate that differs from background — and is that rate accelerating in specific lineages "
        "or geographies?"),
    (0, "Why this needs its own statistical framework, not a simple mutation counter:", {"bold": True, "size": 17, "color": NAVY}),
    (1, "G4 motifs are GC-rich by construction — any raw disruption signal is confounded with GC content "
        "and must be adjusted for it before it means anything."),
    (1, "A handful of tip sequences from one outbreak clade are not independent observations; naive "
        "per-genome counting manufactures false significance."),
    (1, "Not every computationally predicted G4 is real. The framework must separate "
        "\"a motif that scores well\" from \"a motif with genuine evidence of forming and mattering\"."),
    (0, "G4-WATCH is built to answer the question only after these three problems are handled — "
        "never before.", {"bold": True, "color": TEAL, "size": 16.5}),
], size=16)
footer(s, "Framework: G4-WATCH Concept Paper v2 / architecture Sections 5–15", 2)

# ======================================================================
# SLIDE 3 — Workflow diagram (clean)
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Pipeline", "End-to-end workflow — Stage 0 through Stage 6")
x0 = Inches(0.5)
top = Inches(1.55)
bw, bh = Inches(1.28), Inches(1.05)
gap = Inches(0.10)
stage_defs = [
    ("Stage 0", "G4 Reference\nAtlas", NAVY),
    ("Stage 1", "Sequence QC", NAVY),
    ("Stage 1.5", "Recombination\nscreen (PHI)", NAVY),
    ("Stage 2", "Alignment +\nphylogeny", NAVY),
    ("Stage 3", "Variant /\ntip-state calls", NAVY),
    ("Stage 4", "Appendix C\ndata floor", AMBER),
    ("Stage 4.5", "D.H1 gate\n(GC-adjusted)", RED),
    ("Stage 5", "G4-EWS score +\nD.H2–D.H4", TEAL),
    ("Stage 6", "Report card +\ndashboard", NAVY),
]
n = len(stage_defs)
total_w = n * bw + (n - 1) * gap
xs = Emu(int((SW - total_w) / 2))
boxes = []
for i, (tag, label, color) in enumerate(stage_defs):
    bx = Emu(int(xs) + i * int(bw + gap))
    b = diagram_box(s, bx, top, bw, bh, tag, label, fill=color, title_size=11, sub_size=8.5)
    boxes.append((bx, b))
    if i < n - 1:
        arrow(s, Emu(int(bx) + int(bw)), Emu(int(top) + int(bh/2)),
              Emu(int(bx) + int(bw) + int(gap)), Emu(int(top) + int(bh/2)), color=GREY, width=2.25)

# gate branch annotation under stage 4.5
gy = top + bh + Inches(0.35)
add_rect(s, xs, gy, total_w, Inches(1.55), LGREY, line=RGBColor(0xCE, 0xD8, 0xDE))
bullets(s, xs + Inches(0.25), gy + Inches(0.12), total_w - Inches(0.5), Inches(1.35), [
    (0, "The D.H1 gate is the single point of control: Stage 5/6 scoring is refused (exit code 3) "
        "unless D.H1 returns SUPPORTED for this pathogen and its config sets operational_mode: true.", {"size": 14.5}),
    (0, "NOT_SUPPORTED and SIGNAL_EXPLAINED_BY_GC are not failures of the pipeline — they are documented "
        "scientific findings written to the append-only testing ledger.", {"size": 14.5}),
], size=14.5)

add_text(s, xs, gy + Inches(1.65), total_w, Inches(0.35),
         "Stages 1–3 build the evidence base • Stage 4/4.5 decide IF scoring is scientifically valid • Stage 5/6 produce the score, ONLY if permitted",
         size=12.5, color=GREY, italic=True, align=PP_ALIGN.CENTER)
footer(s, "14 Nextflow processes implement these 9 stages end to end (profiles: standard / slurm / awsbatch / docker / singularity)", 3)

# ======================================================================
# SLIDE 4 — Tools per stage
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Tools", "What runs at each stage")
x, y, w, h = content_area()
rows = [
    ("Stage", "Task", "Tool(s)"),
    ("0", "G4 motif scan, two-axis confidence classification", "G4Hunter (custom), pattern-motif scanner, Biopython"),
    ("1", "Sequence QC (completeness, N-content, date precision)", "custom QC (g4watch.qc), Biopython"),
    ("1.5", "Recombination screening (mandatory, gates the floor)", "PhiPack (Phi test), vendored"),
    ("2", "Multiple alignment; ML phylogeny; time-scaled tree", "MAFFT 7.5; IQ-TREE2 2.3.6; TreeTime 0.11.4"),
    ("2/4", "Ancestral state reconstruction on the rooted tree", "R + ape 5.7.1 (Rscript)"),
    ("3", "Variant calling vs. reference; per-genome tip states", "custom (g4watch.variants), reference-pinned alignment"),
    ("4/4.5", "Minimum-data floor; GC-confound-adjusted gate", "SciPy (Fisher exact); pooled logistic regression + LRT"),
    ("5", "Seven-term metrics; G4-EWS-core / integrated score; CUSUM/EWMA", "NumPy/SciPy; custom scoring package"),
    ("6", "Report card, gate-status dashboard, calibration check", "custom reporting; FastAPI + SSE web dashboard"),
]
table_shape = s.shapes.add_table(len(rows), 3, x, y, w, Inches(4.9)).table
table_shape.columns[0].width = Inches(1.1)
table_shape.columns[1].width = Inches(5.9)
table_shape.columns[2].width = Inches(5.25)
for r, row in enumerate(rows):
    for c, val in enumerate(row):
        cell = table_shape.cell(r, c)
        cell.text = val
        cell.margin_left = Pt(6); cell.margin_right = Pt(6); cell.margin_top = Pt(3); cell.margin_bottom = Pt(3)
        cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        para = cell.text_frame.paragraphs[0]
        para.font.size = Pt(13 if r else 13.5)
        para.font.bold = (r == 0)
        para.font.name = "Calibri"
        para.font.color.rgb = WHITE if r == 0 else DARK
        cell.fill.solid()
        cell.fill.fore_color.rgb = NAVY if r == 0 else (LGREY if r % 2 == 0 else WHITE)
footer(s, "Container images: g4watch/core, alignment, phylogenetics, selection, statistics, g4prediction, variants, web-backend/db/proxy, acquisition", 4)

# ======================================================================
# SLIDE 5 — Stage 0: Atlas & two-axis confidence
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stage 0", "The G4 Reference Atlas — two independent confidence axes")
x, y, w, h = content_area()
bullets(s, x, y, Inches(5.9), Inches(3.6), [
    (0, "Axis 1 — Structural Confidence", {"bold": True, "color": NAVY, "size": 16}),
    (1, "Is there real evidence this G4 forms at all? The ONLY axis that gates "
        "eligibility for scoring.", {"size": 13.5}),
    (0, "Axis 2 — Functional Context", {"bold": True, "color": NAVY, "size": 16}),
    (1, "Does it sit in an annotated functional region? Reported for interpretation "
        "— never gates.", {"size": 13.5}),
    (0, "Kept strictly independent by design", {"bold": True, "color": TEAL, "size": 15}),
    (1, "so a well-annotated but weakly-supported motif can never borrow credibility "
        "from its location, and a strong motif in an unannotated region is never "
        "discarded for that reason alone.", {"size": 13.5}),
], size=14, space_after=6)
# confidence ladder
ladder = [
    ("EC", "Experimentally Confirmed", TEAL),
    ("BC", "Biophysically Confirmed", TEAL),
    ("SC", "Strong Computational Candidate", RGBColor(0x6E, 0xA8, 0xC9)),
    ("MC", "Moderate Computational Candidate", AMBER),
    ("WC", "Weak Computational Candidate", RGBColor(0xC9, 0x9A, 0x4E)),
    ("AA", "Algorithm Artefact (gap / low-quality region)", RED),
]
ly = Inches(1.55)
lx = Inches(6.7)
lw, lh = Inches(6.1), Inches(0.66)
for i, (tag, label, color) in enumerate(ladder):
    yy = Emu(int(ly) + i * int(lh + Inches(0.06)))
    shp = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, lx, yy, lw, lh)
    shp.adjustments[0] = 0.18
    shp.fill.solid(); shp.fill.fore_color.rgb = color
    shp.line.fill.background(); shp.shadow.inherit = False
    tf = shp.text_frame; tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = Pt(10)
    p = tf.paragraphs[0]
    r = p.add_run(); r.text = f"{tag}  —  {label}"; r.font.size = Pt(14); r.font.bold = True
    r.font.color.rgb = WHITE; r.font.name = "Calibri"
add_rect(s, x, Inches(5.15), Inches(5.9), Inches(1.85), LGREY, line=RGBColor(0xCE, 0xD8, 0xDE))
add_text(s, x + Inches(0.18), Inches(5.28), Inches(5.55), Inches(1.6),
         "SC rule (relaxed 2026-09; a judgment call, NOT a ROC fit \u2014 see R-11):\n"
         "  ≥ 1 prediction tool\n"
         "  AND |G4Hunter score| ≥ 1.2\n"
         "  AND phylogenetic conservation ≥ 75%\n"
         "Was ≥ 2 tools / 1.5 / 85%: 0% sensitivity to confirmed viral G4s.\n"
         "G4Hunter's sign is strand, not quality — SC/MC compare |score|.",
         size=12.5, font="Consolas", color=NAVY, line_spacing=1.2)
footer(s, "g4watch/atlas/confidence.py — structural_confidence() / functional_context()", 5)

prs.save(os.path.join(HERE, "_partial_1.pptx"))
print("part 1 done:", len(prs.slides.__iter__.__self__._sldIdLst))

# ======================================================================
# SLIDE 6 — Appendix C minimum-data floor
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stage 4", "The Appendix C minimum-data floor — 9 unconditional checks")
x, y, w, h = content_area()
add_text(s, x, y, w, Inches(0.55),
         "Runs BEFORE the D.H1 gate, on every locus, and never short-circuits — every failing check is "
         "reported at once so a pathogen is told INSUFFICIENT_DATA for a specific, complete reason.",
         size=15, color=DARK, italic=True)
rows = [
    ("Check", "Threshold"),
    ("Sequences in window", "≥ 30"),
    ("Sequences per lineage (smallest)", "≥ 20"),
    ("Distinct timepoints", "≥ 3"),
    ("Metadata completeness (date/host/country)", "≥ 90%"),
    ("Informative clades — locus group", "≥ 3"),
    ("Informative clades — control group", "≥ 3"),
    ("Matched control region found", "required"),
    ("Alignment/QC pass fraction", "≥ 50%"),
    ("Recombination screen completed", "required"),
]
t = s.shapes.add_table(len(rows), 2, x, y + Inches(0.7), Inches(7.6), Inches(4.6)).table
t.columns[0].width = Inches(5.5); t.columns[1].width = Inches(2.1)
for r, row in enumerate(rows):
    for c, val in enumerate(row):
        cell = t.cell(r, c); cell.text = val
        cell.margin_left = Pt(6); cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        para = cell.text_frame.paragraphs[0]
        para.font.size = Pt(13.5); para.font.bold = (r == 0)
        para.font.color.rgb = WHITE if r == 0 else DARK
        cell.fill.solid(); cell.fill.fore_color.rgb = NAVY if r == 0 else (LGREY if r % 2 == 0 else WHITE)
bullets(s, Inches(8.35), y + Inches(0.7), Inches(4.4), Inches(4.6), [
    (0, "Only the first three numbers (30 / 20 / 3) are stated in the architecture text as "
        "\"the floor\"; the remaining six are this project's own documented judgment call — "
        "pipeline-completeness preconditions later stages already depend on implicitly.", {"size": 13.5}),
    (0, "INSUFFICIENT_DATA ≠ NOT_SUPPORTED", {"bold": True, "color": RED, "size": 15}),
    (1, "the floor firing means \"cannot test yet\"; D.H1 returning NOT_SUPPORTED means "
        "\"tested, and found no effect\". The dashboard must never conflate the two.", {"size": 13.5}),
], size=13.5)
footer(s, "g4watch/validation/minimum_data_gate.py", 6)

# ======================================================================
# SLIDE 7 — D.H1 gate mechanics
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stage 4.5", "D.H1 — the gating hypothesis test", kicker_color=RED)
x, y, w, h = content_area()
mono_box(s, x, y, Inches(6.0), Inches(1.15),
         'D.H1: "G4 Atlas loci show significantly\n'
         'lower disruption than matched controls,\n'
         'surviving GC-content adjustment."',
         size=13.5, color=NAVY)
bullets(s, x, y + Inches(1.35), Inches(6.0), Inches(4.1), [
    (0, "Step 1 — raw test (always computed, both rules)", {"bold": True, "size": 15, "color": NAVY}),
    (1, "Fisher's exact test on locus-vs-control disruption, per locus (2×2 table).", {"size": 13.5}),
    (0, "Step 2 — GC-confound adjustment", {"bold": True, "size": 15, "color": NAVY}),
    (1, "ONE pooled logistic model across every locus and control simultaneously:", {"size": 13.5}),
    (2, "disruption ~ GC + locus_1_dummy + ... + locus_k_dummy",
        {"size": 12.5, "font": "Consolas", "color": NAVY}),
    (1, "Each locus's dummy tested by likelihood-ratio test (full model vs. that dummy "
        "removed), then Benjamini–Hochberg FDR-corrected across all loci.", {"size": 13.5}),
    (0, "Why pooled, not per-pair:", {"bold": True, "size": 14.5}),
    (1, "a single locus-control pair has only 2 distinct GC values — GC and locus identity "
        "are perfectly collinear and unidentifiable alone.", {"size": 13.5}),
    (0, "Step 3 — direction, checked explicitly", {"bold": True, "size": 15, "color": RED}),
    (1, "Both tests above are TWO-SIDED — they ask \u201cis this locus different?\u201d, never "
        "\u201cis it more conserved?\u201d. D.H1 is directional, so locus_rate < control_rate is "
        "required separately. Without it a locus MORE disrupted than its control opens the gate.",
        {"size": 13}),
], size=14, space_after=6)
# verdict table
vt_rows = [
    ("Verdict", "Meaning"),
    ("SUPPORTED", "Survives GC adjustment AND is LESS disrupted than its control"),
    ("SIGNAL_OPPOSITE_DIRECTION", "Real GC-adjusted effect, but MORE disrupted — against D.H1"),
    ("SIGNAL_EXPLAINED_BY_GC", "Raw effect present, explained away by GC (methods-null)"),
    ("NOT_SUPPORTED", "No significant difference at all (biological null)"),
]
t = s.shapes.add_table(len(vt_rows), 2, Inches(6.85), y + Inches(0.1), Inches(5.9), Inches(2.4)).table
t.columns[0].width = Inches(2.5); t.columns[1].width = Inches(3.4)
for r, row in enumerate(vt_rows):
    for c, val in enumerate(row):
        cell = t.cell(r, c); cell.text = val
        cell.margin_left = Pt(6); cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        para = cell.text_frame.paragraphs[0]
        para.font.size = Pt(12.5); para.font.bold = (r == 0)
        para.font.color.rgb = WHITE if r == 0 else DARK
        cell.fill.solid(); cell.fill.fore_color.rgb = RED if r == 0 else (LGREY if r % 2 == 0 else WHITE)
add_rect(s, Inches(6.85), y + Inches(2.35), Inches(5.9), Inches(2.9), LGREY, line=RGBColor(0xCE,0xD8,0xDE))
bullets(s, Inches(7.1), y + Inches(2.5), Inches(5.4), Inches(2.6), [
    (0, "Two selectable decision rules", {"bold": True, "color": NAVY, "size": 15.5}),
    (1, "conjunction (default): raw Fisher AND GC-adjusted test must both be significant.", {"size": 13.5}),
    (1, "gc_adjusted: the GC-adjusted FDR test alone decides; raw Fisher is still computed and reported.", {"size": 13.5}),
    (0, "Gate is the ONLY switch", {"bold": True, "color": TEAL, "size": 15}),
    (1, "No code downstream may treat a pathogen as scoreable without an actual SUPPORTED Dh1GateResult "
        "passed in explicitly.", {"size": 13.5}),
], size=13.5)
footer(s, "g4watch/validation/dh1_gate.py, gc_confound_gate.py", 7)

# ======================================================================
# SLIDE 8 — D.H2 / D.H3 / D.H4
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stages 5–6", "D.H2, D.H3, D.H4 — the three questions D.H1 unlocks")
x, y, w, h = content_area()
cards = [
    ("D.H2", "Does G4 information add predictive value\nbeyond conventional genomic-epidemiology signals?",
     "Nested model comparison M1→M4 (next slide). Answered conservatively: requires BOTH a significant "
     "likelihood-ratio gain of M4 over M2 AND M4 beating M2 on held-out AUC.", TEAL),
    ("D.H3", "Do G4 state transitions (present→disrupted/gained)\noccur more often in expanding clades?",
     "Fisher's exact test on (transition / no-transition) × (expanding / not-expanding), clade-level — "
     "NOT tip-level, which would treat one outbreak's 200 sequences as 200 independent events.", AMBER),
    ("D.H4", "Same question as D.H2, from the other direction:\ndoes the G4 signal hold once conventional\nsignals are controlled for?",
     "Answered by the same nested LRT + held-out-AUC test as D.H2 — reported together, "
     "since they are the two directions of one comparison.", RGBColor(0x6E, 0xA8, 0xC9)),
]
cw = Inches(3.95)
for i, (tag, q, note, color) in enumerate(cards):
    cx = x + Emu(int(i * (cw + Inches(0.2))))
    diagram_box(s, cx, y, cw, Inches(1.35), tag, "", fill=color, title_size=20)
    add_text(s, cx, y + Inches(1.5), cw, Inches(1.2), q, size=14, color=NAVY, bold=True, line_spacing=1.15)
    add_text(s, cx, y + Inches(2.75), cw, Inches(2.5), note, size=12.5, color=DARK, line_spacing=1.2)
add_text(s, x, Inches(6.55), w, Inches(0.4),
         "Clade trajectory (D.H3) is a sampling-share estimate from tip collection dates — a claim about the sequenced record, not directly about transmission incidence.",
         size=12, color=GREY, italic=True)
footer(s, "g4watch/validation/dh3_test.py, phylo/clade_growth.py, validation/model_comparison.py", 8)

# ======================================================================
# SLIDE 9 — the seven G.2 terms
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stage 5", "The seven G.2 surveillance terms")
x, y, w, h = content_area()
rows = [
    ("Term", "Definition", "Formula (per surveillance window)"),
    ("ΔG4C", "Change in mean G4 conservation vs. rolling baseline", "mean(conservation)ₜ − baseline mean"),
    ("G4D", "Mean weighted disruption frequency across loci", "disrupted genomes ÷ assessed genomes"),
    ("G4G", "Proportion of genomes carrying a novel G4 gain", "genomes with ≥1 novel G4Hunter hit ÷ N"),
    ("ΔG4MB", "Change in G4-associated mutation burden vs. baseline", "Σ G4-locus SNPs ÷ Σ all SNPs, Δ vs. baseline"),
    ("LF", "Fastest-growing lineage frequency change", "max over lineages of (freqₜ − freqₜ₋₁)"),
    ("GE", "Shannon entropy of the geography of G4-changed genomes", "−Σ pᵢ·ln(pᵢ), countries of disrupted/gained only"),
    ("TA", "Temporal acceleration: recent vs. historical substitution rate", "rate(recent) ÷ rate(history), clipped at 10"),
]
t = s.shapes.add_table(len(rows), 3, x, y, w, Inches(5.3)).table
t.columns[0].width = Inches(1.3); t.columns[1].width = Inches(5.3); t.columns[2].width = Inches(5.65)
for r, row in enumerate(rows):
    for c, val in enumerate(row):
        cell = t.cell(r, c); cell.text = val
        cell.margin_left = Pt(6); cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        para = cell.text_frame.paragraphs[0]
        para.font.size = Pt(12.5 if r else 13)
        para.font.bold = (r == 0)
        para.font.name = "Consolas" if c == 2 and r else "Calibri"
        para.font.color.rgb = WHITE if r == 0 else DARK
        cell.fill.solid(); cell.fill.fore_color.rgb = TEAL if r == 0 else (LGREY if r % 2 == 0 else WHITE)
footer(s, "g4watch/metrics/surveillance_metrics.py — ΔG4C/G4D/G4G/ΔG4MB feed model M3; LF/GE/TA are added for M4 only", 9)

prs.save(os.path.join(HERE, "_partial_2.pptx"))
print("part 2 done:", len(prs.slides._sldIdLst))

# ======================================================================
# SLIDE 10 — Baseline & normalization
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stage 5", "How the baseline is set — and why z-scoring is mandatory")
x, y, w, h = content_area()
add_text(s, x, y, w, Inches(0.85),
         "The seven raw terms live on incomparable scales — a proportion, an unbounded entropy, a ratio. "
         "Summing them with equal weight=1 would not be \"equal weight\" at all: it silently gives each "
         "term influence proportional to its own raw numeric range.",
         size=15.5, color=DARK, line_spacing=1.25)
mono_box(s, x, y + Inches(1.0), Inches(7.6), Inches(1.3),
         "z = (raw_value − baseline_mean) / baseline_stdev\n\n"
         "baseline_mean, baseline_stdev = mean, stdev(baseline_values)", size=16, color=NAVY)
bullets(s, x, y + Inches(2.5), Inches(7.6), Inches(2.9), [
    (0, "Baseline values", {"bold": True, "color": NAVY, "size": 16}),
    (1, "the SAME metric's own values across the pre-surveillance-window baseline period, "
        "for THIS pathogen — never a cross-pathogen or literature default."),
    (0, "Hard guards, not conventions", {"bold": True, "color": TEAL, "size": 16}),
    (1, "fewer than 2 baseline values → refuses (cannot estimate a standard deviation)."),
    (1, "zero-variance baseline → refuses, rather than divide by zero or return a meaningless 0."),
    (1, "NormalizedMetric is a distinct type: every scoring function requires it as input, so a raw, "
        "un-normalized float cannot be passed in and silently distort a score."),
], size=15)
add_rect(s, Inches(8.35), y + Inches(1.0), Inches(3.9), Inches(2.6), LGREY, line=RGBColor(0xCE,0xD8,0xDE))
add_text(s, Inches(8.55), y + Inches(1.15), Inches(3.5), Inches(2.3),
         "\"Baseline period\" is a design commitment, not a knob:\n\n"
         "the pre-surveillance window recorded once per pathogen "
         "(architecture Section 9.4), so a later window can never "
         "be scored against a baseline chosen to make it look calmer.",
         size=13, color=NAVY, italic=True, line_spacing=1.25)
footer(s, "g4watch/metrics/normalization.py — z_against_baseline()", 10)

# ======================================================================
# SLIDE 11 — G4-EWS-core (M3)
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stage 5", "G4-EWS-core (model M3) — the G4-only score")
x, y, w, h = content_area()
mono_box(s, x, y, w, Inches(1.1),
         "G4-EWS-core = w1·ΔG4C  +  w2·G4D  +  w3·G4G  +  w4·ΔG4MB*\n"
         "(every term a NormalizedMetric — already z-scored against its own baseline)", size=17, color=NAVY)
bullets(s, x, y + Inches(1.35), w, Inches(2.0), [
    (0, "Exactly 4 terms, by construction", {"bold": True, "size": 16, "color": TEAL}),
    (1, "This fixes a real defect found in Revision 1: its \"G4-only\" formula had quietly mixed in "
        "3 conventional epidemiology terms, which made the whole point of the M1–M4 comparison — "
        "\"does G4 information add value beyond conventional signals\" — not well-posed. "
        "M3 is defined as exactly this function; there is no other \"G4-only\" formula anywhere else "
        "in the codebase for it to contradict."),
    (0, "ΔG4MB*  =  orthogonalized excess mutation burden", {"bold": True, "size": 16, "color": NAVY}),
    (1, "residualized so it captures mutation burden IN EXCESS of what genome-wide mutation rate "
        "already predicts — not raw mutation count, which would double-count background evolution."),
], size=15)
footer(s, "g4watch/scoring/g4_ews_core.py, metrics/mutation_burden_residual.py", 11)

# ======================================================================
# SLIDE 12 — Integrated score (M4) & model family
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stage 5", "Integrated Surveillance Score (M4) and the nested model family")
x, y, w, h = content_area()
mono_box(s, x, y, Inches(7.6), Inches(1.0),
         "Integrated Score  =  G4-EWS-core  +  w5·LF + w6·GE + w7·TA", size=17, color=NAVY)
t_rows = [
    ("Model", "Terms", "Question it answers"),
    ("M1", "LF only (null)", "baseline signal from lineage frequency alone"),
    ("M2", "LF + GE + TA", "conventional genomic-epidemiology signal, no G4"),
    ("M3", "ΔG4C + G4D + G4G + ΔG4MB*", "G4-EWS-core — G4 signal alone"),
    ("M4", "M3 + M2 (all 7 terms)", "does G4 add value beyond conventional signals?"),
]
t = s.shapes.add_table(len(t_rows), 3, x, y + Inches(1.25), Inches(7.6), Inches(2.6)).table
t.columns[0].width = Inches(0.9); t.columns[1].width = Inches(3.0); t.columns[2].width = Inches(3.7)
for r, row in enumerate(t_rows):
    for c, val in enumerate(row):
        cell = t.cell(r, c); cell.text = val
        cell.margin_left = Pt(6); cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        para = cell.text_frame.paragraphs[0]
        para.font.size = Pt(13); para.font.bold = (r == 0)
        para.font.color.rgb = WHITE if r == 0 else DARK
        cell.fill.solid(); cell.fill.fore_color.rgb = NAVY if r == 0 else (LGREY if r % 2 == 0 else WHITE)
add_rect(s, Inches(8.35), y, Inches(3.9), Inches(4.3), LGREY, line=RGBColor(0xCE,0xD8,0xDE))
bullets(s, Inches(8.55), y + Inches(0.15), Inches(3.5), Inches(4.0), [
    (0, "Decision rule (D.H2/D.H4)", {"bold": True, "color": NAVY, "size": 15}),
    (1, "M4 significantly beats M2 by likelihood-ratio test", {"size": 13}),
    (1, "AND M4 beats M2 on held-out AUC", {"size": 13}),
    (0, "Both required — deliberately conservative:", {"bold": True, "color": RED, "size": 14.5}),
    (1, "a likelihood gain that does not survive to held-out discrimination is overfitting, "
        "not evidence G4 information helps.", {"size": 13}),
], size=13.5)
footer(s, "g4watch/scoring/integrated_score.py, validation/model_comparison.py — g4_adds_value()", 12)

prs.save(os.path.join(HERE, "_partial_3.pptx"))
print("part 3 done:", len(prs.slides._sldIdLst))

# ======================================================================
# SLIDE 13 — How the WEIGHTS are calculated
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Weights", "How w1…w4 are actually fitted — regularized, never guessed")
x, y, w, h = content_area()
bullets(s, x, y, Inches(7.3), Inches(2.2), [
    (0, "Regularized logistic MLE, no intercept", {"bold": True, "color": NAVY, "size": 16.5}),
    (1, "no intercept because G4-EWS-core is a pure weighted sum with no constant term — "
        "fitting one would fit a model the score itself cannot represent."),
    (0, "L2 penalty is the DEFAULT, not an option", {"bold": True, "color": RED, "size": 16.5}),
    (1, "the four core terms are correlated by construction (disruption and conservation-change "
        "measure overlapping aspects of the same locus history). An unregularized fit at this "
        "project's sample sizes gives large, unstable, sign-flipping coefficients that LOOK like "
        "findings. Setting strength=0 requires passing allow_unregularized=True explicitly."),
], size=15)
mono_box(s, x, y + Inches(2.35), Inches(7.3), Inches(0.55),
         "objective(w) = −log-likelihood(w)  +  strength · ‖w‖²", size=15.5, color=NAVY)
bullets(s, x, y + Inches(3.05), Inches(7.3), Inches(1.9), [
    (0, "Mandatory prior-sensitivity report", {"bold": True, "color": TEAL, "size": 16}),
    (1, "the fit is re-run at 5 penalty strengths spanning two orders of magnitude "
        "(0.1, 0.5, 1.0, 5.0, 20.0). If the fitted weight DIRECTION barely moves across that "
        "range, the signal is in the data; if it swings, it is in the prior."),
    (1, "prior_dominated = True (cosine deviation > 0.2, or any sign flip) → the weights "
        "must NOT be reported as a finding, by rule, not by judgment call at report time."),
], size=15)
add_rect(s, Inches(8.55), y, Inches(3.7), Inches(4.6), LGREY, line=RGBColor(0xCE,0xD8,0xDE))
add_text(s, Inches(8.75), y + Inches(0.15), Inches(3.3), Inches(0.4), "Design rationale, verbatim:",
         size=13.5, bold=True, color=NAVY)
add_text(s, Inches(8.75), y + Inches(0.6), Inches(3.3), Inches(3.8),
         "\"The optimiser is plain MLE with a penalty term, deliberately simple and readable, "
         "because the point of this module is auditability rather than speed.\"",
         size=13.5, color=DARK, italic=True, line_spacing=1.3)
footer(s, "g4watch/scoring/weight_fitting.py — fit_core_weights(), PriorSensitivityReport", 13)

# ======================================================================
# SLIDE 14 — From score to alarm: CUSUM & EWMA
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Detection", "Turning a score series into an alarm — CUSUM & EWMA")
x, y, w, h = content_area()
mono_box(s, x, y, Inches(6.0), Inches(0.85),
         "CUSUM:  Sₜ = max(0,  Sₜ₋₁ + zₜ − k)     alarm if Sₜ > h", size=15.5, color=NAVY)
mono_box(s, x, y + Inches(1.0), Inches(6.0), Inches(0.85),
         "EWMA:   Eₜ = λ·zₜ + (1−λ)·Eₜ₋₁          alarm if Eₜ > h", size=15.5, color=NAVY)
bullets(s, x, y + Inches(2.05), Inches(6.0), Inches(3.2), [
    (0, "One-sided, upper only", {"bold": True, "size": 15, "color": NAVY}),
    (1, "the surveillance question is directional (an INCREASE in disruption); a two-sided "
        "chart would spend half its false-alarm budget on decreases nobody would act on.", {"size": 13.5}),
    (0, "Why both charts", {"bold": True, "size": 15, "color": NAVY}),
    (1, "CUSUM accumulates evidence — strongest against a sustained step shift. EWMA's decaying "
        "weight tracks a gradual drift — the more plausible shape for a G4 locus eroding "
        "over successive seasons.", {"size": 13.5}),
    (0, "λ = 0.2 (EWMA smoothing constant)", {"bold": True, "size": 15, "color": TEAL}),
    (1, "standard surveillance default balancing responsiveness against noise.", {"size": 13.5}),
], size=14)
add_text(s, Inches(6.5), y, Inches(6.2), Inches(3.4), "", size=1)
s.shapes.add_picture(os.path.join(HERE, "ewma_chart.png"), Inches(6.55), y - Inches(0.05), width=Inches(6.2))
footer(s, "g4watch/scoring/cusum.py, ewma.py — chart is an illustrative synthetic series, not a real result", 14)

# ======================================================================
# SLIDE 15 — Calibrating the threshold (control limit)
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Threshold", "How the alarm THRESHOLD is calibrated — target ARL₀ = 200")
x, y, w, h = content_area()
bullets(s, x, y, Inches(6.0), Inches(3.0), [
    (0, "The problem being corrected", {"bold": True, "color": RED, "size": 16}),
    (1, "surveillance score windows share sequences, share clades, share the same slowly-changing "
        "viral population — consecutive scores are NOT independent draws. Calibrating a control "
        "limit under an i.i.d. assumption sets it far too low; in production the alarm rate comes "
        "out many times the nominal false-positive rate. \"A system that cries wolf is worse "
        "than none.\""),
    (0, "The fix — moving-block bootstrap", {"bold": True, "color": TEAL, "size": 16}),
    (1, "resample the BASELINE series in contiguous blocks (preserving its short-range dependence), "
        "then bisect on the control limit h until the target average run length (ARL₀ = 200 quiet "
        "windows per false alarm) is achieved on 500 such resamples."),
], size=14.5)
add_text(s, x, y + Inches(3.15), Inches(6.0), Inches(1.4),
         "Block length: n^(1/3) rule of thumb, doubled when lag-1 autocorrelation exceeds 0.5; "
         "bounded to [2, n/4]. Needs ≥20 baseline observations — fewer would let the limit be "
         "dominated by the baseline's own sampling noise.",
         size=13, color=GREY, italic=True, line_spacing=1.25)
s.shapes.add_picture(os.path.join(HERE, "calibration_chart.png"), Inches(6.55), y + Inches(0.15), width=Inches(6.2))
footer(s, "g4watch/scoring/cusum.py — calibrate_cusum(), compare_iid_vs_block_calibration()", 15)

prs.save(os.path.join(HERE, "_partial_4.pptx"))
print("part 4 done:", len(prs.slides._sldIdLst))

# ======================================================================
# SLIDE 16 — Warning level ladder
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Detection", "From alarms to a warning level an operator can act on")
x, y, w, h = content_area()
bullets(s, x, y, Inches(6.0), Inches(1.6), [
    (0, "Raw level = consecutive alarm windows", {"bold": True, "size": 15.5, "color": NAVY}),
    (1, "1 → WATCH, 2 → ELEVATED, ≥3 → HIGH. Requiring persistence is the cheapest guard "
        "against acting on one transient in an autocorrelated series.", {"size": 13.5}),
], size=14)
levels = [
    ("HIGH", "≥ 3 consecutive alarm windows", RED),
    ("ELEVATED", "2 consecutive alarm windows", AMBER),
    ("WATCH", "1 alarm window", RGBColor(0xC9, 0x9A, 0x4E)),
    ("NONE", "no current alarm", TEAL),
    ("INSUFFICIENT_EVIDENCE", "gate closed, OR analysis underpowered", GREY),
]
ly = y + Inches(1.7)
for i, (tag, desc, color) in enumerate(levels):
    yy = ly + Emu(int(i * int(Inches(0.62))))
    diagram_box(s, x, yy, Inches(2.6), Inches(0.52), tag, "", fill=color, title_size=13)
    add_text(s, x + Inches(2.8), yy, Inches(3.3), Inches(0.52), desc, size=12.5, color=DARK, anchor=MSO_ANCHOR.MIDDLE)
add_rect(s, Inches(6.85), y, Inches(5.9), Inches(5.0), LGREY, line=RGBColor(0xCE,0xD8,0xDE))
bullets(s, Inches(7.1), y + Inches(0.15), Inches(5.4), Inches(4.7), [
    (0, "Two caps enforced before a level is ever shown", {"bold": True, "color": RED, "size": 16}),
    (1, "Underpowered analysis → capped at INSUFFICIENT_EVIDENCE, always, whatever the chart did. "
        "An alarm from a test that could not have detected the effect is not evidence of it.", {"size": 14}),
    (1, "Structural-confidence cap → EC/BC can reach HIGH; SC caps at ELEVATED; MC and WC cap at "
        "WATCH; AA can never warn at all (NONE).", {"size": 14}),
    (0, "gate_permitted = False (D.H1 not SUPPORTED)", {"bold": True, "color": NAVY, "size": 15.5}),
    (1, "→ INSUFFICIENT_EVIDENCE unconditionally, before either cap is even checked.", {"size": 14}),
    (0, "Every WarningAssessment states the binding cap by name", {"bold": True, "color": TEAL, "size": 15}),
    (1, "\"chart alone indicated ELEVATED; capped by structural confidence WC (cap: WATCH)\" — "
        "never a bare level with no way to see why.", {"size": 14}),
], size=14)
footer(s, "g4watch/warning/classifier.py — classify_warning(), CONFIDENCE_CAP", 16)

# ======================================================================
# SLIDE 17 — Structural confidence funnel (visualization)
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Visualization", "Worked example — how many candidates actually reach SC")
x, y, w, h = content_area()
s.shapes.add_picture(os.path.join(HERE, "funnel_chart.png"), x, y, width=Inches(7.6))
add_rect(s, Inches(8.35), y, Inches(3.9), Inches(4.6), LGREY, line=RGBColor(0xCE,0xD8,0xDE))
bullets(s, Inches(8.55), y + Inches(0.15), Inches(3.5), Inches(4.3), [
    (0, "Reading this funnel", {"bold": True, "color": NAVY, "size": 15.5}),
    (1, "Most raw G4Hunter hits (window=25, t=1.2) never reach SC eligibility — this is expected, "
        "not a defect: the SC bar is deliberately conservative.", {"size": 13.5}),
    (0, "Cross-pathogen honesty check", {"bold": True, "color": RED, "size": 15}),
    (1, "This funnel varies enormously by pathogen. FMDV's reference genome (53.6% GC) yields only "
        "4 raw G4Hunter hits genome-wide — it may genuinely lack canonical G4 structure.", {"size": 13.5}),
    (1, "Running this funnel before committing to a full pipeline run is now standard practice: "
        "it tells you in minutes whether a pathogen can ever reach a scoreable Atlas.", {"size": 13.5}),
], size=13.5)
footer(s, "Illustrative worked example (EBV reference genome, NC_007605.1); numbers from this project's own runs", 17)

# ======================================================================
# SLIDE 18 — Stage 6: the report card
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Stage 6", "The pathogen report card — one document per pathogen")
x, y, w, h = content_area()
sections = ["Identity", "Dataset", "Quality (QC)", "Genomic metrics (G.2 terms)",
            "D.H1–D.H4 verdicts", "Early-warning indicators", "Score", "Confidence",
            "Validation status", "Limitations", "Interpretation"]
cols = 3
cw2, ch2 = Inches(3.95), Inches(0.62)
for i, name in enumerate(sections):
    r, c = divmod(i, cols)
    bx = x + Emu(int(c * int(cw2 + Inches(0.1))))
    by = y + Emu(int(r * int(ch2 + Inches(0.12))))
    diagram_box(s, bx, by, cw2, ch2, name, "", fill=NAVY, title_size=13.5)
bullets(s, x, y + Inches(3.0), w, Inches(2.6), [
    (0, "Three rules enforced in code, not left to whoever assembles the report", {"bold": True, "color": RED, "size": 16}),
    (1, "Every section carries status reported / blocked / unavailable + a reason — a blank section "
        "and an absent one are never conflated."),
    (1, "Scored sections (score, warning level, EW indicators) appear ONLY when the D.H1 gate "
        "permits scoring; otherwise shown as blocked with the gate's own explanation."),
    (1, "A card built from synthetic/demonstration data is stamped non-authoritative at the top "
        "level — visible to every consumer without inspecting the underlying data."),
], size=14.5)
footer(s, "g4watch/reporting/report_card.py — build_report_card()", 18)

prs.save(os.path.join(HERE, "_partial_5.pptx"))
print("part 5 done:", len(prs.slides._sldIdLst))

# ======================================================================
# SLIDE 19 — Web dashboard / visualization surface
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Visualization", "The web dashboard — what a user actually sees")
x, y, w, h = content_area()
bullets(s, x, y, Inches(6.4), Inches(4.9), [
    (0, "FastAPI + SSE, session-authenticated", {"bold": True, "color": NAVY, "size": 16}),
    (1, "PBKDF2-HMAC-SHA256 (600,000 iterations) local auth; secure, signed session cookies; "
        "a route allowlist gates everything else behind login."),
    (0, "Per-pathogen dashboard page", {"bold": True, "color": NAVY, "size": 16}),
    (1, "gate status (never blocked — shows WHY scoring is or isn't permitted)"),
    (1, "the Atlas, both confidence axes shown independently"),
    (1, "recombination screen result, variant summary, molecular-clock diagnostics"),
    (1, "score trend, warning-level history, calibration status"),
    (0, "Live updates over Server-Sent Events", {"bold": True, "color": TEAL, "size": 16}),
    (1, "a long pipeline run streams stage-by-stage progress to the page rather than the user "
        "polling or watching a terminal."),
], size=14.5)
add_rect(s, Inches(7.1), y, Inches(5.65), Inches(4.9), LGREY, line=RGBColor(0xCE,0xD8,0xDE))
add_text(s, Inches(7.3), y + Inches(0.15), Inches(5.25), Inches(0.4), "Route surface (abridged):",
         size=14, bold=True, color=NAVY)
mono_box(s, Inches(7.3), y + Inches(0.6), Inches(5.4), Inches(4.1),
         "GET /\n"
         "  overview, all pathogens\n"
         "GET /pathogen/{name}\n"
         "  per-pathogen dashboard page\n"
         "GET /api/pathogen/{name}\n"
         "  full JSON payload\n"
         "GET /api/pathogen/{name}/gate\n"
         "  D.H1 gate status only\n"
         "POST /login   POST /logout   GET /whoami\n\n"
         "TLS + HSTS + CSP at the nginx proxy\n"
         "rate limits: 5 req/min login,\n"
         "             30 req/s general", size=13, color=NAVY, line_spacing=1.12)
footer(s, "web/backend/app.py, web/workstation/tracks.py, web/proxy/nginx.conf", 19)

# ======================================================================
# SLIDE 20 — Validation status (honest)
# ======================================================================
s = add_slide(); set_bg(s)
header(s, "Validation", "Current validation status — reported plainly, including the negative", kicker_color=RED)
x, y, w, h = content_area()
bullets(s, x, y, w, Inches(2.6), [
    (0, "The D.H1 gate was opening on evidence that contradicted D.H1", {"bold": True, "size": 16, "color": RED}),
    (1, "D.H1 claims G4 loci are LESS disrupted than matched controls. Both underlying tests — Fisher "
        "and the pooled likelihood-ratio test — are two-sided, so neither could tell 'more conserved' "
        "from 'more variable'. Direction was never checked."),
    (1, "On the FMDV 2026 corpus the ONLY locus returning SUPPORTED, FMDV2026-G4-004, was more "
        "disrupted than its control (0.775 vs 0.545). It carried the entire pathogen verdict, and that "
        "verdict was what set operational_mode true. Fixed: direction is now required explicitly, and "
        "a wrong-direction effect gets its own verdict (SIGNAL_OPPOSITE_DIRECTION) rather than being "
        "reported as no effect."),
    (1, "Two unit tests asserted the inverted behaviour, so the suite agreed with the code and both "
        "disagreed with the hypothesis. Corrected together."),
], size=14)
bullets(s, x, y + Inches(2.75), w, Inches(2.1), [
    (0, "The SC threshold was relaxed on cross-species evidence — a judgment call, not a ROC fit", {"bold": True, "size": 16, "color": AMBER}),
    (1, "Every confirmed viral G4 available (HIV-1 LTR and nef) classified as WC under the old "
        "≥2-tool / ≥1.5 / ≥85% rule — 0% sensitivity against its own ground truth. Now ≥1 tool, "
        "≥1.2, ≥75%. The confirmed set holds 3 loci from 1 virus with derived coordinates; "
        "calibration.py requires 30 across 4 families and refuses to report an operating point."),
    (0, "Both findings are the framework working, not failing", {"bold": True, "size": 16, "color": TEAL}),
    (1, "A gate that blocks scoring, and a gate defect caught before it authorised one, are documented "
        "results. The FMDV 2026 verdict is now SIGNAL_OPPOSITE_DIRECTION and scoring stays blocked."),
], size=14)
footer(s, "g4watch/validation/dh1_gate.py (direction); calibration.py — 30 positives / 4 families required, not met", 20)

# ======================================================================
# SLIDE 21 — Summary
# ======================================================================
s = add_slide(); set_bg(s, NAVY)
add_text(s, Inches(0.6), Inches(0.5), Inches(11), Inches(0.7), "Summary", size=34, color=WHITE, bold=True)
add_rect(s, Inches(0.6), Inches(1.25), Inches(3.2), Pt(3), TEAL)
bullets(s, Inches(0.6), Inches(1.6), Inches(12.1), Inches(5.4), [
    (0, "A gate, not a formula, is the core design decision", {"bold": True, "size": 18, "color": RGBColor(0x9F,0xC7,0xD8)}),
    (1, "D.H1 must return SUPPORTED before ANY score is produced for a pathogen — an unsupported "
        "pathogen gets a documented negative result, never a number that looks like one.", {"color": WHITE, "size": 15}),
    (0, "Every number in the score traces to a cited threshold or a fitted, auditable procedure",
     {"bold": True, "size": 18, "color": RGBColor(0x9F,0xC7,0xD8)}),
    (1, "z-scored against a pathogen's own baseline; weights regularized with a mandatory prior-"
        "sensitivity check; alarm thresholds calibrated to a stated target ARL by block-bootstrap.", {"color": WHITE, "size": 15}),
    (0, "Confidence caps make the warning level honest about its own evidence quality",
     {"bold": True, "size": 18, "color": RGBColor(0x9F,0xC7,0xD8)}),
    (1, "a weak computational candidate cannot raise a HIGH alarm however cleanly its chart fires.", {"color": WHITE, "size": 15}),
    (0, "Open, stated item: SC threshold recalibration against confirmed viral G4s",
     {"bold": True, "size": 18, "color": AMBER}),
    (1, "the operating point was relaxed on 3 confirmed loci from 1 virus. A defensible threshold "
        "needs 30 across 4 families with stated coordinates — the next required validation step, "
        "not yet performed.", {"color": WHITE, "size": 15}),
    (0, "No pathogen currently has an open gate", {"bold": True, "size": 18, "color": AMBER}),
    (1, "FMDV: INSUFFICIENT_DATA. FMDV 2026: SIGNAL_OPPOSITE_DIRECTION. Every surveillance score in "
        "this deck is machinery described, not a result claimed.", {"color": WHITE, "size": 15}),
], size=15, line_spacing=1.15, space_after=10)
footer_box = add_text(s, Inches(0.6), Inches(7.1), Inches(11), Inches(0.3),
         "G4-WATCH — ICAR-NIVEDI  ·  technical briefing  ·  all figures traced to the implemented codebase",
         size=11, color=RGBColor(0x8A,0x9F,0xAA), italic=True)

prs.save(DECK)
print("FINAL slide count:", len(prs.slides._sldIdLst))
print("wrote", DECK)
