"""Figures for the G4-WATCH overview deck.

Every number plotted here is read from a file this repository produced:
the Atlas TSVs under data/atlases/, the append-only testing ledger, and
the pathogen configs. Nothing is simulated, smoothed or shaped to look
like an expected result.

That distinction matters because a sibling script,
``ppt_assets/make_charts.py``, deliberately draws ILLUSTRATIVE series for
the technical briefing -- synthetic data shaped to show what a control
chart looks like. Those figures are labelled as such on their slides and
are not reused here. A deck about what the pipeline found may only plot
what it found.
"""
from __future__ import annotations

import collections
import csv
import glob
import os
from pathlib import Path

import matplotlib
import yaml

matplotlib.use("Agg")  # headless: these run in CI with no display
import matplotlib.pyplot as plt  # noqa: E402 - must follow matplotlib.use

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

NAVY, TEAL, AMBER, RED, GREY = "#1b3a5c", "#1f8a70", "#d98c1a", "#c0392b", "#9aa5ad"
INK = "#222a30"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.edgecolor": "#8a949c", "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": "#5a6670", "ytick.color": "#5a6670",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})


def _atlas_rows(path: Path) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _ledger_latest_run(pathogen: str) -> list[dict]:
    """The ledger is append-only, so the current state is the latest run.

    Reading every row would mix a superseded verdict with the one that
    replaced it -- exactly what happened to FMDV2026-G4-004, whose
    SUPPORTED row from before the directionality fix still sits in the
    file above its SIGNAL_OPPOSITE_DIRECTION replacement.
    """
    with open(ROOT / "data/atlases/testing_ledger.tsv", newline="") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["pathogen"] == pathogen]
    if not rows:
        return []
    latest = max(r["timestamp"] for r in rows)
    return [r for r in rows if r["timestamp"] == latest]


# ── figure 1: the Atlas, and how much of it can be scored ───────────
def figure_atlas() -> Path:
    """Total loci, and the SC-eligible subset that can actually be scored.

    A single stacked bar of totals is dominated by EBV's 1,410 loci and
    says nothing useful: the quantity that gates everything downstream is
    how many loci reach SC, and for most pathogens that is zero. Plotting
    the total alongside it shows that the shortfall is in the evidence
    tier, not in how many G4s were predicted.
    """
    data = []
    for path in sorted(glob.glob(str(ROOT / "data/atlases/G4_Reference_Atlas_v*.tsv"))):
        base = os.path.basename(path)
        if "preconservation" in base.lower():
            continue                      # a superseded build of the same Atlas
        rows = _atlas_rows(Path(path))
        counts = collections.Counter(r["structural_confidence"] for r in rows)
        data.append((base.split(".")[-2].upper(), len(rows), counts.get("SC", 0)))
    data.sort(key=lambda d: (-d[2], -d[1]))

    fig, ax = plt.subplots(figsize=(10.4, 4.6))
    idx = list(range(len(data)))
    ax.bar([i - 0.19 for i in idx], [d[1] for d in data], width=0.38, color=GREY,
           label="predicted loci")
    ax.bar([i + 0.19 for i in idx], [d[2] for d in data], width=0.38, color=TEAL,
           label="scoring-eligible (SC)")
    # Log scale: EBV's 1,410 against BTV's 0 is three orders of magnitude,
    # and a linear axis renders every bar but the first as a flat line.
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(0, 3000)
    ax.set_yticks([0, 1, 10, 100, 1000])
    ax.set_yticklabels(["0", "1", "10", "100", "1,000"])
    for i, (_, total, sc) in enumerate(data):
        ax.text(i - 0.19, total * 1.35 + 0.3, str(total), ha="center", fontsize=9, color="#5a6670")
        ax.text(i + 0.19, max(sc, 0) * 1.35 + 0.3, str(sc), ha="center", fontsize=9,
                color=TEAL, weight="bold")
    # FMDV2026 is wide enough to collide with its neighbour at this figure
    # width, so the labels tilt rather than the figure growing.
    ax.set_xticks(idx)
    ax.set_xticklabels([d[0] for d in data], rotation=28, ha="right")
    ax.set_ylabel("Atlas loci (log scale)")
    # Above the plot: inside it, the legend sits on top of EBV's bar.
    ax.legend(frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.02))
    ax.set_title("The Atlas, and how much of it can be scored", loc="left",
                 fontsize=13, color=NAVY, weight="bold", pad=34)
    # No caption baked into the figure: the slide that carries it says this
    # in its own body text, and two copies collide at slide scale.
    fig.tight_layout()
    out = HERE / "fig_atlas.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


# ── figure 2: the D.H1 test, locus rate against control rate ────────
def figure_dh1() -> Path:
    """Every locus that cleared the data floor, plotted against its controls.

    D.H1 predicts points BELOW the diagonal: G4 loci less disrupted than
    matched controls. The plot is the result -- the cloud sits on the
    diagonal, and the one locus far off it is above.
    """
    rows = [r for r in _ledger_latest_run("FMDV2026") if r["locus_disruption_rate"]]
    x = [float(r["control_disruption_rate"]) for r in rows]
    y = [float(r["locus_disruption_rate"]) for r in rows]
    verdicts = [r["verdict"] for r in rows]

    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.plot([0, 1], [0, 1], color="#c3ccd2", lw=1.2, zorder=1)
    ax.text(0.62, 0.575, "no difference", rotation=45, fontsize=9.5,
            color="#8a949c", ha="center", va="bottom")
    ax.fill_between([0, 1], [0, 1], [0, 0], color=TEAL, alpha=0.05, zorder=0)
    ax.text(0.72, 0.09, "D.H1 predicts\nloci fall here", fontsize=10, color=TEAL, ha="center")

    for xi, yi, v in zip(x, y, verdicts):
        if v == "SIGNAL_OPPOSITE_DIRECTION":
            ax.scatter([xi], [yi], s=150, color=RED, zorder=4, edgecolor="white", linewidth=1.5)
            ax.annotate("FMDV2026-G4-004\nq = 1.1e-09",
                        (xi, yi), textcoords="offset points", xytext=(-14, -42),
                        fontsize=10, color=RED, ha="right", weight="bold")
        elif v == "SIGNAL_EXPLAINED_BY_GC":
            ax.scatter([xi], [yi], s=62, color=AMBER, zorder=3, edgecolor="white", linewidth=1)
        else:
            ax.scatter([xi], [yi], s=48, color=NAVY, alpha=0.55, zorder=2,
                       edgecolor="white", linewidth=0.8)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Disruption rate, matched control regions")
    ax.set_ylabel("Disruption rate, G4 locus")
    ax.set_title("D.H1 on FMDV2026 — 18 loci past the data floor", loc="left",
                 fontsize=13, color=NAVY, weight="bold")
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=NAVY, alpha=.6, label="not supported (15)"),
        plt.Line2D([], [], marker="o", ls="", color=AMBER, label="explained by GC (2)"),
        plt.Line2D([], [], marker="o", ls="", color=RED, label="opposite direction (1)"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left", fontsize=10)
    fig.tight_layout()
    out = HERE / "fig_dh1.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


# ── figure 3: corpus scale ──────────────────────────────────────────
def figure_corpus() -> Path:
    data = []
    for path in sorted(glob.glob(str(ROOT / "config/*.yaml"))):
        cfg = yaml.safe_load(Path(path).read_text()) or {}
        if cfg.get("provisioned") is False:
            continue
        md = (cfg.get("corpus") or {}).get("metadata_tsv")
        if not md or not (ROOT / md).is_file():
            continue
        corpus = sum(1 for _ in open(ROOT / md)) - 1
        aligned_dir = (ROOT / md).parent / "aligned"
        hits = sorted(aligned_dir.glob("*aligned_to_ref.fasta"))
        aligned = sum(1 for line in open(hits[0]) if line.startswith(">")) if hits else 0
        data.append((cfg["pathogen"], corpus, aligned))
    data.sort(key=lambda d: -d[1])

    fig, ax = plt.subplots(figsize=(10, 4.4))
    labels = [d[0] for d in data]
    idx = range(len(data))
    ax.bar([i - 0.19 for i in idx], [d[1] for d in data], width=0.38,
           color=GREY, label="fetched")
    ax.bar([i + 0.19 for i in idx], [d[2] for d in data], width=0.38,
           color=NAVY, label="passed QC and aligned")
    ax.set_xticks(list(idx))
    ax.set_xticklabels(labels, rotation=28, ha="right")
    ax.set_ylabel("Genomes")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.set_title("Corpora built, and what survived QC", loc="left",
                 fontsize=13, color=NAVY, weight="bold")
    fig.tight_layout()
    out = HERE / "fig_corpus.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


if __name__ == "__main__":
    for fn in (figure_atlas, figure_dh1, figure_corpus):
        print(f"wrote {fn()}")
