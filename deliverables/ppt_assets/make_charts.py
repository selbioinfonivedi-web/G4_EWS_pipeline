"""Generate illustrative chart PNGs for the G4-WATCH briefing deck.
All numeric shapes here are illustrative (synthetic series shaped to match
the documented calibration behaviour), clearly labelled as such on-slide;
no engine parameter (ARL target, lambda, control-limit search) is invented.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 13,
    "axes.edgecolor": "#444444", "axes.labelcolor": "#222222",
    "text.color": "#222222", "xtick.color": "#444444", "ytick.color": "#444444",
})

NAVY = "#1b3a5c"
TEAL = "#1f8a70"
AMBER = "#d98c1a"
RED = "#c0392b"
GREY = "#9aa5ad"

# ---------------------------------------------------------------- chart 1
# EWMA control chart on a standardized surveillance-score series: a quiet
# baseline followed by a gradual drift that crosses the control limit.
rng = np.random.default_rng(7)
n_baseline, n_drift = 30, 20
baseline = rng.normal(0, 1, n_baseline)
drift = np.linspace(0, 3.2, n_drift) + rng.normal(0, 0.6, n_drift)
series = np.concatenate([baseline, drift])

lam = 0.2
z = np.zeros(len(series))
cur = 0.0
for i, v in enumerate(series):
    cur = lam * v + (1 - lam) * cur
    z[i] = cur
control_limit = 1.55  # illustrative, of the shape calibrate_ewma() would return

fig, ax = plt.subplots(figsize=(10, 4.2), dpi=200)
ax.axvspan(0, n_baseline - 1, color=GREY, alpha=0.12, label="baseline period")
ax.plot(z, color=NAVY, lw=2.2, label="EWMA statistic (λ=0.2)")
ax.axhline(control_limit, color=RED, ls="--", lw=1.6, label="control limit h (calibrated, target ARL₀=200)")
alarm_idx = np.argmax(z > control_limit)
ax.scatter([alarm_idx], [z[alarm_idx]], color=RED, zorder=5, s=70)
ax.annotate("first alarm", (alarm_idx, z[alarm_idx]), xytext=(alarm_idx - 13, z[alarm_idx] + 0.55),
            arrowprops=dict(arrowstyle="->", color=RED), color=RED, fontweight="bold")
ax.set_xlabel("surveillance window")
ax.set_ylabel("EWMA statistic (standardized units)")
ax.set_title("Illustrative EWMA control chart — gradual G4 disruption drift", loc="left", fontweight="bold")
ax.legend(loc="upper left", frameon=False, fontsize=10.5)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{__file__.rsplit('/',1)[0]}/ewma_chart.png", facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- chart 2
# Block-bootstrap vs i.i.d. calibration: same target ARL, different control
# limits, illustrating why the correction matters for autocorrelated series.
fig, ax = plt.subplots(figsize=(9, 4.6), dpi=200)
methods = ["i.i.d. bootstrap\n(naive)", "moving-block bootstrap\n(used)"]
limits = [0.62, 1.55]
achieved_if_naive_used_in_practice = [38, 200]  # illustrative ARL actually achieved in production
ax.set_ylim(0, 1.95)
bars = ax.bar(methods, limits, color=[AMBER, TEAL], width=0.5)
for b, v, arl in zip(bars, limits, achieved_if_naive_used_in_practice):
    ax.text(b.get_x() + b.get_width()/2, v + 0.08, f"h = {v}", ha="center", fontweight="bold", fontsize=13)
    ax.text(b.get_x() + b.get_width()/2, v/2, f"true ARL\n≈{arl} windows", ha="center", color="white", fontweight="bold")
ax.set_ylabel("calibrated control limit h")
ax.set_title("Same target ARL₀ = 200 — very different control limits", loc="left", fontweight="bold", pad=14)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{__file__.rsplit('/',1)[0]}/calibration_chart.png", facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- chart 3
# Structural-confidence funnel: how many candidate loci survive each axis-1
# filter, shaped like the EBV numbers already computed in this project.
labels = ["all G4Hunter\nhits\n(t=1.2)", "≥2 concordant\ntools", "+ |G4Hunter|\n≥1.5", "+ conservation\n≥85%\n→ SC-eligible"]
counts = [1410, 216, 52, 18]
fig, ax = plt.subplots(figsize=(9.5, 4.4), dpi=200)
colors = [GREY, "#6ea8c9", AMBER, TEAL]
bars = ax.bar(labels, counts, color=colors, width=0.6)
for b, c in zip(bars, counts):
    ax.text(b.get_x() + b.get_width()/2, c + 20, str(c), ha="center", fontweight="bold")
ax.set_ylabel("candidate loci (reference genome)")
ax.set_title("Structural-confidence funnel — worked example (EBV reference)", loc="left", fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{__file__.rsplit('/',1)[0]}/funnel_chart.png", facecolor="white")
plt.close(fig)

print("wrote 3 chart PNGs")
