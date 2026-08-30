"""Native G4Hunter implementation (Bedrat, Amina & Mergny, 2016 — Nucleic
Acids Research, "Prediction of putative G-quadruplexes by G4Hunter...").

Why native rather than a subprocess wrapper: G4Hunter is a short, fully
published, deterministic sliding-window scoring algorithm (no trained model,
no hidden parameters) — a different category from a tool like BEAST2 or a
trained ML classifier, where reimplementation would risk unverifiable
divergence from the validated original. The algorithm here was read
line-by-line from the original author's reference script (vendored at
`vendor/g4hunter_reference/G4Hunter.py`, fetched from
github.com/AnimaTardeb/G4Hunter, Python 2, GPLv3) and reimplemented for
Python 3 with an equivalent public API. The tests in
`tests/unit/g4prediction/test_g4hunter.py` hand-verify scores against that
same reference logic.

One deliberate, documented deviation from the reference script: when
merging adjacent qualifying windows into one hit region, the reference
script (`WriteSeq`) re-slices the raw sequence and recomputes per-base
scores fresh from the start of the merged span — which can under-count the
first few bases of a merged region if they are the tail of a G/C run that
began before the span. This implementation instead reuses per-base scores
computed once over the FULL input sequence and slices into them, which is
more correct (a run's cap-at-4 score at a given position never depends on
where you started slicing) and avoids that edge case entirely.
"""

from __future__ import annotations

from dataclasses import dataclass

DEFAULT_WINDOW = 25
DEFAULT_THRESHOLD = 1.2
_RUN_CAP = 4


@dataclass(frozen=True)
class G4HunterHit:
    """One merged region where the windowed G4Hunter score met the
    threshold. Coordinates are 0-based, `end` exclusive (Python slice
    convention) — callers mapping to 1-based genome coordinates must add 1
    to `start`."""

    start: int
    end: int
    score: float
    strand: str  # "+" (G-rich, forms on this strand) or "-" (C-rich, forms on the complementary strand)


def per_base_scores(sequence: str) -> list[int]:
    """One score per base: +min(run_length, 4) inside a run of consecutive
    G's, -min(run_length, 4) inside a run of consecutive C's, 0 otherwise.
    Case-insensitive; any non-G/C character (including T/U/A/N/ambiguity
    codes) scores 0, matching the reference implementation exactly (its
    U/T/A-specific branch is present only as dead/commented-out code)."""
    seq = sequence.upper()
    n = len(seq)
    scores = [0] * n
    i = 0
    while i < n:
        base = seq[i]
        if base == "G" or base == "C":
            j = i
            while j < n and seq[j] == base:
                j += 1
            run_len = j - i
            magnitude = min(run_len, _RUN_CAP)
            value = magnitude if base == "G" else -magnitude
            for k in range(i, j):
                scores[k] = value
            i = j
        else:
            i += 1
    return scores


def sliding_window_scores(sequence: str, window: int = DEFAULT_WINDOW) -> list[float]:
    """Mean per-base score in each window of length `window`, one entry per
    window start position (step 1). Returns an empty list if the sequence
    is shorter than `window`."""
    scores = per_base_scores(sequence)
    n = len(scores)
    if n < window:
        return []
    return [sum(scores[start : start + window]) / window for start in range(n - window + 1)]


def predict(
    sequence: str,
    window: int = DEFAULT_WINDOW,
    threshold: float = DEFAULT_THRESHOLD,
) -> list[G4HunterHit]:
    """Finds every maximal run of consecutive qualifying windows (windowed
    score >= threshold, or <= -threshold) and reports each as one merged
    hit region, scored by the mean per-base score across the FULL merged
    span (not the mean of the constituent window means — matches the
    reference tool's own merged-region reporting convention).

    Non-obvious, verified-against-the-reference-tool behavior: a merged
    region spans strictly MORE bases than any single qualifying window (a
    run of k consecutive qualifying window-starts covers k-1+window bases,
    versus one window's own `window` bases), so its own recomputed mean can
    legitimately fall BELOW the per-window threshold that triggered it —
    e.g. a 25nt-window hit spanning 3 consecutive qualifying window starts
    covers 27 bases, and if the 2 extra boundary bases are weakly-scoring,
    the reported region score can end up under 1.2 even at the default
    threshold. This was confirmed empirically (not assumed) against a real
    FMDV genome region during Sprint 2 and traced to this exact mechanism —
    it is a property of the region-averaging convention, not a scoring bug.
    Callers that need "was ANY window in this hit >= threshold" rather than
    "is the merged region's own mean >= threshold" should use
    `sliding_window_scores` directly instead of relying on `.score` here."""
    if threshold <= 0:
        raise ValueError("threshold must be positive; sign is handled internally")

    scores = per_base_scores(sequence)
    window_scores = sliding_window_scores(sequence, window)
    qualifying = [i for i, s in enumerate(window_scores) if s >= threshold or s <= -threshold]

    hits: list[G4HunterHit] = []
    idx = 0
    while idx < len(qualifying):
        run_start = qualifying[idx]
        run_end = run_start
        j = idx
        while j + 1 < len(qualifying) and qualifying[j + 1] == qualifying[j] + 1:
            j += 1
            run_end = qualifying[j]
        idx = j + 1

        span_start = run_start
        span_end = run_end + window  # exclusive; last merged window covers [run_end, run_end+window)
        region = scores[span_start:span_end]
        region_mean = sum(region) / len(region)
        strand = "+" if region_mean > 0 else "-"
        hits.append(
            G4HunterHit(start=span_start, end=span_end, score=round(region_mean, 3), strand=strand)
        )
    return hits
