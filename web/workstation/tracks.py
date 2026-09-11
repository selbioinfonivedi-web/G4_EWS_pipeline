"""Derived views: genome tracks, ordination, alignment slices, spectra.

Everything here is computed from artifacts already on disk -- the
alignment, the ML distance matrix, the ancestral-state reconstruction --
none of which the interface was previously reading. Where an artifact is
absent the function returns ``None`` and the view reports itself
unavailable rather than drawing something plausible.

Results are cached per pathogen: the alignment is 848 x 8,206 and the
distance matrix is 847 x 847, so recomputing them per request would make
the interface feel broken.
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

from g4watch.io.fasta import read_fasta

REPO_ROOT = Path(__file__).resolve().parents[2]
GAP = set("-.Nn")


def _paths(pathogen: str) -> dict[str, Path]:
    """Artifact locations, derived from the pathogen's own config.

    This used to build them from ``data/reference_genomes/<pathogen>/corpus``,
    which assumed every pathogen keeps its corpus in a directory named after
    itself. A second FMDV corpus in ``corpus_2026`` was therefore invisible
    to every track view even though the files were there: the payload
    reported the alignment MISSING and the tracks silently returned nothing.
    The config already records where the corpus lives, so ask it.
    """
    from g4watch.config import ConfigError, load_config

    try:
        config = load_config(pathogen)
        base = Path(config.corpus_metadata_tsv).parent
        stem = pathogen.lower()
    except (ConfigError, TypeError, AttributeError):
        base = REPO_ROOT / "data" / "reference_genomes" / pathogen / "corpus"
        stem = pathogen

    def _first(*candidates: Path) -> Path:
        """The first candidate that exists, else the first (for the error)."""
        for candidate in candidates:
            if candidate.is_file():
                return candidate
        return candidates[0]

    return {
        "alignment": _first(
            base / "aligned" / f"{stem}_qc_passed_aligned_to_ref.fasta",
            base / "aligned" / f"{pathogen}_qc_passed_aligned_to_ref.fasta",
        ),
        "mldist": _first(
            base / "phylogenetics" / f"{stem}_iqtree.mldist",
            base / "phylogenetics" / f"{pathogen}_iqtree.mldist",
        ),
        "ancestral": base / "phylogenetics" / "ancestral_serotype_reconstruction.tsv",
        "variants": _first(
            base / "variants" / f"{stem}_variants_all.tsv",
            base / "variants" / f"{pathogen}_variants_all.tsv",
        ),
        "g4_variants": _first(
            base / "variants" / f"{stem}_g4_variant_intersection.tsv",
            base / "variants" / f"{pathogen}_g4_variant_intersection.tsv",
        ),
        "recombination": base / "recombination" / "recombination_screen.log",
    }


@lru_cache(maxsize=2)
def _matrix(pathogen: str) -> tuple[np.ndarray, list[str]] | None:
    """The alignment as a uint8 matrix, so column statistics are vectorised."""
    path = _paths(pathogen)["alignment"]
    if not path.is_file():
        return None
    records = read_fasta(path)
    ids = list(records)
    arr = np.frombuffer("".join(records[i].upper() for i in ids).encode("ascii", "replace"), dtype=np.uint8)
    return arr.reshape(len(ids), -1), ids


# ── genome tracks ───────────────────────────────────────────────────
@lru_cache(maxsize=2)
def genome_tracks(pathogen: str, window: int = 40) -> dict | None:
    """Diversity, GC and gap fraction in fixed windows along the genome.

    Diversity is 1 - (frequency of the commonest base), computed per
    column over ungapped bases only, then averaged per window. That is a
    plain, assumption-free measure: it does not model substitution and
    makes no claim about selection. It exists so a reader can see whether
    the G4 loci sit anywhere unusual before any statistics are run.
    """
    got = _matrix(pathogen)
    if got is None:
        return None
    mat, _ids = got
    n_seq, n_col = mat.shape

    codes = {b: ord(b) for b in "ACGT"}
    counts = np.stack([(mat == c).sum(axis=0) for c in codes.values()])  # 4 x n_col
    covered = counts.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        major = counts.max(axis=0) / np.where(covered == 0, 1, covered)
    diversity = np.where(covered == 0, np.nan, 1.0 - major)
    gc_col = np.where(covered == 0, np.nan, (counts[1] + counts[2]) / np.where(covered == 0, 1, covered))
    gap_frac = 1.0 - covered / n_seq

    n_win = int(np.ceil(n_col / window))
    pad = n_win * window - n_col

    def bin_mean(v: np.ndarray) -> list[float | None]:
        padded = np.concatenate([v, np.full(pad, np.nan)])
        with np.errstate(invalid="ignore"):
            m = np.nanmean(padded.reshape(n_win, window), axis=1)
        return [None if np.isnan(x) else round(float(x), 5) for x in m]

    return {
        "window": window,
        "n_windows": n_win,
        "genome_length": n_col,
        "n_sequences": n_seq,
        "starts": [i * window + 1 for i in range(n_win)],
        "diversity": bin_mean(diversity),
        "gc": bin_mean(gc_col),
        "gaps": bin_mean(gap_frac),
    }


def locus_vs_background(pathogen: str, loci: list[dict]) -> list[dict]:
    """Mean diversity inside each locus against the genome-wide mean.

    Reported as a plain ratio with no p-value attached. The statistical
    test is D.H1's job, and it is gated; this is a descriptive number.
    """
    tracks = genome_tracks(pathogen)
    if tracks is None:
        return []
    div = np.array([np.nan if d is None else d for d in tracks["diversity"]], dtype=float)
    window = tracks["window"]
    background = float(np.nanmean(div))
    out = []
    for locus in loci:
        lo = max((locus["start"] - 1) // window, 0)
        hi = min((locus["end"] - 1) // window + 1, len(div))
        seg = div[lo:hi]
        value = float(np.nanmean(seg)) if seg.size and not np.all(np.isnan(seg)) else None
        out.append(
            {
                "id": locus["id"],
                "diversity": None if value is None else round(value, 5),
                "background": round(background, 5),
                "ratio": None if value is None or background == 0 else round(value / background, 3),
            }
        )
    return out


# ── ordination ──────────────────────────────────────────────────────
@lru_cache(maxsize=2)
def ordination(pathogen: str) -> dict | None:
    """Classical MDS of the IQ-TREE maximum-likelihood distance matrix.

    Double-centre, eigendecompose, keep the first two axes. The percentage
    of variance each axis explains is reported, because an MDS plot whose
    axes explain almost nothing is a Rorschach test, not a result.
    """
    path = _paths(pathogen)["mldist"]
    if not path.is_file():
        return None
    lines = path.read_text().split("\n")
    n = int(lines[0].strip())
    ids, rows = [], []
    for line in lines[1 : n + 1]:
        parts = line.split()
        if not parts:
            continue
        ids.append(parts[0])
        rows.append([float(x) for x in parts[1:]])
    d = np.array(rows, dtype=float)
    if d.shape[0] != d.shape[1]:
        return None

    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d**2) @ j
    vals, vecs = np.linalg.eigh(b)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pos = np.clip(vals, 0, None)
    coords = vecs[:, :2] * np.sqrt(pos[:2])
    total = pos.sum() or 1.0
    return {
        "ids": ids,
        "x": [round(float(v), 5) for v in coords[:, 0]],
        "y": [round(float(v), 5) for v in coords[:, 1]],
        "explained": [round(float(pos[0] / total * 100), 2), round(float(pos[1] / total * 100), 2)],
        "n": n,
    }


# ── ancestral states ────────────────────────────────────────────────
@lru_cache(maxsize=2)
def ancestral_states(pathogen: str) -> dict | None:
    path = _paths(pathogen)["ancestral"]
    if not path.is_file():
        return None
    states: dict[str, dict] = {}
    tally: Counter = Counter()
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            state = row.get("most_likely_state", "")
            probs = {k: float(v) for k, v in row.items() if k not in {"node_id", "most_likely_state"} and v}
            states[row["node_id"]] = {
                "state": state,
                "support": round(max(probs.values()), 4) if probs else None,
            }
            tally[state] += 1
    return {"nodes": states, "distribution": dict(tally.most_common()), "n": len(states)}


# ── alignment slice ─────────────────────────────────────────────────
def alignment_slice(pathogen: str, start: int, end: int, max_rows: int = 300) -> dict | None:
    """Columns of the alignment for a region, as difference-from-reference.

    Rows are returned as compact strings: '.' where the base matches the
    first (reference) row, the base letter where it differs, '-' for a
    gap. That is the conventional MSA rendering and it keeps the payload
    small enough to send.
    """
    got = _matrix(pathogen)
    if got is None:
        return None
    mat, ids = got
    start = max(1, start)
    end = min(end, mat.shape[1])
    if end < start:
        return None
    block = mat[:, start - 1 : end]
    ref = block[0]
    rows = []
    for i in range(min(len(ids), max_rows)):
        row = block[i]
        chars = np.where(row == ref, ord("."), row).astype(np.uint8)
        rows.append({"id": ids[i], "seq": chars.tobytes().decode("ascii", "replace")})
    variable = int((np.array([(block == c).sum(axis=0) for c in map(ord, "ACGT")]).astype(bool).sum(axis=0) > 1).sum())
    return {
        "start": start,
        "end": end,
        "reference": ids[0],
        "n_shown": len(rows),
        "n_total": len(ids),
        "variable_columns": variable,
        "rows": rows,
    }


# ── mutation spectrum ───────────────────────────────────────────────
@lru_cache(maxsize=2)
def mutation_spectrum(pathogen: str) -> dict | None:
    """Counts of each substitution type against the reference row.

    Descriptive only. A transition/transversion ratio near 0.5 would mean
    no mutational bias at all, which for a real RNA virus would be a sign
    something is wrong with the alignment rather than a discovery.
    """
    got = _matrix(pathogen)
    if got is None:
        return None
    mat, _ids = got
    ref = mat[0]
    bases = "ACGT"
    codes = {ord(b): b for b in bases}
    spectrum: Counter = Counter()
    for i in range(1, mat.shape[0]):
        row = mat[i]
        diff = row != ref
        for r, q in zip(ref[diff], row[diff]):
            if r in codes and q in codes:
                spectrum[f"{codes[r]}>{codes[q]}"] += 1
    transitions = sum(spectrum[k] for k in ("A>G", "G>A", "C>T", "T>C"))
    transversions = sum(v for k, v in spectrum.items() if k not in ("A>G", "G>A", "C>T", "T>C"))
    return {
        "spectrum": dict(sorted(spectrum.items(), key=lambda kv: -kv[1])),
        "transitions": transitions,
        "transversions": transversions,
        "ti_tv": round(transitions / transversions, 3) if transversions else None,
        "compared_to": "reference row of the alignment",
    }


# ── geography ───────────────────────────────────────────────────────
# Approximate country centroids, to two decimal places. Enough to place a
# dot on a world map; never used for any distance or spatial statistic.
CENTROIDS = {
    "India": (22.35, 78.67),
    "Pakistan": (30.38, 69.35),
    "United Kingdom": (54.00, -2.00),
    "Viet Nam": (16.00, 106.00),
    "Thailand": (15.00, 101.00),
    "China": (35.00, 103.00),
    "Turkey": (39.00, 35.00),
    "Iran": (32.00, 53.00),
    "Kenya": (0.15, 37.90),
    "Nigeria": (9.08, 8.68),
    "Egypt": (26.00, 30.00),
    "South Africa": (-29.00, 24.00),
    "Zimbabwe": (-19.02, 29.15),
    "Botswana": (-22.33, 24.68),
    "Zambia": (-13.13, 27.85),
    "Uganda": (1.37, 32.29),
    "Ethiopia": (8.62, 39.62),
    "Tanzania": (-6.37, 34.89),
    "Israel": (31.50, 34.75),
    "Saudi Arabia": (24.00, 45.00),
    "Nepal": (28.17, 84.25),
    "Bangladesh": (24.00, 90.00),
    "Sri Lanka": (7.00, 81.00),
    "Myanmar": (22.00, 98.00),
    "Malaysia": (2.50, 112.50),
    "Indonesia": (-5.00, 120.00),
    "Philippines": (13.00, 122.00),
    "Republic of Korea": (37.00, 127.50),
    "South Korea": (37.00, 127.50),
    "Japan": (36.00, 138.00),
    "Russia": (60.00, 100.00),
    "Kazakhstan": (48.00, 68.00),
    "Mongolia": (46.00, 105.00),
    "Argentina": (-34.00, -64.00),
    "Brazil": (-10.00, -55.00),
    "Colombia": (4.00, -72.00),
    "Ecuador": (-2.00, -77.50),
    "Peru": (-10.00, -76.00),
    "Venezuela": (8.00, -66.00),
    "Bolivia": (-17.00, -65.00),
    "Paraguay": (-23.00, -58.00),
    "Uruguay": (-33.00, -56.00),
    "France": (46.00, 2.00),
    "Germany": (51.00, 9.00),
    "Italy": (42.83, 12.83),
    "Spain": (40.00, -4.00),
    "Greece": (39.00, 22.00),
    "Netherlands": (52.50, 5.75),
    "Bulgaria": (43.00, 25.00),
    "Albania": (41.00, 20.00),
    "Algeria": (28.00, 3.00),
    "Morocco": (32.00, -5.00),
    "Tunisia": (34.00, 9.00),
    "Libya": (25.00, 17.00),
    "Sudan": (15.00, 30.00),
    "Somalia": (10.00, 49.00),
    "Eritrea": (15.00, 39.00),
    "Cameroon": (6.00, 12.00),
    "Ghana": (8.00, -2.00),
    "Malawi": (-13.50, 34.00),
    "Mozambique": (-18.25, 35.00),
    "Namibia": (-22.00, 17.00),
    "Angola": (-12.50, 18.50),
    "Afghanistan": (33.00, 65.00),
    "Iraq": (33.00, 44.00),
    "Kuwait": (29.34, 47.66),
    "Bahrain": (26.00, 50.55),
    "United Arab Emirates": (24.00, 54.00),
    "Oman": (21.00, 57.00),
    "Yemen": (15.00, 48.00),
    "Taiwan": (23.50, 121.00),
    "Hong Kong": (22.25, 114.17),
    "Laos": (18.00, 105.00),
    "Cambodia": (13.00, 105.00),
    "Bhutan": (27.50, 90.50),
}


def geography(samples: list[dict]) -> dict:
    """Placeable and unplaceable counts, stated separately.

    A map that silently drops the records it cannot place misrepresents
    the corpus, so the unplaced count is returned alongside.
    """
    placed: dict[str, dict] = {}
    unplaced: Counter = Counter()
    for s in samples:
        country = s.get("c")
        if country in CENTROIDS:
            lat, lon = CENTROIDS[country]
            entry = placed.setdefault(country, {"lat": lat, "lon": lon, "n": 0, "lineages": {}})
            entry["n"] += 1
            entry["lineages"][s["l"]] = entry["lineages"].get(s["l"], 0) + 1
        elif country and country != "—":
            unplaced[country] += 1
    return {
        "places": placed,
        "n_placed": sum(p["n"] for p in placed.values()),
        "unplaced": dict(unplaced.most_common()),
        "n_unplaced": sum(unplaced.values()),
    }


# ── Stage 1.5 and Stage 3, surfaced ─────────────────────────────────
def recombination(pathogen: str) -> dict | None:
    """The PHI screen's verdict, or None if the screen has not run.

    Surfaced because Stage 1.5 gates the validity of everything
    phylogenetic downstream: reconstructing ancestral states across a
    recombinant alignment reconstructs a history that never happened. A
    reader looking at a tree deserves to see whether that check passed,
    rather than having to know it exists and go find the log.
    """
    path = _paths(pathogen)["recombination"]
    if not path.is_file():
        return None
    text = path.read_text()
    found: dict = {"log_path": str(path), "raw": text.strip()[-800:]}
    for pattern, key, cast in (
        (r"([0-9]+) sequences", "n_sequences", int),
        (r"([0-9]+) informative sites", "n_informative_sites", int),
        (r"p\s*=\s*([0-9.eE+-]+)", "p_value", float),
        (r"tier:\s*(\w+)", "tier", str),
        (r"significant:\s*(\w+)", "significant", lambda v: v.strip().lower() == "true"),
    ):
        match = re.search(pattern, text)
        if match:
            try:
                found[key] = cast(match.group(1))
            except (ValueError, TypeError):
                pass
    found["interpretation"] = (
        "Significant recombination detected. A single tree does not describe this corpus; "
        "phylogenetic results downstream are not trustworthy until it is partitioned."
        if found.get("significant")
        else "No significant recombination detected, so a single tree is a defensible "
             "description of this corpus."
    )
    return found


def variant_summary(pathogen: str) -> dict | None:
    """Counts from Stage 3, including the share falling inside Atlas loci.

    Deliberately a summary and not the table: the FMDV 2026 variant file is
    over a million rows, and shipping that to a browser would be a way of
    appearing to show the data while making it unreadable.
    """
    paths = _paths(pathogen)
    if not paths["variants"].is_file():
        return None

    total = 0
    by_type: dict[str, int] = {}
    genomes: set[str] = set()
    with paths["variants"].open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            by_type[row.get("variant_type", "?")] = by_type.get(row.get("variant_type", "?"), 0) + 1
            genomes.add(row.get("accession", ""))

    in_g4 = 0
    per_locus: dict[str, int] = {}
    if paths["g4_variants"].is_file():
        with paths["g4_variants"].open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                in_g4 += 1
                per_locus[row.get("atlas_id", "?")] = per_locus.get(row.get("atlas_id", "?"), 0) + 1

    return {
        "n_variants": total,
        "n_genomes": len(genomes),
        "by_type": by_type,
        "n_in_atlas_loci": in_g4,
        "fraction_in_atlas_loci": (in_g4 / total) if total else None,
        "per_locus": dict(sorted(per_locus.items(), key=lambda kv: -kv[1])[:20]),
        "mean_per_genome": (total / len(genomes)) if genomes else None,
    }


def molecular_clock(pathogen: str) -> dict | None:
    """TreeTime's root-to-tip regression, or None if dating has not run.

    Surfaced because a weak clock quietly undermines anything that reads
    calendar time, and the number lives in a file nobody opens. On the
    FMDV 2026 corpus the fit is r^2 = 0.06 with 55 outliers: the rate
    estimate is not usable for dating, and any statement of the form "this
    lineage expanded over N years" inherits that uncertainty. Reporting it
    beside the tree is the difference between a caveat and a surprise.
    """
    base = _paths(pathogen)["mldist"].parent
    path = base / "treetime_output" / "molecular_clock.txt"
    if not path.is_file():
        return None
    text = path.read_text()
    found: dict = {"path": str(path)}
    for pattern, key in ((r"--rate:\s*([0-9.eE+-]+)", "rate"),
                         (r"--r\^2:\s*([0-9.eE+-]+)", "r_squared")):
        match = re.search(pattern, text)
        if match:
            try:
                found[key] = float(match.group(1))
            except ValueError:
                pass
    outliers = base / "treetime_output" / "outliers.tsv"
    if outliers.is_file():
        found["n_outliers"] = max(0, sum(1 for _ in outliers.open()) - 1)

    r2 = found.get("r_squared")
    found["usable_for_dating"] = bool(r2 is not None and r2 >= 0.5)
    found["interpretation"] = (
        f"Weak temporal signal (r^2 = {r2}). Root-to-tip divergence barely tracks "
        "sampling date, so the time-scaling is unreliable and any date-dependent "
        "result inherits that uncertainty. Topology and rooting are unaffected."
        if r2 is not None and r2 < 0.5
        else f"Clock fit r^2 = {r2}."
    )
    return found
