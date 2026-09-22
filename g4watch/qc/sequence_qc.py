"""Per-sequence QC gate (concept paper Part J.2 / architecture Section 5's
H.2.2 table): genome completeness, N-content, and date precision.

Every threshold below is the concept paper's own stated MINIMUM (not the
preferred value) — records failing the minimum are excluded; records between
minimum and preferred are included but the distinction is preserved in
`QcResult` so a completeness/quality report can still show it.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

DEFAULT_MIN_COMPLETENESS_FRACTION = 0.90
PREFERRED_MIN_COMPLETENESS_FRACTION = 0.95
DEFAULT_MAX_N_CONTENT_FRACTION = 0.05
PREFERRED_MAX_N_CONTENT_FRACTION = 0.001

_VALID_BASES = frozenset("ACGTU")
_YEAR_RE = re.compile(r"(?<!\d)(1[6-9]\d{2}|20\d{2})(?!\d)")  # 1600-2099, defensive bound


@dataclass(frozen=True)
class QcResult:
    passed: bool
    completeness_fraction: float
    n_content_fraction: float
    has_year_precision_date: bool
    has_month_precision_date: bool
    reasons_failed: tuple[str, ...] = field(default_factory=tuple)


def genome_completeness_fraction(seq_length: int, reference_length: int) -> float:
    if reference_length <= 0:
        raise ValueError("reference_length must be positive")
    return seq_length / reference_length


def n_content_fraction(sequence: str) -> float:
    """Fraction of bases that are NOT one of the unambiguous A/C/G/T/U
    (case-insensitive). Catches N and all other IUPAC ambiguity codes,
    matching the concept paper's own "N-content: ambiguous bases" framing
    (not literally only the character 'N')."""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    ambiguous = sum(1 for base in seq if base not in _VALID_BASES)
    return ambiguous / len(seq)


def extract_year(date_str: str | None) -> int | None:
    if not date_str:
        return None
    match = _YEAR_RE.search(date_str)
    return int(match.group(1)) if match else None


_MONTH_NAMES = (
    "jan", "feb", "mar", "apr", "may", "jun",
    "jul", "aug", "sep", "oct", "nov", "dec",
)


def has_month_precision(date_str: str | None) -> bool:
    """True if the date string has both a year AND a month component.
    Recognizes GenBank's collection_date conventions: 'DD-Mon-YYYY',
    'Mon-YYYY', 'YYYY-MM', or a bare 'YYYY-MM-DD'."""
    if not date_str or extract_year(date_str) is None:
        return False
    lowered = date_str.lower()
    if any(month in lowered for month in _MONTH_NAMES):
        return True
    return bool(re.search(r"\b\d{4}-\d{2}(-\d{2})?\b", date_str))


def extract_month(date_str: str | None) -> int | None:
    """The calendar month (1-12), when the date has month precision.

    Recognizes exactly the formats ``has_month_precision`` does, and
    returns ``None`` on anything it does not -- a caller should always
    check the two together, or just treat ``None`` as "not precise
    enough", which is what every existing caller of
    ``has_month_precision`` already does.

    A GenBank collection_date is sometimes a range
    ('18-Jul-2024/21-Jul-2024', a collector unsure of the exact day). Only
    the first date in the range is read, on the assumption a range short
    enough to be given as a range does not cross a month boundary -- true
    of every real example seen in this corpus, and if it is ever wrong
    for some record, the record still gets a real month, just from the
    wrong half of a range that already could not be narrowed to a single
    day.
    """
    if not has_month_precision(date_str):
        return None
    first = date_str.split("/")[0].strip()
    lowered = first.lower()
    for index, name in enumerate(_MONTH_NAMES):
        if name in lowered:
            return index + 1
    numeric = re.search(r"\b\d{4}-(\d{2})(?:-\d{2})?\b", first)
    if numeric:
        month = int(numeric.group(1))
        if 1 <= month <= 12:
            return month
    return None


def evaluate_sequence_qc(
    seq_length: int,
    reference_length: int,
    sequence_for_n_content: str,
    collection_date: str | None,
    min_completeness: float = DEFAULT_MIN_COMPLETENESS_FRACTION,
    max_n_content: float = DEFAULT_MAX_N_CONTENT_FRACTION,
) -> QcResult:
    completeness = genome_completeness_fraction(seq_length, reference_length)
    n_content = n_content_fraction(sequence_for_n_content)
    year_ok = extract_year(collection_date) is not None
    month_ok = has_month_precision(collection_date)

    reasons: list[str] = []
    if completeness < min_completeness:
        reasons.append(
            f"genome completeness {completeness:.1%} below minimum {min_completeness:.0%}"
        )
    if n_content > max_n_content:
        reasons.append(f"N-content {n_content:.1%} above maximum {max_n_content:.0%}")
    if not year_ok:
        reasons.append("collection date missing year-level precision")

    return QcResult(
        passed=(len(reasons) == 0),
        completeness_fraction=completeness,
        n_content_fraction=n_content,
        has_year_precision_date=year_ok,
        has_month_precision_date=month_ok,
        reasons_failed=tuple(reasons),
    )
