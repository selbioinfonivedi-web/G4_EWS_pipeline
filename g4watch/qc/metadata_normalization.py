"""Metadata field normalization. Split out from sequence_qc.py because this
is cleanup of free-text submitter conventions, not a pass/fail QC check.

Motivating real finding (Sprint 3, real FMDV corpus): serotype was recorded
as "Asia 1", "Asia1", "Asia-1", and "FMDV-Asia 1" for the SAME serotype
across different submissions, and "SAT2" vs "SAT 2" similarly — normalizing
whitespace/hyphenation/prefix variants is necessary before any honest
per-serotype breakdown, or a report would undercount every serotype that
happens to have inconsistent submitter formatting.
"""

from __future__ import annotations

import re

_FMDV_PREFIX_RE = re.compile(r"^FMDV[\s-]*", re.IGNORECASE)
_WHITESPACE_HYPHEN_RE = re.compile(r"[\s-]+")


def normalize_serotype(raw: str) -> str:
    """Uppercases and strips whitespace/hyphens/an optional leading 'FMDV'
    prefix, so 'Asia 1', 'Asia1', 'Asia-1', and 'FMDV-Asia 1' all normalize
    to 'ASIA1'. Deliberately does NOT attempt to resolve lineage-vs-serotype
    confusion in the raw data (e.g. a record recorded as "Pan Asia O" mixes
    a lineage name into the serotype field) — that would require domain
    knowledge this function cannot verify, so such values are normalized
    only at the whitespace/case level and reported as their own distinct
    (if messy) category rather than silently reassigned."""
    if not raw:
        return ""
    without_prefix = _FMDV_PREFIX_RE.sub("", raw.strip())
    return _WHITESPACE_HYPHEN_RE.sub("", without_prefix).upper()
