"""Metadata field normalization — cleanup of free-text submitter
conventions, not a pass/fail QC check.

WHY THIS IS CONFIG-DRIVEN. The lineage vocabulary of a virus is a
scientific fact about that virus, so it lives in the pathogen YAML beside
every other such fact. This module previously hard-coded one pathogen's
seven serotypes, its species-name prefixes and its subtype grammar. That
made every other pathogen either unsupported or, worse, silently scored
against the wrong vocabulary the moment its config set
``lineage_field: serotype``.

Nothing here knows the name of any pathogen. Give it a vocabulary and it
normalises against that; give it none and it resolves nothing rather than
guessing.

THE PROBLEM IT SOLVES. On a real corpus one lineage is recorded four ways
by four submitters -- "Asia 1", "Asia1", "Asia-1", and the same value with
the species name glued to the front -- and a second lineage differs only
by an interior space. Normalising whitespace, hyphenation and species
prefix is necessary before any honest per-lineage breakdown, or a report
undercounts every lineage whose submitters formatted it inconsistently.
Which names are canonical, and which prefixes to strip, are declared per
pathogen in its config; this module supplies only the mechanism.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

_WHITESPACE_HYPHEN_RE = re.compile(r"[\s-]+")
_TYPE_PREFIX_RE = re.compile(r"^TYPE", re.IGNORECASE)


@dataclass(frozen=True)
class LineageVocabulary:
    """One pathogen's lineage naming rules, read from its config.

    ``canonical`` is the closed set a record may resolve to. ``aliases``
    maps a verified alternative spelling or lineage name onto its parent.
    ``name_prefixes`` are species-name forms that appear glued to the
    front of a lineage value. ``subtype_pattern`` is an optional regex
    whose first group is the parent lineage, for historical subtype
    labels such as O1 or A22.
    """

    canonical: tuple[str, ...] = ()
    aliases: dict[str, str] = field(default_factory=dict)
    name_prefixes: tuple[str, ...] = ()
    subtype_pattern: str | None = None

    @classmethod
    def from_config(cls, config) -> LineageVocabulary:
        """Build from ``corpus.lineage_vocabulary`` in a pathogen config.

        An absent section yields an empty vocabulary, which resolves
        nothing. That is the correct behaviour for a pathogen whose
        lineage names are already clean, and for one nobody has curated
        yet — the alternative, guessing, is what this module exists to
        avoid.
        """
        raw = (config.raw.get("corpus") or {}).get("lineage_vocabulary") or {}
        return cls(
            canonical=tuple(str(x).upper() for x in raw.get("canonical") or ()),
            aliases={str(k).upper(): str(v).upper() for k, v in (raw.get("aliases") or {}).items()},
            name_prefixes=tuple(str(x).upper() for x in raw.get("name_prefixes") or ()),
            subtype_pattern=raw.get("subtype_pattern"),
        )

    # ── internals ──
    def _prefix_re(self) -> re.Pattern[str] | None:
        if not self.name_prefixes:
            return None
        joined = "|".join(re.escape(p) for p in self.name_prefixes)
        return re.compile(rf"^(?:{joined})", re.IGNORECASE)

    def _subtype_re(self) -> re.Pattern[str] | None:
        return re.compile(self.subtype_pattern) if self.subtype_pattern else None

    # ── public ──
    def normalize(self, raw: str) -> str:
        """Uppercase and strip whitespace, hyphens and a species prefix.

        Deliberately does NOT resolve lineage-versus-parent confusion —
        a value that mixes a lineage name into the lineage field is
        normalised only at the whitespace and case level and reported as
        its own distinct, messy category rather than silently reassigned.
        """
        if not raw:
            return ""
        text = raw.strip()
        prefix = self._prefix_re()
        if prefix is not None:
            text = prefix.sub("", text).lstrip(" -")
        return _WHITESPACE_HYPHEN_RE.sub("", text).upper()

    def _token(self, raw: str) -> str:
        """Reduce a free-text field to the token that may name a lineage.

        Strips the species prefix and an optional "type" infix, removes
        whitespace and hyphens, uppercases, then takes the text before
        the first "/" — GenBank isolate names follow
        LINEAGE/COUNTRY/ID/YEAR, so the lineage is in the first field.
        """
        token = _WHITESPACE_HYPHEN_RE.sub("", (raw or "").strip()).upper()
        prefix = self._prefix_re()
        if prefix is not None:
            token = prefix.sub("", token)
        token = _TYPE_PREFIX_RE.sub("", token)
        return token.split("/")[0]

    def match(self, token: str) -> str:
        """Resolve one token, or "" if it does not name a lineage.

        Deliberately strict: only exact canonical names, declared aliases
        and the configured subtype grammar resolve. A prefix match would
        map a country code like "AUS" onto a lineage "A", which is the
        kind of silent misassignment that corrupts the per-lineage floor.
        """
        if not token:
            return ""
        if token in self.canonical:
            return token
        if token in self.aliases:
            return self.aliases[token]
        subtype = self._subtype_re()
        if subtype is not None:
            found = subtype.match(token)
            if found and found.group(1) in self.canonical:
                return found.group(1)
        return ""

    def resolve(self, *fields: str) -> str:
        """Resolve a record to a canonical lineage, or "".

        Tries each field in order and returns the first that resolves, so
        a caller passes them most-authoritative-first::

            vocabulary.resolve(row["serotype"], row["organism"], row["isolate"])

        **Why a secondary field is trusted as a fallback.** On a real
        corpus, records frequently carry the lineage in the organism or
        isolate name while the dedicated field is blank or less specific.
        Consulting them recovered 266 of 269 sequences that appeared to
        have no lineage at all — including 93 that alone moved one lineage
        from below the Appendix C per-lineage floor to well above it.
        Those sequences were never missing a lineage; the parser was
        reading one of the several fields that record it.
        """
        for value in fields:
            resolved = self.match(self._token(value))
            if resolved:
                return resolved
        return ""


#: Resolves nothing. Used where no pathogen context is available, so the
#: absence of a vocabulary is visible as unresolved lineages rather than
#: as one pathogen's names quietly applied to another's data.
EMPTY_VOCABULARY = LineageVocabulary()


def normalize_lineage(raw: str, vocabulary: LineageVocabulary = EMPTY_VOCABULARY) -> str:
    """Module-level convenience wrapper around :meth:`LineageVocabulary.normalize`."""
    return vocabulary.normalize(raw)


def canonical_lineage(*fields: str, vocabulary: LineageVocabulary = EMPTY_VOCABULARY) -> str:
    """Module-level convenience wrapper around :meth:`LineageVocabulary.resolve`."""
    return vocabulary.resolve(*fields)
