"""Pathogen configuration loading and validation.

Every parameter that can change a scientific result lives in
``config/<pathogen>.yaml`` rather than in code, so that a run is fully
described by the triple (git commit, config file, accession list). This
module is the single place those files are read, and it validates them
strictly: an unknown key, a missing required key, or an out-of-range
value raises :class:`ConfigError` rather than being silently defaulted.
That is the same no-silent-fallback rule the rest of the package follows
(Build Architecture Section 3) applied to configuration.

Two properties are deliberately NOT configurable, because making them
configurable would make them bypassable:

* the Appendix C minimum-data floor (``validation/minimum_data_gate.py``)
* the D.H1 gate's role as a hard precondition for Stage 5 — see
  :func:`g4watch.gating.assert_scoring_permitted`.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from .phylo.recombination_screen import RecombinationTier

REPO_ROOT = Path(__file__).resolve().parents[1]


def config_dir() -> Path:
    """Where named pathogens (``"fmdv"``) are looked up.

    Resolved at call time, in this order:

    1. ``$G4WATCH_CONFIG_DIR``, when set — the explicit answer, and the one
       a deployment should use.
    2. ``<package parent>/config`` — correct in a repository checkout, where
       the package sits beside the config directory.
    3. ``./config`` relative to the working directory — correct for an
       INSTALLED package, where the package lives in site-packages and its
       parent holds no config at all.

    Rule 3 exists because of a real failure: inside the core container the
    package is installed, so rule 2 resolved to
    ``/usr/local/lib/python3.12/site-packages/config`` and every named
    pathogen failed with "Available pathogens: none". The image could not
    run a single stage for any pathogen, even with the repository mounted
    and the working directory set to it.
    """
    override = os.environ.get("G4WATCH_CONFIG_DIR")
    if override:
        return Path(override)
    packaged = REPO_ROOT / "config"
    if packaged.is_dir():
        return packaged
    return Path.cwd() / "config"


#: Kept for callers that import it directly. Prefer :func:`config_dir`,
#: which re-resolves rather than freezing the value at import time.
CONFIG_DIR = REPO_ROOT / "config"

# Top-level sections every pathogen config must define. Presence is
# required even when the pathogen is unprovisioned and the values are
# null, so that an incomplete config is visibly incomplete rather than
# quietly missing a whole stage's settings.
_REQUIRED_SECTIONS = frozenset(
    {
        "pathogen",
        "display_name",
        "genome_type",
        "reference",
        "corpus",
        "qc",
        "alignment",
        "recombination",
        "phylogenetics",
        "g4_prediction",
        "control_regions",
        "atlas",
        "dh1_gate",
        "operational_mode",
    }
)

_KNOWN_SECTIONS = _REQUIRED_SECTIONS | {"taxonomy_id", "provisioned", "provisioning_notes"}


class ConfigError(ValueError):
    """Raised for any malformed, incomplete or out-of-range config.

    Never caught-and-defaulted inside this package: a bad config must stop
    a run, not silently change what it computes.
    """


def _require(mapping: dict, key: str, where: str) -> Any:
    if key not in mapping:
        raise ConfigError(f"{where}: missing required key {key!r}")
    return mapping[key]


def _fraction(value: Any, where: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ConfigError(f"{where}: expected a number in [0, 1], got {value!r}")
    if not 0.0 <= float(value) <= 1.0:
        raise ConfigError(f"{where}: expected a fraction in [0, 1], got {value!r}")
    return float(value)


@dataclass(frozen=True)
class PathogenConfig:
    """A validated pathogen configuration.

    ``raw`` keeps the whole parsed document so that stages can read
    pathogen-specific sections without this class having to grow a field
    per parameter; the attributes below are the ones the gating and
    orchestration layers depend on, and those are validated eagerly.
    """

    path: Path
    raw: dict
    pathogen: str
    display_name: str
    provisioned: bool
    operational_mode: bool
    recombination_tier: RecombinationTier
    recombination_enabled: bool
    dh1_alpha: float
    repo_root: Path

    # --- resolved paths (None when the pathogen is unprovisioned) ---

    def _resolve(self, *keys: str) -> Path | None:
        node: Any = self.raw
        for key in keys:
            if not isinstance(node, dict) or node.get(key) is None:
                return None
            node = node[key]
        return self.repo_root / str(node)

    @property
    def reference_fasta(self) -> Path | None:
        return self._resolve("reference", "fasta")

    @property
    def reference_accession(self) -> str | None:
        return self.raw["reference"].get("accession")

    @property
    def corpus_metadata_tsv(self) -> Path | None:
        return self._resolve("corpus", "metadata_tsv")

    @property
    def corpus_sequences_fasta(self) -> Path | None:
        return self._resolve("corpus", "sequences_fasta")

    @property
    def atlas_path(self) -> Path | None:
        return self._resolve("atlas", "path")

    @property
    def ledger_path(self) -> Path:
        # The ledger is study-wide and always required — a run that cannot
        # record its own result must not start.
        raw = self.raw["dh1_gate"].get("ledger")
        if not raw:
            raise ConfigError(f"{self.path}: dh1_gate.ledger must be set (the testing ledger is mandatory)")
        return self.repo_root / str(raw)

    @property
    def genome_type(self) -> str:
        """Declared genome type, e.g. ``ssRNA_positive``.

        A real property because consumers were reaching for it with
        ``getattr(config, "genome_type", "")``. It is a required config
        section but not a dataclass field, so that guess silently
        returned the empty string everywhere and the workstation and the
        report card both displayed a blank genome type.
        """
        return str(self.raw.get("genome_type") or "")

    @property
    def genome_length(self) -> int | None:
        value = (self.raw.get("reference") or {}).get("genome_length")
        return int(value) if value is not None else None

    @property
    def atlas_version(self) -> str:
        """The Atlas version the config declares.

        Read from the config rather than parsed out of the Atlas
        filename: a file called ``atlas.tsv`` is perfectly valid and its
        name says nothing about the version it holds.
        """
        return str((self.raw.get("atlas") or {}).get("version") or "")

    @property
    def lineage_field(self) -> str:
        return self.raw["corpus"].get("lineage_field") or "lineage"

    @property
    def lineage_fallback_fields(self) -> tuple[str, ...]:
        """Secondary metadata columns that may also carry the lineage.

        Consulted in order when the primary ``lineage_field`` is blank or
        less specific. The defaults are generic GenBank columns, not
        anything pathogen-specific; a config may override the list.
        """
        raw = (self.raw.get("corpus") or {}).get("lineage_fallback_fields")
        if raw is None:
            return ("organism", "isolate", "strain")
        if not isinstance(raw, list) or not all(isinstance(x, str) for x in raw):
            raise ConfigError(
                f"{self.path}:corpus.lineage_fallback_fields must be a list of strings, got {raw!r}"
            )
        return tuple(raw)

    @property
    def exclude_lineages(self) -> tuple[str, ...]:
        """Lineages dropped before any statistic is computed.

        Configured rather than passed on the command line because an
        exclusion changes what every downstream number means, and a run
        must be fully described by (commit, config, accession list). A
        lineage that no longer circulates cannot be the subject of early
        warning, and carrying it only to fail the per-lineage floor blocks
        the analysis of everything that does. The CLI can still add to
        this for an exploratory run; it cannot silently replace it.
        """
        raw = self.raw["corpus"].get("exclude_lineages") or []
        if isinstance(raw, str):
            raise ConfigError(
                f"{self.path}:corpus.exclude_lineages must be a list, not a string "
                f"(got {raw!r}); a bare string would silently exclude one character at a time"
            )
        if not isinstance(raw, list) or not all(isinstance(x, str) for x in raw):
            raise ConfigError(f"{self.path}:corpus.exclude_lineages must be a list of strings, got {raw!r}")
        return tuple(x.strip() for x in raw if x.strip())

    def require_provisioned(self) -> None:
        """Fail closed on a pathogen whose reference/corpus are not filled in.

        Scaffold configs (LSDV, PPRV, NDV, CSFV) carry the decisions the
        architecture already fixes but no fabricated accessions, so any
        attempt to run one stops here with the curator's own checklist.
        """
        if self.provisioned:
            return
        notes = self.raw.get("provisioning_notes", "").strip()
        raise ConfigError(
            f"{self.pathogen} is not provisioned (config/{self.pathogen.lower()}.yaml sets "
            f"provisioned: false), so there is no reference genome or corpus to run against.\n{notes}"
        )


def _validate(raw: dict, path: Path, repo_root: Path) -> PathogenConfig:
    if not isinstance(raw, dict):
        raise ConfigError(f"{path}: expected a YAML mapping at the top level")

    missing = sorted(_REQUIRED_SECTIONS - raw.keys())
    if missing:
        raise ConfigError(f"{path}: missing required section(s): {', '.join(missing)}")
    unknown = sorted(raw.keys() - _KNOWN_SECTIONS)
    if unknown:
        raise ConfigError(
            f"{path}: unknown top-level key(s): {', '.join(unknown)}. "
            "Config keys are validated strictly so a typo cannot silently disable a stage."
        )

    recombination = _require(raw, "recombination", str(path))
    tier_raw = _require(recombination, "tier", f"{path}:recombination")
    try:
        tier = RecombinationTier(tier_raw)
    except ValueError as exc:
        valid = ", ".join(t.value for t in RecombinationTier)
        raise ConfigError(f"{path}:recombination.tier: {tier_raw!r} is not one of {valid}") from exc

    if not recombination.get("enabled", False):
        # Section 11 makes the screen mandatory for every pathogen; a
        # config cannot opt out of it.
        raise ConfigError(
            f"{path}:recombination.enabled must be true — recombination screening is mandatory "
            "for every pathogen (Build Architecture Section 11), not a per-pathogen option."
        )

    operational_mode = _require(raw, "operational_mode", str(path))
    if not isinstance(operational_mode, bool):
        raise ConfigError(f"{path}:operational_mode must be a boolean, got {operational_mode!r}")

    qc = _require(raw, "qc", str(path))
    _fraction(_require(qc, "min_completeness_fraction", f"{path}:qc"), f"{path}:qc.min_completeness_fraction")
    _fraction(_require(qc, "max_n_content_fraction", f"{path}:qc"), f"{path}:qc.max_n_content_fraction")

    dh1 = _require(raw, "dh1_gate", str(path))
    alpha = _fraction(_require(dh1, "alpha", f"{path}:dh1_gate"), f"{path}:dh1_gate.alpha")
    if alpha <= 0.0:
        raise ConfigError(f"{path}:dh1_gate.alpha must be > 0, got {alpha}")

    pathogen = str(_require(raw, "pathogen", str(path)))
    provisioned = bool(raw.get("provisioned", True))

    if provisioned:
        # A config claiming to be runnable must actually point at things.
        for section, key in (("reference", "fasta"), ("reference", "accession"), ("corpus", "metadata_tsv")):
            if not raw.get(section, {}).get(key):
                raise ConfigError(
                    f"{path}: {section}.{key} must be set when provisioned: true "
                    f"(set provisioned: false for a scaffold config)"
                )

    return PathogenConfig(
        path=path,
        raw=raw,
        pathogen=pathogen,
        display_name=str(_require(raw, "display_name", str(path))),
        provisioned=provisioned,
        operational_mode=operational_mode,
        recombination_tier=tier,
        recombination_enabled=True,
        dh1_alpha=alpha,
        repo_root=repo_root,
    )


def _infer_repo_root(config_path: Path) -> Path:
    """The directory a config's relative paths are resolved against.

    A config in ``<project>/config/<name>.yaml`` describes ``<project>``,
    so the root is two levels up; a config sitting anywhere else describes
    the directory it lives in. Inferring from the config file rather than
    from this package's own location is what lets a config outside the
    installed tree — a test fixture, a deployment holding its data
    elsewhere — resolve its paths correctly.
    """
    parent = config_path.parent
    return parent.parent if parent.name == "config" else parent


def load_config(pathogen_or_path: str | Path, *, repo_root: Path | None = None) -> PathogenConfig:
    """Load and validate a pathogen config.

    Accepts either a pathogen name (``"fmdv"``, resolved against
    ``config/``) or an explicit path to a YAML file. Relative paths inside
    the config resolve against ``repo_root``, which defaults to the root
    inferred from the config's own location (see :func:`_infer_repo_root`).
    """
    candidate = Path(pathogen_or_path)
    if candidate.suffix in {".yaml", ".yml"}:
        path = candidate if candidate.is_absolute() else (Path.cwd() / candidate)
    else:
        path = config_dir() / f"{str(pathogen_or_path).lower()}.yaml"
    repo_root = repo_root or _infer_repo_root(path.resolve())

    if not path.exists():
        available = ", ".join(sorted(p.stem for p in config_dir().glob("*.yaml"))) or "none"
        raise ConfigError(f"No config at {path}. Available pathogens: {available}")

    try:
        with open(path) as handle:
            raw = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        # Surface a malformed config as a ConfigError like every other
        # config problem, so callers that validate many files can report
        # one bad file and carry on instead of dying on a raw traceback.
        raise ConfigError(f"{path}: not valid YAML — {exc}") from exc
    return _validate(raw, path, repo_root)


def available_pathogens() -> list[str]:
    """Pathogen names with a config file present, provisioned or not."""
    return sorted(p.stem for p in CONFIG_DIR.glob("*.yaml"))
