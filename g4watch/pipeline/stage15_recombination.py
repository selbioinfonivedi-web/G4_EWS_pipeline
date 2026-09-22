"""Stage 1.5 — recombination screening (Build Architecture Section 11).

Mandatory for every pathogen and run *before* phylogenetics, because
ancestral-state reconstruction on a recombinant alignment reconstructs a
history that never happened. The per-pathogen tier comes from config:
``high_priority`` (LSDV) routes a significant PHI result to a non-tree-
based estimate rather than merely flagging it.

This module only wires config to :func:`screen_recombination`; the PHI
test itself, its parsing and its tiering live in
``g4watch/phylo/recombination_screen.py`` and are unit-tested there.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..config import ConfigError, PathogenConfig
from ..phylo.recombination_screen import (
    PhiExecutionError,
    RecombinationScreenResult,
    default_phi_binary,
    screen_recombination,
)


@dataclass(frozen=True)
class Stage15Result:
    screen: RecombinationScreenResult
    alignment_path: Path
    completed: bool

    def summary(self) -> str:
        phi = self.screen.phi
        p = "undefined" if phi.phi_p_value is None else f"{phi.phi_p_value:.4g}"
        lines = [
            f"PHI test on {self.alignment_path.name}: {phi.n_sequences} sequences, "
            f"{phi.n_informative_sites} informative sites, p = {p}",
            f"  tier: {self.screen.tier.value}",
            f"  significant: {self.screen.significant}",
        ]
        if phi.undefined_reason:
            lines.append(f"  undefined because: {phi.undefined_reason}")
        if self.screen.route_to_non_tree_based_estimate:
            lines.append(
                "  ROUTING: significant recombination in a high-priority pathogen — downstream "
                "analysis must use a non-tree-based estimate, not the ML tree."
            )
        elif self.screen.significant:
            lines.append(
                "  FLAGGED: significant recombination detected at standard tier — recorded and "
                "carried forward; interpret tree-based results with this in mind."
            )
        return "\n".join(lines)


def run_stage15_recombination(
    config: PathogenConfig,
    *,
    alignment_path: Path,
    phi_binary: Path | None = None,
) -> Stage15Result:
    if not alignment_path.exists():
        raise ConfigError(f"Alignment not found: {alignment_path}")

    # A second, independently-frozen copy of this constant used to live
    # here, computed from THIS file's own __file__ — same bug as the one
    # documented on default_phi_binary(), duplicated rather than shared,
    # so fixing one copy silently left the other broken. This is the path
    # `g4watch recombination` (the CLI, and the only thing Nextflow's
    # RECOMBINATION_SCREEN actually calls) runs through; the other copy
    # was fixed but never exercised by a real containerised run.
    binary = phi_binary or default_phi_binary()
    if not Path(binary).exists():
        # No silent skip: an unrun mandatory screen must stop the pipeline,
        # because minimum_data_gate is about to be told it completed.
        raise PhiExecutionError(
            f"PhiPack binary not found at {binary}. Stage 1.5 is mandatory (Section 11) and cannot "
            "be skipped — build it with `make vendor` (see docs/installation.md) or pass --phi-binary."
        )

    recombination = config.raw.get("recombination") or {}
    screen = screen_recombination(
        alignment_path,
        tier=config.recombination_tier,
        alpha=float(recombination.get("alpha", 0.05)),
        phi_binary=binary,
    )
    return Stage15Result(screen=screen, alignment_path=alignment_path, completed=True)
