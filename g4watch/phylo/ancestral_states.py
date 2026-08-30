"""ML ancestral state reconstruction for a discrete character (e.g. a G4
Atlas locus's Conserved/Disrupted/Gained/... state per tip), via
scripts/R/ancestral_state_reconstruction.R (ape::ace(), ER model).

Substitution note: see that R script's own header comment for the full
reasoning — this uses ape/phangorn's ML reconstruction (an explicitly
approved alternative in the architecture) rather than phytools' stochastic
character mapping, because phytools is not installed in this environment.

This module is the "required upstream dependency" architecture Section 6
calls for: g4c_phylo/g4d_phylo (Sprint 6 scope, not yet built) will consume
an AncestralReconstructionRun's per-node most-likely states via
phylo/clade_collapse.py, rather than working from raw tip proportions.
"""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

DEFAULT_R_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "R" / "ancestral_state_reconstruction.R"
DEFAULT_RESOLVE_ROOT_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "R" / "resolve_root_polytomy.R"


@dataclass(frozen=True)
class AncestralNodeResult:
    node_id: int
    state_probabilities: dict[str, float]
    most_likely_state: str
    # The set of tip labels descending from this node -- a content-based,
    # tool-independent identifier. Callers correlating this result against
    # a tree built independently in another library (e.g. Bio.Phylo) MUST
    # match on tip_set, never on node_id -- there is no guarantee ape's
    # internal node numbering agrees with another tool's own convention
    # (see ancestral_state_reconstruction.R's tip_set comment).
    tip_set: frozenset[str]


@dataclass(frozen=True)
class AncestralReconstructionRun:
    nodes: tuple[AncestralNodeResult, ...]
    log_likelihood: float

    def node(self, node_id: int) -> AncestralNodeResult:
        for n in self.nodes:
            if n.node_id == node_id:
                return n
        raise KeyError(f"No reconstructed node with id {node_id}")

    def node_by_tip_set(self, tip_set: frozenset[str]) -> AncestralNodeResult:
        """The robust lookup — see AncestralNodeResult.tip_set's docstring
        on why this, not `node()`, is the correct way to correlate results
        against a tree built independently in another tool."""
        for n in self.nodes:
            if n.tip_set == tip_set:
                return n
        raise KeyError(f"No reconstructed node with tip set {tip_set}")


class AncestralReconstructionError(RuntimeError):
    """Raised on any Rscript failure — never silently swallowed."""


def _parse_r_float(value: str) -> float:
    """R's write.table() writes missing values as the literal string 'NA'
    (not 'NaN', which Python's float() would accept directly) -- found
    running this against real data where ace() failed to converge at some
    nodes (see ancestral_state_reconstruction.R's NaN/NA handling)."""
    if value == "NA":
        return float("nan")
    return float(value)


def _parse_output_tsv(text: str) -> AncestralReconstructionRun:
    lines = text.strip().splitlines()
    header = lines[0].split("\t")
    # Column order (fixed, written by ancestral_state_reconstruction.R):
    # node_id, <state columns...>, most_likely_state, loglik, tip_set
    state_columns = header[1:-3]
    nodes = []
    loglik = 0.0
    for line in lines[1:]:
        fields = line.split("\t")
        node_id = int(fields[0])
        n_states = len(state_columns)
        probs = {state: _parse_r_float(value) for state, value in zip(state_columns, fields[1 : 1 + n_states])}
        most_likely_state = fields[1 + n_states]
        loglik = _parse_r_float(fields[2 + n_states])
        tip_set = frozenset(fields[3 + n_states].split(","))
        nodes.append(
            AncestralNodeResult(
                node_id=node_id,
                state_probabilities=probs,
                most_likely_state=most_likely_state,
                tip_set=tip_set,
            )
        )
    return AncestralReconstructionRun(nodes=tuple(nodes), log_likelihood=loglik)


def resolve_treetime_root_polytomy(
    treetime_tree_path: str | Path,
    output_newick_path: str | Path,
    rscript_path: str | Path = DEFAULT_RESOLVE_ROOT_SCRIPT,
) -> None:
    """Reformats a TreeTime-rooted tree's trifurcating root into the
    bifurcating form ape::ace() requires (see resolve_root_polytomy.R's
    header for exactly why this is safe for TreeTime output specifically
    and NOT a general-purpose "make any tree rooted" utility). Callers must
    only use this on a tree they know already reflects a real rooting
    decision (TreeTime's, an outgroup's, etc.) — never on a raw, genuinely
    unrooted ML tree, which this would silently root arbitrarily."""
    result = subprocess.run(
        ["Rscript", str(rscript_path), str(treetime_tree_path), str(output_newick_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise AncestralReconstructionError(
            f"resolve_root_polytomy.R failed:\n{result.stdout}\n{result.stderr}"
        )


def reconstruct_ancestral_states(
    tree_newick_path: str | Path,
    tip_states: dict[str, str],
    rscript_path: str | Path = DEFAULT_R_SCRIPT,
) -> AncestralReconstructionRun:
    if len(set(tip_states.values())) < 2:
        raise AncestralReconstructionError(
            "Need at least 2 distinct states across tips for ancestral reconstruction to be meaningful "
            f"(got: {set(tip_states.values())})"
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        states_path = Path(tmpdir) / "tip_states.tsv"
        output_path = Path(tmpdir) / "output.tsv"
        with states_path.open("w") as fh:
            for tip, state in tip_states.items():
                fh.write(f"{tip}\t{state}\n")

        result = subprocess.run(
            ["Rscript", str(rscript_path), str(tree_newick_path), str(states_path), str(output_path)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise AncestralReconstructionError(
                f"ancestral_state_reconstruction.R failed:\n{result.stdout}\n{result.stderr}"
            )
        return _parse_output_tsv(output_path.read_text())
