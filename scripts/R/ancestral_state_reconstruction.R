#!/usr/bin/env Rscript
# ML ancestral state reconstruction for a discrete character, via ape::ace().
#
# Substitution note (Sprint 4): the architecture (Concept Paper v2 Section
# G.4) names phytools stochastic character mapping as the primary method.
# phytools is NOT installed in this environment; ape and phangorn ARE, and
# the architecture's own text explicitly lists "ML ancestral reconstruction
# (ape/phangorn in R)" as an accepted alternative method, not something
# invented for this substitution. This script implements that alternative:
# a single best-fit ML reconstruction (equal-rates Mk model, Pagel 1994) via
# ape::ace(), not phytools' stochastic mapping (which samples many possible
# character histories, not just the marginal ancestral state probabilities
# at each node). If phytools is installed later, a stochastic-mapping
# script should be added alongside this one, not replace it -- the two
# answer related but different questions.
#
# Usage: Rscript ancestral_state_reconstruction.R <tree.nwk> <tip_states.tsv> <output.tsv>
#   tip_states.tsv: two columns, tip_label \t state, NO header row.
#
# PRECONDITION: the input tree must already be rooted and fully bifurcating
# (ace() requires this and fails with a cryptic error otherwise -- found the
# hard way in Sprint 4 running this against IQ-TREE's raw ML output, which
# is unrooted by convention since GTR is a reversible model). This script
# does NOT auto-root the tree: which rooting to use (midpoint, outgroup,
# temporal/TreeTime) is a real modeling decision, and choosing one silently
# on the caller's behalf would hide that choice, the opposite of this
# project's no-silent-fallback rule. Root the tree upstream (this pipeline's
# intended source is TreeTime's time-rooted output, not an arbitrary
# midpoint rooting of the raw ML tree) and pass the rooted tree in.

suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: ancestral_state_reconstruction.R <tree.nwk> <tip_states.tsv> <output.tsv>")
}
tree_path <- args[1]
states_path <- args[2]
output_path <- args[3]

tree <- read.tree(tree_path)

if (!is.rooted(tree)) {
  stop(paste(
    "Tree is not rooted. ace() requires a rooted, bifurcating tree.",
    "Root it upstream (e.g. via TreeTime's time-calibrated output, or",
    "ape::root()/midpoint rooting if no temporal rooting is available)",
    "and pass the rooted tree to this script -- it will not choose a",
    "rooting method for you."
  ))
}
if (!is.binary(tree)) {
  stop(paste(
    "Tree is not fully bifurcating (contains a polytomy). ace() requires a",
    "fully resolved tree -- resolve polytomies upstream (e.g.",
    "ape::multi2di(), which arbitrarily but deterministically resolves",
    "them with zero-length branches) before calling this script."
  ))
}
states_df <- read.table(states_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                         col.names = c("tip", "state"))

tip_state <- setNames(states_df$state, states_df$tip)
missing_tips <- setdiff(tree$tip.label, names(tip_state))
if (length(missing_tips) > 0) {
  stop(paste("Tree tips with no assigned state:", paste(missing_tips, collapse = ", ")))
}
tip_state <- tip_state[tree$tip.label]  # order to match tree$tip.label exactly

state_factor <- factor(tip_state)
levels_used <- levels(state_factor)

if (length(levels_used) < 2) {
  stop("Need at least 2 distinct states across tips for ancestral reconstruction to be meaningful")
}

fit <- ace(state_factor, tree, type = "discrete", method = "ML", model = "ER")

anc_probs <- fit$lik.anc
colnames(anc_probs) <- levels_used

n_nan_rows <- sum(apply(anc_probs, 1, function(row) any(is.nan(row)) || any(is.na(row))))
if (n_nan_rows > 0) {
  warning(paste(
    n_nan_rows, "of", nrow(anc_probs),
    "internal nodes have NaN/NA ancestral state probabilities",
    "(ace() likely did not fully converge -- with many states and many",
    "tips, the ER model's optimizer can fail at some nodes without",
    "erroring the whole fit). These are reported as 'UNRESOLVED' rather",
    "than silently picking an arbitrary state."
  ))
}

n_tips <- length(tree$tip.label)
internal_node_ids <- (n_tips + 1):(n_tips + tree$Nnode)

most_likely_state <- vapply(seq_len(nrow(anc_probs)), function(i) {
  row <- anc_probs[i, ]
  if (any(is.nan(row)) || any(is.na(row))) {
    return("UNRESOLVED")
  }
  levels_used[which.max(row)]
}, character(1))

out <- data.frame(node_id = internal_node_ids, anc_probs, check.names = FALSE)
out$most_likely_state <- most_likely_state
out$loglik <- fit$loglik

# tip_set: the sorted, comma-joined tip labels descending from each internal
# node. Sprint 6 rationale: a caller matching this output against a tree
# built independently in another tool (e.g. Bio.Phylo in Python) must NOT
# assume ape's internal node_id numbering matches that other tool's own
# node ordering -- there is no guarantee they agree, and a silent mismatch
# would corrupt every downstream clade-collapse result without any obvious
# symptom. tip_set is a content-based, tool-independent join key: match on
# WHICH tips a node subtends, never on a numeric ID assigned by a different
# library's internal convention.
out$tip_set <- vapply(internal_node_ids, function(node_id) {
  paste(sort(extract.clade(tree, node_id)$tip.label), collapse = ",")
}, character(1))

write.table(out, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("OK\n")
