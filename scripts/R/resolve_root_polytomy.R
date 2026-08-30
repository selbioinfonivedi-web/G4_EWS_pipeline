#!/usr/bin/env Rscript
# Resolves a trifurcating (or otherwise polytomous) ROOT into a bifurcating
# one via ape::multi2di() (deterministic, adds only zero-length branches --
# does not change any branch's biological length, does not move any tip).
#
# When this is safe to use: on a tree from a tool that already computed and
# applied a real root (e.g. TreeTime's --reroot), but writes it in the
# standard Newick "unrooted" trifurcating-root display convention. Applying
# multi2di() there just reformats an ALREADY-MADE rooting decision into the
# shape ape::ace() requires.
#
# When this is NOT safe to use: on a genuinely unrooted tree (e.g. raw
# IQ-TREE ML output) with no real root decision behind it -- multi2di()
# would then silently fabricate a root at an arbitrary point (wherever the
# first taxon happens to be listed in the Newick), which is exactly the
# kind of hidden scientific decision this project's house rule forbids.
# Callers must know which case they are in; this script does not guess.
#
# Usage: Rscript resolve_root_polytomy.R <input_tree> <output_newick>
#   <input_tree> may be Newick or Nexus (auto-detected by ape::read.tree /
#   read.nexus based on file content).

suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: resolve_root_polytomy.R <input_tree> <output_newick>")
}
input_path <- args[1]
output_path <- args[2]

first_line <- readLines(input_path, n = 1)
if (grepl("^#NEXUS", first_line, ignore.case = TRUE)) {
  tree <- read.nexus(input_path)
  if (is.list(tree) && !inherits(tree, "phylo")) {
    tree <- tree[[1]]  # multiple trees in file -- take the first
  }
} else {
  tree <- read.tree(input_path)
}

tree <- multi2di(tree)

if (!is.rooted(tree) || !is.binary(tree)) {
  stop("Tree is still not rooted/binary after multi2di() -- unexpected; inspect the input tree manually.")
}

write.tree(tree, output_path)
cat("OK\n")
