"""Config-driven pipeline stages.

Each module here is one stage of the Build Architecture Section 6 order,
exposed as a plain function that takes a validated
:class:`g4watch.config.PathogenConfig` plus explicit paths and returns a
result dataclass. Printing, file writing and exit codes belong to
:mod:`g4watch.cli`; the Nextflow modules in ``workflow/`` call that CLI.

Keeping the stages importable and side-effect-light is what makes the
pipeline testable without Nextflow, and what lets ``tests/integration/``
exercise the real stage order on small fixtures.

Stage order (Section 6)::

    Stage 0    atlas construction            stage0_atlas
    Stage 1    acquisition + QC + alignment  stage1_qc
    Stage 1.5  recombination screening       stage15_recombination
    Stage 2    phylogenomics + ancestral     (external: IQ-TREE2/TreeTime)
    Stage 3    variant analysis              (g4watch.variants)
    Stage 4    G4 surveillance metrics       stage45_dh1 (per-locus clades)
    Stage 4.5  GC-confound gate + D.H1       stage45_dh1
    --------- GATE: g4watch.gating ----------
    Stage 5    scoring                       stage5_scoring  (fail-closed)
    Stage 6    reporting                     stage6_reporting (fail-closed)
"""
