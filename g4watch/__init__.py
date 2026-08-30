"""G4-WATCH — G-quadruplex genomic early-warning framework for livestock viruses.

ICAR-NIVEDI. See ``docs/`` for methods, and ``G4_WATCH_Build_Architecture.md``
for the design this package implements.

Two conditions must both hold before any deployment claim is made for a
pathogen (Build Architecture Section 12):

1. ``operational_mode: true`` in ``config/<pathogen>.yaml``; and
2. a ``SUPPORTED`` D.H1 verdict recorded in the testing ledger.

:mod:`g4watch.gating` enforces both. Stage 5 (scoring) and Stage 6
(scored reporting) refuse to run on real data until they hold.
"""

from __future__ import annotations

__version__ = "1.0.0"

__all__ = ["__version__"]
