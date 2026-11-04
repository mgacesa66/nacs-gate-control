# Publication plan status — 2026-11-02

The authoritative six-step plan remains **frozen without modification** in [`PUBLICATION_PLAN_FROZEN.md`](PUBLICATION_PLAN_FROZEN.md). This file is a versioned progress record, not a replacement plan.

- **Step 1 — completed**: conditional molecular-parameter/recalibration study (`v0.2.0.dev1`).
- **Step 2 — completed**: systematic screened-basis convergence audit up to 25 states per molecule; not a proof of full 128² pair-space convergence (`v0.2.0.dev2`).
- **Step 3 — completed as an exploratory matched-constraints, bounded-budget study**: four control objectives under identical fields, pulse bounds and finite phase smoothing, two global seeds each, equal-count 11-state refinement, independent 25-state check and paired held-out error probes (`v0.2.0.dev3`). See [`reports/PAPER_STEP3_MATCHED_OBJECTIVES.md`](reports/PAPER_STEP3_MATCHED_OBJECTIVES.md). The local eight-state Powell searches have equal budget ceilings but unequal actual evaluation counts; **no certified equal-compute global-optimum comparison or experimental control bandwidth validation** is claimed.
- **Step 4 — next**: literature/novelty audit.
- Step 5 — manuscript and figures, pending.
- Step 6 — freeze results and submit, pending.

Public GitHub v0.1.0, frozen spectral exports, Step 1 archive, and parent Step 2 archive are unchanged.
