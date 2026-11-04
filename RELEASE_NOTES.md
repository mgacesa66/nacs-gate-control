## v0.2.0.dev3 — arXiv Step 3 (research extension only)

- Four physically constrained objective functions: isolated-molecule access, selected-aux51 penalty, full logical-channel fidelity, and full channel with total-spectator penalty.
- Same fixed microwave orientation, field amplitude, two-segment transfer controls, 2-µs finite phase ramps, Gaussian echo, and 1.30-ms duration bound.
- Saved 8 global-search seed checkpoints, four 8-state local refinements, four identical-budget 11-state refinements, nominal controls, corrected-basis nominal checks, and 24 paired held-out draws for each of two hypothetical uncertainty envelopes.
- 25-state union-selected basis cross-check with 625-dimensional pair state vectors; column propagator verified against original dense exponential to <3e-11 in two cases.
- Core conclusion: isolated transfer is not predictive of complete-gate fidelity. Guarded full-channel improvement over other full-gate objectives is modest, and broad-noise robustness remains poor. Earlier strong aux51 cancellation depended on idealized sharp phase steps, not reproduced by finite-phase-ramp search.
- Source/release limits and actual differing Powell early-termination function-call counts disclosed. No figures, no changes to prior archives.

## v0.2.0.dev2 — arXiv Step 2 (research development only)

- Systematic microwave-coupling-ranked 11/13/15/17/19/21/23/25-state pair-basis gate propagation.
- Independent resonant dipole-dipole interaction screening and 25-state cross-screen selection.
- All Step 1 parameter probes and larger-basis ellipticity follow-ups; full complex logical map and phase-alignment checks.
- Central qualitative results survive; intermediate truncations are nonmonotonic, and the study does NOT certify complete 128² convergence.
- No figures. No public v0.1.0 or frozen Step 1 content was modified.

# Research extension v0.2.0.dev1 — 2026-11-02

This is a **separate research extension**, not a revision of the frozen GitHub v0.1.0 baseline.

- Frozen six-step arXiv preparation plan, with Step 1 explicitly scoped.
- Recomputed 128-level NaCs Hamiltonian and tracked eigenvectors for 11 hypothetical physical-parameter probes.
- Compared three microwave calibration regimes and four unchanged stored control protocols (132 gate outputs).
- Distinguished fixed nominal logical-Z calibration from per-static-scenario Z recalibration.
- Added leakage redistribution, spectral shifts, dressed-state label-overlap diagnostics, and CSV/JSON results.
- Six new regression and invariant tests, in addition to nine from the original repository.

**Important qualification:** frequency and target-Rabi tracking are idealized diagnostics, not power-dressed experimental spectroscopic calibrations. Physical-parameter variations are hypothetical; complete larger-basis convergence for the perturbed Hamiltonians is deferred to Step 2. No figures.

# Release notes — v0.1.0

**Type:** GitHub research-software release and technical report; no manuscript figures.

## Included

- Directly importable Python 3 package from the chronological research code, without nested ZIP extraction at runtime.
- Immutable, SHA-256 checked NaCs spectroscopy v1.0 dressed spectrum, raw complex microwave dipole matrices, and archived control vectors.
- Canonical **11-state** molecular-gate model, with corrected strongly coupled `aux51` and `aux79` spectators.
- Square echo, truncated-Gaussian echo, stage-resolved leakage, and spectator-guarded two-piece phase-transfer control.
- Historical spectroscopy, 49-point stabilization, and leakage-redistribution numerical datasets as text/CSV/JSON.
- New no-figure release audit: 9-/11-/13-/ranked-15-state fixed-control checks, seven spectroscopy residuals, seven physical-parameter perturbations with fixed absolute microwave carriers, paired bootstrap of archived simulated errors.
- Tests and GitHub Actions workflow.

## What is NOT included

Experimental source data from Harvard Dataverse, previously generated plots, earlier draft PRA PDFs, or all chronological research ZIP archives. Their source hashes are preserved. The old nine-state code is not used as the release's canonical model.

## Limitations requiring work before an arXiv paper

- Effective optical and microwave polarizations are not experimentally confirmed.
- Spectroscopic agreement obtained by inferred microwave-power shifts is not independently verified.
- Selected basis enlargement is not complete 128² pair convergence.
- Fixed-carrier gate fidelity can be very sensitive to optical ellipticity; see `docs/MODEL_VALIDITY.md`.
- No experimental noise distribution, molecular motion, decoherence or irreversible loss.
- Local rather than certified globally optimal controls.

## Verification

The release tests were executed against the complete source tree with one-thread BLAS, including exact pair evolution and ranked basis enlargement. See the generated `benchmarks/data/release_audit/` outputs and run commands in `docs/REPRODUCIBILITY.md`.
