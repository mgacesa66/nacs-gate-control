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
