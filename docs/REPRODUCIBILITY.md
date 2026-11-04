# Reproducing the release

## Requirements and environment

The code targets Python 3.10+ with NumPy, SciPy and SymPy. `pytest` is included in the development extra. The release has been validated using single-threaded OpenBLAS; other BLAS implementations may change floating-point roundoff and eigenvector phases in near-degenerate subspaces.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
OPENBLAS_NUM_THREADS=1 python -m pytest -q
OPENBLAS_NUM_THREADS=1 python scripts/quickstart.py
OPENBLAS_NUM_THREADS=1 python scripts/release_audit.py --mode all --output outputs/release_audit
python scripts/uncertainty_summary.py --output outputs/uncertainty_summary.json
```

No web access, experimental datasets or proprietary code is required to run the stored simulations. The paper's figures are intentionally excluded.

## Which source computes what?

| Command | Actual computation | Typical output |
|---|---|---|
| `pytest -q` | Immutable export checksum, Hermiticity, unitarity, actual Gaussian gate and aux51 propagation, selected-state convergence | Test results |
| `scripts/quickstart.py` | Three complete 11-state-per-molecule gate propagations | Terminal fidelity/leakage |
| `release_audit.py --mode spectroscopy` | Transparent comparison against seven archived measured lines | `spectroscopy_accounting.csv` |
| `release_audit.py --mode basis` | Fixed-controls 9/11/13/ranked-15 propagation for four stored controls | `basis_check.csv` |
| `release_audit.py --mode sensitivity` | Seven fresh 128-state microscopic re-diagonalizations and fixed-carrier full-gate evaluations | `parameter_sensitivity.csv` |
| `uncertainty_summary.py` | Paired bootstrap of stored 30-shot held-out errors, seeded reproducibly | JSON |

## Unit and phase conventions

- Hamiltonians \(H/h\) are in **cyclic kHz**, not angular frequency.
- Times are in **milliseconds**; propagation is \(U(t)=\exp[-2\pi i (H/h)t]\).
- Dipole transition matrices use **Debye**.
- The pair basis is the direct product of selected dressed eigenstates; the logical input states are \(|00\rangle,|01\rangle,|10\rangle,|11\rangle\).
- The target is the echo-containing logical operation stored in `src/nacs/gate/metrics.py`, rather than a bare iSWAP matrix.
- One nominal virtual-Z phase frame is used per protocol. It is not optimized individually for each noisy gate realization.

## Frozen artifacts and provenance

`benchmarks/frozen/gate_ready.npz` contains the same byte sequence as the previous v1.0 baseline. Its SHA-256 is

```text
e8c99f59e0cbab339f0f48ed23b21acba5fb81a86130f9297f89891b2aecdeb5
```

The other entries in `frozen/manifests/source_archives_sha256.json` provide checksums for chronological research ZIPs. The GitHub source tree is a flattened/renamed refactor and should not be represented as byte-identical to those archived modules.

Stored `benchmarks/data` CSV and JSON files are **archival numerical records**. They are immutable for comparison; the reproduction scripts put new outputs under `outputs/`.

## Basis-selection restrictions

Do not report the historical nine-state optimum as final. The canonical basis contains `aux51` and `aux79`. The 13-state and 15-state tests extend it, but this is a controlled subset of the full 128 single-molecule dressed levels.

The physical-parameter check uses a **global one-to-one overlap assignment**, rather than assuming identical eigenenergy indices, and aligns complex eigenvector phases. A small minimum overlap indicates the assignment requires scrutiny near a level mixing region; it is not proof of a spectroscopic level crossing.

## Sampling and uncertainty semantics

The 49-point stabilization matrix summarizes previously optimized controls under hypothetical **uniform** amplitude, detuning and exchange-error distributions. The source seeds and held-out datasets are preserved. Values reported as `p10` and `mean` describe those sampled distributions, not confidence in the underlying molecular constants.

The paired bootstrap resamples the archived *same* noise coordinates for competing protocols. Its interval is a measure of finite numerical ensemble sampling, not experimental uncertainty. The spectral/calibration uncertainty and hardware uncertainty are not covered.

## What this release does not reproduce automatically

- Every earlier local optimizer search across its original multi-archive research tree. Previously optimized control vectors are stored and directly propagated instead; reproducing searches may yield different local optima.
- Full 128×128 one-molecule states **for both** interacting molecules: the pair is truncated and checked at 11/13/15 selected levels per molecule.
- Physical finite-temperature motion, irreversible decoherence, loss, or actual apparatus noise.
- Experimental multi-power microwave resonance measurements (not present in the inspected public Dataverse files).

## Releasing a future version

Keep `benchmarks/frozen` immutable. Add new `v2` arrays and a new manifest if the spectroscopy changes. Create a separate benchmark directory and new test thresholds when controls or uncertainty assumptions change; do not overwrite historical CSV/JSON results to make the record look better.


## arXiv preparation — Step 2 screened-basis validation

See [`../reports/PAPER_STEP2_BASIS_CONVERGENCE.md`](../reports/PAPER_STEP2_BASIS_CONVERGENCE.md). The fast profile is for installation smoke tests; the frozen numerical comparisons were generated with `--profile full --steps 64`. Independent dipolar-screen, union and logical-map phase-alignment scripts and CSV/NPZ data are in `benchmarks/data/paper_step2/`. Never overwrite frozen benchmark data with reruns; the scripts default to ignored `outputs/` directories.

## arXiv preparation — Step 3 matched-control objectives

The [Step 3 technical report](../reports/PAPER_STEP3_MATCHED_OBJECTIVES.md) defines all control bounds, model assumptions, four scores, objective evaluation budgets, calibration conventions and error distributions. The frozen stage-by-stage search checkpoints, nominal results, paired 24-sample noise realizations and 25-state logical operators are at `benchmarks/data/paper_step3/`.

Commands (with `OPENBLAS_NUM_THREADS=1`):

```bash
python scripts/paper_step3_optimize.py --objective channel_guarded --seed 2026 --iterations 7 --population 30
python scripts/paper_step3_refine.py --objective channel_guarded --input outputs/paper_step3 --output outputs/paper_step3_refine/channel_guarded.json
python scripts/paper_step3_refine_11.py --objective channel_guarded --input outputs/paper_step3 --output outputs/paper_step3_refine_11/channel_guarded.json
python scripts/paper_step3_validate.py --steps 64 --output outputs/paper_step3_final_validation
python scripts/paper_step3_basis_columns.py
python scripts/paper_step3_noise.py --protocol channel_guarded --scenario tight --samples 24
```

The 25-state `column_propagator` is independently checked against the original dense pair propagator for two finite-ramp controls. The work uses the same frozen spectroscopy data and is *not* a calibration of real microwave hardware or full 128² basis proof.

### Refit the effective optical parameters

From an editable installation, run:

    OPENBLAS_NUM_THREADS=1 python -m nacs.spectroscopy.calibrate --output outputs/optical_fit

This recovers the fitter from NaCs_spectroscopy_frozen_v1.zip (SHA-256 in
frozen/manifests/source_archives_sha256.json). Only its obsolete follow-up
CLI message was removed. The unweighted fit uses three bright transition
offsets at 864.1 G and reports the 1->e line separately. The expected result is near -3526.155 kHz and 54.817 degrees.
Tests allow 0.002 kHz and 0.001 degrees respectively, and independently
require the three fitted residuals to agree with the static export within
0.1 kHz. Finite-difference fitting of MHz-scale eigenvalues can produce
small platform-dependent parameter differences. The fitter reports line centers
at full fitted precision. Table I uses the static seven-line export at the
rounded parameters in benchmarks/frozen/baseline_manifest.json. The common
reference correction enters the table comparison only.
