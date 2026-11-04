# Publication plan — Step 1: model-parameter sensitivity and recalibration

**Status:** computed and tested, 2026-11-02; this is a conditional theoretical sensitivity audit, **not** experimental calibration or a proof of broad Hamiltonian-parameter robustness. **No figures in this extension.** Full raw data and source are included.

## Question and scientific boundaries

The first GitHub release's physical-parameter audit varied the magnetic field and fitted tweezer parameters but held the *absolute* microwave frequencies and amplitudes fixed. It reported large apparent sensitivity, especially to the fitted ellipticity. This follow-up separates (1) frequency miscalibration, (2) intended-transition Rabi variation, and (3) changes in the actual multilevel spectator dynamics.

The six-step research publication plan is frozen in [`../PUBLICATION_PLAN_FROZEN.md`](../PUBLICATION_PLAN_FROZEN.md). The original public v0.1.0 archive and its frozen spectroscopy NPZ are not modified.

### Model and input assumptions

- Same 128-dimensional **single-molecule** static NaCs Hamiltonian and constants as v0.1.0, re-diagonalized at each parameter point. Eigenstates are assigned by global one-to-one maximum-overlap matching against the baseline, with complex phase gauge aligned to the original eigenvectors.
- Corrected selected **11-state-per-molecule** (121 pair states) coherent model; microwave and dipole–dipole interactions stay active during the square hyperfine transfers and Gaussian rotational echo.
- Fixed physical microwave polarization at the original hypothetical 35° orientation, fixed stored control durations, phase segments, dwell times, dipole separation, and echo Gaussian envelope.
- Four archived controls: `unmodulated`, `gate_only`, `aux51`-targeted, and `all_spectators`-guarded. **No gate reoptimization under perturbed Hamiltonians**; comparisons are of the SAME stored controllers under parameter change.
- Only *illustrative*, individually varied probes: magnetic field ±0.1 G; effective optical tensor ±0.5%; fitted optical ellipticity ±0.1°, ±0.25°, and ±0.5°. These are NOT experimental uncertainty intervals.
- The calculated molecular eigenvalue differences are used as *idealized tracked bare resonances*; actual microwave line centers include power-induced Stark shifts. This procedure is a theoretical calibration diagnostic, not a demonstrated laboratory calibration algorithm.

### Three recalibration modes

| Mode | Absolute carrier frequencies | MW electric field | Interpretation |
|---|---|---|---|
| `fixed` | Unchanged from the frozen baseline | Unchanged | Slow physical drift with **no** frequency feedback |
| `frequency` | Each transfer and echo carrier follows its tracked dressed transition, retaining the frozen control detuning | Unchanged | Optimistic two-frequency carrier tracking; not line-center power calibration |
| `frequency_rabi` | As above | Transfer and echo peak fields separately rescaled to retain the nominal *target* Rabi magnitudes | Optimistic measured Rabi recalibration; retains all other off-resonant dipole changes |

For each scenario, report both (a) **fixed nominal local-Z** fidelity and (b) a **refitted static local-Z** fidelity. Only the latter appears in the summary tables below. Locally recalibrating Z is not allowed separately for every random shot in a noise ensemble; here the model perturbation is a fixed hypothetical change. This audit does not include stochastic noise draws.

## 1. Calibration strongly changes the answer

The guarded controller's **unconditional average logical fidelity**, with scenario-specific static local-Z correction, is:

| Perturbation | Fixed carriers, fixed E field | Tracked carriers, fixed E field | Tracked carriers + target-Rabi calibration |
|---|---:|---:|---:|
| Frozen reference | 94.64% | 94.64% | 94.64% |
| B −0.1 G | 91.82% | 94.00% | 93.97% |
| B +0.1 G | 92.51% | 94.26% | 94.21% |
| Optical tensor −0.5% | 92.37% | 92.54% | 93.60% |
| Optical tensor +0.5% | 92.38% | 92.56% | 93.41% |
| Ellipticity −0.1° | 5.86% | 88.08% | 90.87% |
| Ellipticity +0.1° | 1.09% | 86.58% | 89.31% |
| Ellipticity −0.25° | 21.15% | 86.37% | 81.55% |
| Ellipticity +0.25° | 1.56% | 88.28% | 81.20% |
| Ellipticity −0.5° | 33.21% | 75.13% | 79.88% |
| Ellipticity +0.5° | 15.93% | 74.17% | 78.27% |

**Interpretation:** The catastrophic fixed-carrier ellipticity outcomes largely measure a deterministic microwave frequency mismatch, not an intrinsic 1–6%-fidelity limitation of the molecular operation. Following the *bare* transition centers restores much of the fidelity. But substantial residual multilevel error remains, especially for large ellipticity changes. The coherent response is **nonmonotonic**: separately restoring target Rabi frequencies can improve one perturbation and worsen another; an isolated transfer Rabi is not the sole correct calibration metric.

### Measured observables versus model predictions

For a ±0.1 G field change, the calculated transfer center moves about ∓0.225 kHz while the echo center barely moves. For a ±0.1° fitted optical-ellipticity change, **both** target transition centers move roughly ±3.47 kHz; at ±0.5° they move ±17.3 kHz. These are **theoretical undriven eigenvalue shifts**, not experimentally measured microwave line shifts.

More tellingly, the calculated relative `aux51`–`e` detuning changes approximately ∓6.94 kHz for ±0.1° ellipticity, ∓17.3 kHz for ±0.25°, and ∓34.7 kHz for ±0.5°. Correcting the `1`↔`e` carrier frequency therefore cannot hold *all* spectator detunings fixed. This is a physical mechanism for residual sensitivity even under ideal carrier tracking.

## 2. Phase calibration is a separate requirement

An important discrepancy shows up for ±0.1 G after retuning the carrier. For the guarded gate, the **locally refitted** Z phase gives 94.00% / 94.26% fidelity, whereas **reusing the nominal Z angles** gives only about 38.06% / 38.34%. This is a deterministic coherent-phase difference under the particular rotating-frame/target-unitary convention, not additional population loss: logical leakage is ~5.5%.

A frequency servo alone therefore does not guarantee the desired logical gate phase. A protocol needs an independently specified, experimentally feasible *phase calibration*. The refitted Z result is the more optimistic scenario, not a free improvement without laboratory calibration.

## 3. The leakage-redistribution claim must be narrowed

At the nominal frozen point:

| Protocol | Fidelity | `aux51` final occupation | Other spectator occupation |
|---|---:|---:|---:|
| Unmodulated | 93.75% | 2.72% | 0.76% |
| `aux51`-targeted | 92.87% | 0.22% | 2.81% |
| All-spectator guarded | 94.64% | 2.43% | 0.73% |

This is a valid *specific counterexample* to the idea that suppressing one spectator guarantees improvement in complete-channel fidelity. **It does not imply that the targeted pulse always reduces fidelity.** Under calibrated ellipticity −0.1°, the targeted pulse reaches ~93.26% while unmodulated reaches ~88.58%, with static local-Z recalibration. Under tensor ±0.5%, however, the targeted solution generally loses much of its suppression advantage and can become drastically worse.

Some parameter cases give `aux51` population *greater* than the unmodulated protocol. Thus the local optimization and cancellation mechanism are strongly model-dependent. We should formulate the paper around the robust logical statement **"single-channel suppression does not guarantee improved logical fidelity"**, not around a universal ranking between the specific stored pulse families.

## 4. What convergence has and has not been checked

- The new propagation is verified at the **complex amplitude** level against the archived baseline `gate_output` for all four stored controllers, within 2×10⁻⁹, and the nominal fidelity and leakage results are recovered.
- The 128-level Hamiltonian is re-diagonalized for each of the **11 probe points**; complex gauge is aligned and eigenvectors are matched one-to-one.
- For the *corrected 11-state* labels, the smallest matched overlap in these probes is ~0.9986. But if one includes four ranked 15-state candidates, the matched overlap at +0.5° ellipticity can be only **0.618**, indicating an ambiguity in extended spectator identity. The 11-state tracked result should therefore **not** be promoted to a fully converged or labeling-independent prediction.
- Complete paired-basis dynamics for **all** perturbed cases were **not** established here. This belongs to Step 2, including an avoided-crossing/continuation-label analysis. Do not infer 128² convergence from microscopic diagonalization.
- The dynamic Gaussian echo is midpoint integrated with 64 steps. The tests include unitary and partition invariants and exact nominal regression; more extensive time-step and basis convergence should accompany a paper figure.

## 5. Reproducibility

From repository root, Python 3.10+ with `numpy`, `scipy`, `sympy`, and `pytest`:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pip install -e '.[dev]'
python -m pytest -q
python scripts/paper_step1_parameter_recalibration.py --output outputs/paper_step1
```

Artifacts: [`gate_comparison.csv`](../benchmarks/data/paper_step1/gate_comparison.csv) (132 matched protocol cases), [`spectral_tracking.csv`](../benchmarks/data/paper_step1/spectral_tracking.csv) (11 microscopic probes), [`hypothesis_summary.csv`](../benchmarks/data/paper_step1/hypothesis_summary.csv) (model sensitivity of protocol ranking), [`run_metadata.json`](../benchmarks/data/paper_step1/run_metadata.json). These are frozen release data; reruns should go to ignored `outputs/` and not overwrite benchmark CSVs.

## Publication decision after Step 1

**Step 1 is complete as a theoretical comparative audit, not as experimental validation.** Recalibration removes most fixed-frequency errors, but the quantitative and comparative outcomes still vary with the fitted Hamiltonian parameters. The strong model dependence of the specific `aux51`-cancellation pulse actually reinforces the methodological motivation for *complete-channel* evaluation, while limiting what we may claim about that pulse's general performance.

**Next:** Step 2 must establish convergence with coupling-ranked additional states across the varied parameters, handle label mixing near avoided crossings, and check that the counterexample and recalibration conclusions survive selected-basis enlargement. No figures are created at this stage.
