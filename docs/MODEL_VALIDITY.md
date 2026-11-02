# Model validity and release checklist

This page separates **numerical correctness**, **spectroscopic calibration**, and **experimental usefulness**. They are not interchangeable.

## 1. Canonical model and deprecated results

- [x] Preserved frozen `gate_ready.npz` and SHA-256 validation.
- [x] The public entry points default to the corrected eleven-state molecular basis with `aux51` and `aux79`.
- [x] The original nine-state calculation is excluded from final headline performance comparisons.
- [x] Eleven, thirteen and selected fifteen-state bases are propagated at identical fixed controls.
- [ ] Complete 128-state selected-dynamics convergence is **not established**; 128×128 one-molecule diagonalization does not imply 128² pair evolution has been converged.

At the corrected guarded controls, the nine-state model gives ~1.5% fidelity, whereas the 11-state model gives ~94.6%. This does **not** mean the earlier optimized nine-state result would have been ~1.5%: it means the old truncation is catastrophically inaccurate when evaluating the **same corrected control settings**. Do not confuse fixed-control convergence with comparing separate optimizations.

The `aux51` state has a strong microwave dipole and must be included even though it is not the nearest resonant state. Additional ranked states were selected from the existing frozen spectrum; all source indices and fidelity values are in the audit CSV.

## 2. Spectroscopy and independent data

- [x] 128-state hyperfine–Zeeman calculation reproduces the frozen export.
- [x] Seven measured microwave lines are listed with raw residuals.
- [x] A common rotational reference offset is documented separately.
- [x] Inferred power-dependent shifts are explicitly labeled **post hoc**, not independently validated.
- [ ] Finite-power line-center scans at ≥2 independently known microwave powers per hyperfine transition are unavailable in the inspected public data.
- [ ] Weak unfit rotational line discrepancy ~6.9 kHz remains.
- [ ] Trap tensor and microwave polarization at the molecules have not been independently calibrated.

The theoretical hyperfine constants were not refitted to explain the three ground-hyperfine discrepancies. Finite-power dressing is plausible but still a hypothesis.

## 3. Physical-parameter sensitivity

The release audit perturbs the actual 128-state Hamiltonian, re-diagonalizes it and tracks eigenstates by one-to-one overlap assignment. It then propagates the same guarded gate with absolute microwave carrier frequencies and field amplitudes held constant (no recalibration). Small illustrative changes give:

| Probe | Nominal guarded gate score |
|---|---:|
| Frozen nominal optical/field parameters | ~94.64% |
| Magnetic field −0.1 G / +0.1 G | ~91.82% / ~92.51% |
| Optical tensor −0.5% / +0.5% | ~92.37% / ~92.38% |
| Optical ellipticity −0.5° / +0.5° | ~33.21% / ~15.93% |

One enlarged-spectrum eigenstate overlap drops to ~0.62 for +0.5° ellipticity. This may indicate a near-degenerate labeling issue or substantial mixing; it requires investigation. We therefore **cannot claim model-parameter robustness** despite numerical agreement among several selected-state truncations.

These ranges are illustrative probes, not measured experimental uncertainties; fixed-carrier results are deliberately a stringent test and do not represent feedback-calibrated apparatus operation.

## 4. Experimentally available controls

- [x] Fixed microwave polarization is used for every pulse in Model A, including the echo.
- [x] Transfers and echo are included in end-to-end fidelity.
- [x] Microwave interaction remains active during access pulses.
- [x] Peak power and maximum total gate duration are constrained in archived optimization studies.
- [ ] Rotation of Ni's actual microwave antenna by 35° has not been demonstrated.
- [ ] Actual microwave field geometry/polarization purity and phase-transition impulse response are unknown.
- [ ] The precise experimental truncated-Gaussian IQ waveform has not been digitized.
- [ ] Molecular motion, finite temperature, dissipative loss and SPAM errors are excluded.

## 5. Uncertainty/statistical treatment

- [x] Each simulated protocol has one virtual-Z calibration shared across noise samples.
- [x] Comparisons use paired held-out noise draws in the archived Monte Carlo files.
- [x] The conditional bootstrap script states exactly which variability it measures.
- [ ] Noise amplitudes are hypothetical; the method cannot infer experimental uncertainties from them.
- [ ] Global optimality of the pulse parameters has not been demonstrated.
- [ ] A larger, multi-seed reoptimization in each enlarged basis is still necessary for a comparative performance claim.

## Release claim boundary

A technically reproducible GitHub report is supported. A quantitative high-fidelity experimental gate claim, a proven polarization sweet spot, or a definitive prediction of Ni's apparatus performance is **not** supported. The current interesting results are a **reduced-basis failure**, a **single-channel leakage-redistribution counterexample**, and a **conditional control-stability map**.
