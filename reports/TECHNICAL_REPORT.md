# Technical report: Complete-gate errors in NaCs molecular qubits

**Version:** 0.1.0 (2026-11-01)  
**Scope:** model and software validation; **not** a journal manuscript or an experimental claim.  
**Model:** frozen NaCs spectroscopy v1.0; two trapped molecules; stationary sites; coherent dynamics; same fixed microwave polarization for transfer, echo and retrieval (Model A).

## Executive summary

This repository contains a reproducible numerical investigation of *complete logical gate* errors, rather than isolated microwave transition fidelities. The main observation is a failure of optimizing one part of a multilevel gate in isolation: a pulse can almost eliminate population in one off-resonant rotational state and yet reduce complete-gate fidelity by redirecting population to other spectators. A second important outcome is a basis-selection failure: comparing two small molecular bases is not enough to establish convergence if both exclude a strongly coupled state.

The work includes both positive and negative outcomes. In particular, a smooth Gaussian rotational echo performs almost identically to a square echo when the full operation is reoptimized. Simple amplitude shaping and restricted phase programming did **not** produce a reliable high-fidelity gate across the previously assumed broad error range. Lower noise improves the gate, but exposes a coherent-error floor around 93–94% for the present pulse family and calibrated assumptions.

This is a **model-specific numerical result**. We have not verified experimental microwave polarization, linewidth, trap tensor shift or noise distribution to the accuracy needed for quantitative apparatus predictions.

## Computational approach

### 1. Single-molecule states and spectroscopy

- Molecular basis: \(|N,m_N,m_{I,\mathrm{Na}},m_{I,\mathrm{Cs}}\rangle\), with \(N=0,1\), \(I_\mathrm{Na}=3/2\), \(I_\mathrm{Cs}=7/2\). There are 128 basis states.
- Terms: rotation, hyperfine quadrupole, nuclear Zeeman, scalar/tensor nuclear spin–spin, spin–rotation, and an effective tensor optical shift.
- Scalar/molecular hyperfine constants are the theoretical NaCs values attributed to Aldegunde and Hutson. An effective optical tensor magnitude and ellipticity were fitted against three selected rotational lines. These are *not* independently calibrated trap settings.
- Diagonalization returns complex dressed microwave dipole matrices in the same eigenvector gauge. The gate uses the exact immutable `benchmarks/frozen/gate_ready.npz` export (SHA-256 checked).

Seven published spectroscopy transitions are included in `benchmarks/data/spectroscopy`. The three fitted bright rotational splittings agree internally after a common ~5.88 kHz offset. The measured hyperfine-changing lines disagree with undriven eigenvalue differences by ~11–27 kHz. Off-resonant microwave dressing **can** generate shifts of suitable signs and magnitudes, but those powers were inferred individually using the observed residuals. The corrected agreement is therefore *not independent prediction*. One unfit weak rotational line still misses by ~6.87 kHz after the common offset. The accessible Harvard Dataverse archive includes a hyperfine-transfer Rabi trace, but no verified multi-power line-center series for the disputed transitions.

### 2. Two-molecule dynamics

The pair Hamiltonian is built from the same complex dressed dipoles using the electric dipole–dipole tensor. Resonant, rotational-excitation-number-conserving terms are retained under the pair RWA, with the molecular separation set by the archived 1.9 μm model. The gate contains transfer, dwell, rotational echo, dwell and reverse transfer. Microwave–dipolar coupling stays active during each finite control pulse.

The projected logical operator \(M=P_\mathrm{logical} U P_\mathrm{logical}\) is generally nonunitary because population can remain in nonlogical molecular levels. For four logical input states, the unconditional average fidelity is

\[
F_{\mathrm{avg}}(M,U_{\mathrm{tar}})=\frac{\operatorname{Tr}(M^\dagger M)+|\operatorname{Tr}(U_{\mathrm{tar}}^\dagger M)|^2}{20},
\]

with a single nominal output virtual-\(Z\) calibration allowed for each protocol. **No shot-by-shot phase recalibration** is permitted in robust comparisons. Logical leakage is \(1-\operatorname{Tr}(M^\dagger M)/4\). The intended \(|e\rangle\) manifold is not counted as a spectator while operating, though residual occupation there *at the end* is logical leakage.

### 3. Basis selection

Canonical gate calculations use eleven molecular states per molecule: the original nine labeled states plus `aux51` and `aux79`. The 11-state gate is compared at identical controls with 13-state and selected 15-state enlargements. A historical nine-state run is retained solely as a counterexample to naive convergence; it cannot be used as a final control benchmark.

The independent release audit propagates the same four stored protocols across four bases, and reveals very large errors in the nine-state truncation at the corrected controls. The maximal 11-to-ranked-15 change in unconditional fidelity across the four checked controls is approximately **0.00022 in absolute fidelity** (about 0.022 percentage points). This is useful evidence of local convergence for this particular added-state list, **not** complete 128-state convergence.

## Reproduced numerical observations

### Echo envelope: a negative result

At the fixed 35° orientation, a joint Gaussian echo gives approximately 93.75% nominal unconditional fidelity, compared with approximately 93.68% for the archived optimized square echo. The difference (~0.075 percentage points) is not a meaningful performance improvement under the present physical-model uncertainty. The Gaussian echo's final average logical leakage is ~5.05%, versus ~5.36% for the square echo. This is a near-equivalence result, not an endorsement of more elaborate pulses.

### Stage-resolved leakage

The Gaussian reference has ~1.60% unwanted spectator occupation after the first transfer, ~2.05% after the echo and ~3.48% after retrieval, averaged over four computational inputs. Free-evolution intervals produce little **net** growth in the spectator sector; phases accumulated there can still affect later retrieval. Temporary spectator occupation during microwave pulses reaches nearly 10% on average but largely returns coherently. The final ~5.05% logical leakage is split into ~3.48% spectators and ~1.57% intended-rotational states that did not return to logical space.

### Single-state suppression redistributes leakage

| Protocol | `aux51` occupation | Other spectator occupation | Overall logical fidelity |
|---|---:|---:|---:|
| Unmodulated Gaussian reference | 2.72% | 0.76% | 93.75% |
| `aux51`-targeted phase control | 0.22% | 2.81% | 92.87% |
| Guarded, all-spectator objective | 2.43% | 0.73% | 94.64% |

The aggressive `aux51` pulse suppresses that particular state by ~92%, but redirects excitation into others, particularly `f_y` and `e_prime` in the selected basis. The total complete-gate fidelity becomes **worse**. A constrained all-spectator penalty can give a modest nominal improvement but does not resolve the broad-error robustness problem. A large, fast phase step can be sensitive to finite microwave I/Q transition times; retaining 5 μs transitions changes the outcome of the aggressive control substantially.

These are local optimization outputs, not certificates of global control optimality.

### Microwave stability map

A 7×7 matrix of amplitude and independent detuning uncertainty envelopes was previously computed with fixed dipolar uncertainty (±1%), and its original numerical values are retained in `benchmarks/data/stability/map_49_cells.csv`. The locally reselected 35° gate achieved about 93.67% mean under ±0.02% amplitude and ±5 Hz residual detuning on a separate 120-shot test. At ±2% and ±500 Hz, the corresponding complete-gate score remained poor (about 40.4% on that 120-shot set). A single carrier recalibration after a slow ±1% amplitude setpoint shift recovered ~92–93% in the selected simulated cases.

These are specifications **conditional on hypothetical uniform noise distributions**; the source does not establish Ni's laboratory fluctuation amplitudes. Calibrating the carrier after slow amplitude drift is not the same as correcting stochastic detuning on each gate shot.

## New release-specific validation

The release audit in `scripts/release_audit.py` independently recomputes the canonical 11-/13-/15-state comparison rather than simply copying old fidelity numbers. It also recomputes seven physical-parameter perturbations:

- magnetic field ±0.1 G;
- effective optical tensor amplitude ±0.5%;
- effective optical ellipticity ±0.5°.

The calculation tracks dressed eigenstates by a one-to-one, maximum-overlap assignment and aligns their complex phases before evaluating the gate. **Absolute** transfer and echo carrier frequencies and microwave electric fields remain fixed; no refitting or recalibration is permitted in this audit.

The nominal guarded gate is ~94.64% in the 11-state model. With ±0.1 G changes, its fixed-carrier fidelity drops to roughly 91.8–92.5%; with ±0.5% tensor changes, to about 92.4%. In the ±0.5° ellipticity tests, the score drops much further (~16–33%); one selected level's overlap with its original eigenstate can fall to ~0.62, indicating that label tracking itself is fragile near a mixing region. This warns strongly against using the model to prescribe hardware parameters without improved optical calibration and further eigenstate-assignment checks.

**This sensitivity calculation does not provide an experimental uncertainty interval.** The ± ranges are deliberately selected probes, not measured parameter distributions; their effect is evaluated without carrier recalibration or pulse reoptimization.

## Reproducibility and statistical caveats

- Use `python scripts/quickstart.py` for nominal gate benchmarks; `scripts/release_audit.py` for spectral, basis and physical-parameter checks; `scripts/uncertainty_summary.py` for paired uncertainty summaries of archived Monte Carlo shots.
- Monte Carlo bootstraps resample **simulated noise realizations**, not alternative molecular Hamiltonians. A 95% interval from this procedure does not capture model uncertainty or experimental repeatability.
- Selected-state convergence is not a guarantee of full Hilbert-space convergence, and a local optimizer does not establish the best achievable gate.
- The full calculations omit motional wavepackets, dephasing, spontaneous loss and SPAM errors. Apparatus hardware availability and polarization controllability were not independently verified.

## What is ready to reuse

The source exports the frozen NaCs microscopic energies/dipoles, selected pair Hamiltonians, exact square-pulse and time-sliced Gaussian propagation, a complete logical-fidelity metric, and stage-resolved, per-spectator population accounting. It is suitable for independent theoretical exploration and for extending a research statement. It is **not yet** a quantitative experimental design tool or a claim of PRA-level model validation.

## Next validation priorities

1. Correct and independently calibrate the optical polarization and effective trap tensor, including level-tracking near avoided crossings.
2. Check strong spectator selection across the full 128-level microwave spectrum for multiple polarizations and pulse bandwidths, then test larger pair bases.
3. Reevaluate historical square, smooth and composite control objectives in a single corrected basis and a truly matched control budget.
4. Calibrate the physical detuning/amplitude noise distributions with experimental data if available.
5. Add molecular motion and loss before interpreting absolute measured gate fidelities.

## Primary sources

See the linked Picard, Aldegunde–Hutson, Ni–Rosenband–Grimes and Bergonzoni references in the repository [README](../README.md). The raw Harvard Dataverse archive is not redistributed in this release.
