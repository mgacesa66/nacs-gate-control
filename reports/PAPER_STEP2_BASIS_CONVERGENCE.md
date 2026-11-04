# arXiv publication plan: Step 2 -- systematic basis enlargement

**Status (2026-11-02): screened selected-state convergence audit completed.** This is not a proof of full 128 x 128 two-molecule convergence. No figures are included. The frozen six-step plan and Step 1 remain unchanged.

## Key results

Our 11-state molecular basis gives a nominal all-spectator-guarded gate fidelity of **94.6380%**. With 25 molecular eigenstates per molecule (625 pair states), selected by the *union* of two independent microwave and dipolar coupling screens, it is **94.6372%**. The `aux51`-targeted control still has worse nominal complete-gate fidelity, even though it strongly suppresses that one spectator. Thus the qualitative leakage-redistribution counterexample survives.

**The important caveat: convergence is nonmonotonic.** In the microwave-selected sequence, adding intermediate states 49, 9, 63 and 44 lowers nominal guarded fidelity to ~94.3585% at 19 states. Subsequent states restore it to ~94.637% at 21-25 states. A simple 11-versus-15-state agreement would conceal this transient basis dependence. Agreement of final scalar fidelities is also not automatically agreement of the full complex logical operator.

## Step 2A. Physics, frozen inputs and selection

* Same fixed 128-state single-molecule NaCs hyperfine-Zeeman model and dipole matrices as the published v0.1.0 GitHub package; no new spectroscopy fits.
* Both microwave and resonant dipole-dipole interactions are retained in every propagated selected pair basis; the same stored complete-gate pulses are used, without gate-level reoptimization.
* The selected single-molecule space is enlarged from 11 to at most 25 states, implying 121 to 625 pair states. The whole 16,384-state pair Hilbert space is **not** propagated.
* Nominal Gaussian echo uses 64 midpoint intervals, while transfer phases follow the archived piecewise-constant controls. Model parameter probes use their previously defined *hypothetical* variations, with distinct fixed-carrier, carrier-tracked and frequency-plus-Rabi calibration cases.
* The model tracks all 128 single-molecule eigenvectors by global maximum-overlap assignment. The labels are conventions; they are not physical invariants when eigenvectors mix near an avoided crossing.

Two independent candidate rankings reduce selection bias:

1. **Microwave coupling screen.** Rank omitted ground and rotational eigenstates by the largest `drive coupling / max(carrier detuning, 10 kHz)` across the transfer and echo electric fields and the core 0, 1 and e states. Use all 11 molecular-parameter probes. Main omitted candidates: original indices 38, 69, 40, 72, 49, 9, 63, 44, 50 and 48.
2. **Pair-dipolar screen.** Rank omitted states by the largest excitation-conserving `|<j,c|Hdd|a,b>| / max(pair energy mismatch, 0.1 kHz)` with a,b,c drawn from the corrected active 11-state basis. This independent metric prioritizes original indices 34, 2, 48, 32, 40, 38 and 57, some of which are not prominent in the microwave-only screen.

**Union check:** The first ten microwave-ranked omissions plus previously unselected states among the eight strongest dipolar candidates add 14 unique states to the original 11. The combined space has 25 states per molecule. The screening heuristics are not certified truncation-error bounds: they do not include all pulse phase-step sidebands, multistep virtual pathways, molecular motion or a complete spectral bandwidth analysis.

## Step 2B. Nominal fidelity versus selected basis

All entries below are percentages, with independently optimized **static local Z corrections** for each fixed basis (not shot-by-shot Z correction). Controls are unchanged between rows. The final row uses the independently cross-screened basis.

| States per molecule | Unmodulated | aux51-targeted | All-spectator-guarded |
|---|---:|---:|---:|
| 11 | 93.75245 | 92.86550 | 94.63802 |
| 15, microwave | 93.75096 | 92.85050 | 94.63629 |
| 17, microwave | 93.60175 | 92.70244 | 94.46688 |
| 19, microwave | 93.49951 | 92.58228 | 94.35854 |
| 21, microwave | 93.74996 | 92.84836 | 94.63715 |
| 23, microwave | 93.75010 | 92.85030 | 94.63735 |
| 25, microwave | 93.75000 | 92.85019 | 94.63742 |
| **25, microwave + dipolar** | **93.74987** | **92.84778** | **94.63725** |

The corrected 11-state and cross-screened 25-state results differ by about **-0.0008 percentage points** for guarded fidelity. The targeted control is more basis-sensitive, but its nominal fidelity stays substantially below the unmodulated and guarded choices. Total logical leakage and final aux51/spectator populations at the cross-screen endpoint also change only slightly.

These are *specific model benchmarks*, not a universal statement that spectator targeting never works. Step 1 demonstrated that protocol rankings can change after modifying and recalibrating the molecular Hamiltonian.

## Step 2C. Parameter perturbations and carrier recalibration

The full screened run includes all 11 Step 1 molecular-parameter configurations, with 11/15 states, plus higher-state comparisons for ellipticity and multiple calibration modes. The 21/23-state check includes the four ellipticity probes of +/-0.1 and +/-0.5 degrees. The 25-state microwave/dipolar cross-screen also includes +0.5 degrees.

**All-spectator-guarded fidelity with ideal frequency and target-Rabi calibration:**

| Illustrative ellipticity change | 11 states | Enlarged result | Change (percentage points) |
|---|---:|---:|---:|
| -0.1 degrees | 90.8686% | 90.8622% (23) | -0.0064 |
| +0.1 degrees | 89.3083% | 89.3135% (23) | +0.0052 |
| -0.5 degrees | 79.8845% | 79.7551% (23) | -0.1294 |
| +0.5 degrees | 78.2702% | **78.4418%** (25, union) | +0.1715 |

The impact of intentional carrier retuning in Step 1 is substantially larger than these final selected-basis fidelity corrections. The calibration conclusion remains: correcting the wanted transition frequency can recover substantial apparent infidelity but cannot restore every spectator detuning.

At +0.5 degrees, the old spectral eigenstates 69 and 72 mix, and the best individual overlap is only ~0.618. We therefore retain **both states together** rather than relying on a physical interpretation of either isolated label.

## Step 2D. What remains uncertain

- **Full Hilbert-space convergence is not proven.** We use complete 128-state microscopic eigenvectors for ranking and tracking, but selected two-molecule propagation only. More distant many-step, multi-photon and dipolar channels could still matter.
- **Scalar logical fidelity is not full channel convergence.** The raw complex 4 x 4 logical map varies more than the local-Z-optimized scalar fidelity, but most of the change is a **global phase**. In a 32-step supplementary check, the nominal guarded 11-to-25-state map has raw Frobenius difference **0.0384**, reduced to **0.00075** after individually allowed local-Z calibration and a common global phase. At +0.5 degrees, the corresponding difference is **0.0449**, reduced to **0.00316** after alignment. The latter remains a nonzero physical-channel correction; these are not full-Hilbert-space bounds. Both raw and phase-aligned comparisons are exported with the corresponding complex 4 x 4 matrices.
- **One-body identity near crossings is ambiguous.** An overlap-based permutation gives unique indices but does not guarantee that an individual dressed level is a smooth eigenvector branch.
- **Hardware validation is incomplete.** Microwave polarization, experimental power-dependent line centers and exact phase-transition bandwidth are not independently established. Motional degrees and irreversible loss are omitted.
- **No new optimizations were performed.** Equal-budget, bandwidth-limited comparisons remain publication Step 3.

**Publication implication:** The central qualitative results survive all tested enlargements, including the two independent spectator screens. This supports using the corrected numerical model for a preprint *with stated convergence limitations*. It does **not** justify labeling the original 11-state or the present 25-state model as exactly converged to the full NaCs pair Hilbert space.

## Reproducibility

From the repository root:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pip install -e '.[dev]'
python -m pytest -q
python scripts/paper_step2_basis_convergence.py --profile fast --steps 32 --output outputs/step2_fast
python scripts/paper_step2_basis_convergence.py --profile full --steps 64 --output outputs/step2_full
python scripts/paper_step2_extended.py --steps 64 --output outputs/step2_extended
python scripts/paper_step2_perturbed_highbasis.py --steps 64 --output outputs/step2_perturbed
python scripts/paper_step2_dipolar_screen.py --steps 48 --output outputs/step2_dipolar
python scripts/paper_step2_union_check.py --steps 64 --output outputs/step2_union
python scripts/paper_step2_operator_alignment.py
```

Frozen output CSV/JSON/NPZ files are in `benchmarks/data/paper_step2/`; reruns go into ignored `outputs/` by default. Inputs remain the prior frozen spectroscopy, v0.1.0, and Step 1 release. No figures are included.
