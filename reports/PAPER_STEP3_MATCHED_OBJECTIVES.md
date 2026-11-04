# Step 3 — Matched microwave control objectives (numerical research report)

**Status:** completed exploratory, bounded-optimization comparison; 2026-11-02.  **Version:** `v0.2.0.dev3`.

**No new measured data or experimental validations.** All results use the frozen NaCs spectroscopy model, one fixed hypothetical microwave antenna orientation of 35 degrees, the same illustrative 1.9-µm molecule separation, and one optimized Gaussian-echo gate design family. The study does **not** show a new high-fidelity robust gate. Step 2's limitations on full-Hilbert-space convergence and Step 1's parameter-calibration uncertainty remain in force. The original GitHub release and Step 1/2 archives are immutable.

## Takeaway

Optimizing two isolated microwave transfers is an inadequate substitute for optimizing the complete, interacting gate. The best *isolated-transfer* candidate has a calculated forward/reverse single-molecule transfer probability of ~98.47% each, but only **82.31%** complete-gate fidelity in the corrected 11-state model, partly because its objective is blind to the interaction dwell and phases. Another independent optimization seed, with a barely worse isolated transfer (~98.44%), achieves ~93.66% complete-gate fidelity. The strong seed-to-seed spread illustrates the lack of complete-gate identifiability in an isolated objective; it should not be interpreted as a uniquely poor isolated optimum.

By contrast, the selected-`aux51`, logical-channel-only, and total-spectator-guarded objectives give *similar*, approximately **94.50–94.64%** nominal complete-gate fidelities on the 11-state basis, after 2-µs phase-transition smoothing. The guarded objective gives a modest nominal improvement but **no meaningful high-fidelity or broad-noise robustness breakthrough**. The earlier near-total suppression of `aux51` using idealized instantaneous phase jumps does not carry over to this physically smoother pulse family. This is an important *negative result*, not a contradiction of the prior phase-jump counterexample.

## 1. Fixed physical resources and candidate parametrization

All four optimization arms are allowed the identical six scalar control variables:

| Control | Bounds |
|---|---|
| Forward/reverse transfer duration, each | 0.206–0.288 ms |
| Free dipolar-evolution dwell | 0.450–0.720 ms |
| Transfer carrier detuning | 17.8–22.2 kHz, relative to frozen bare transition |
| Forward midpoint phase step | −π to +π rad |
| Reverse midpoint phase step | −π to +π rad |
| Reverse baseline phase | −0.5 to +0.5 rad |

Both transfers use two equal-duration phase segments with identical, **2-µs raised-cosine I/Q phase transitions** (six midpoint propagator slices per transition). The microwave **electric-field amplitude**, not the target Rabi frequency, is fixed across methods; the fixed 35° antenna polarization is shared across access, echo and retrieval. No instantaneous control jumps are allowed. The 0.3T-width truncated-Gaussian rotational echo is held fixed for all protocols; the same dipole–dipole Hamiltonian is active during every driven pulse. Maximum total operation duration is **1.30 ms**. Virtual local-Z corrections are optimized **once per nominal protocol** and remain fixed for held-out errors. All simulated nonlogical population counts as gate failure.

The frozen model's microwave polarization axis is *inferred*, not measured in the experiment. The 2-µs phase transition is an assumed bandwidth constraint, **not** a calibrated instrumental impulse response.

## 2. Objective definitions

Define unconditional projected logical-map fidelity `F`, final `aux51` occupation `P_aux`, total spectator occupation `P_spec` and isolated single-molecule forward/reverse transfer probabilities `P_f`, `P_r`. The predeclared objectives are:

| Arm | Minimized score | What it can miss |
|---|---|---|
| Isolated transfer | `1 - (P_f + P_r)/2` | Gate dwell, two-molecule exchange, echo, logical phases, spectators after retrieval |
| Selected spectator (`aux51`) | `1 - F + 2 P_aux` | Leakage redirected into **other** dressed spectator states |
| Complete logical channel | `1 - F` | Hardware/control robustness and explicit nonlogical-state penalty outside fidelity |
| Guarded logical channel | `1 - F + 2 P_spec` | Can trade coherent phase error against leakage and is not globally optimal |

The factor of two is a **declared heuristic weight**, not a derived physical constant. The selected-spectator arm retains a fidelity penalty to avoid the trivial solution of never performing the gate. The channel-only arm is a fourth reference, supplementing the three planned comparisons.

## 3. Search effort and safeguards

1. **Global search:** each arm receives the exact same four historical starting candidates, two predetermined random seeds (2026, 2027), two independent differential-evolution searches, **240 objective evaluations per seed** (**480 per arm**) on the *eight*-state search model. The training echo uses 16 Gaussian slices. No uncertainties are included in optimization.
2. **Exploratory local polishing on eight states:** each arm is allowed 240 additional Powell evaluations; runs terminated early after 160, 78, 134 and 59 calls (isolated, selected, logical-only and guarded respectively). These *differing actual counts must be disclosed*. Search-effort ceilings were equal, but **actual total function-call counts were not precisely equal**. Therefore, these are matched constrained *exploratory* searches, not proof of a strict equal-compute global-optimum comparison.
3. **Corrected-basis refinement:** one shared deterministic coordinate/random local-search schedule, the same seed (9184), and exactly **181 objective evaluations per arm** on the *11-state* canonical model (24 Gaussian midpoint slices); the fixed frozen control bounds and phase smoothing remain unchanged. All checkpoints and original candidates are preserved.
4. **Independent validation:** no new optimization at validation time. Nominal 11-state results use 64 Gaussian slices. Selected 15-state 32-slice and union-selected 25-state 24-slice computations test basis bias. The 25-state, 625-dimensional pair calculation uses Hermitian-eigenvector column propagation; direct matrix-exponential calculations for two cases agree within `2.7e-11` in all complex final amplitudes.
5. **Held-out error evaluation:** 24 shared quasi-static error draws per scenario, common to all protocols, generated with seed 32094. Once-per-protocol nominal local-Z calibration is reused for every draw. Noise is not measured experimental noise.

Search is nonconvex. None of the methods is certified globally optimal, and different penalty weights, reference hardware constraints or optimization starts could alter the rankings.

## 4. Nominal results, corrected 11-state model

All entries are averages over four logical computational inputs where applicable. These are **unconditional logical-channel** fidelities, not postselected conditional fidelities. Times include both transfers, complete dwell and fixed echo.

| Protocol | Fidelity | Total logical leakage | `aux51` | Other spectators (exclusive) | All spectators | Isolated transfer forward / reverse |
|---|---:|---:|---:|---:|---:|---:|
| Isolated-transfer winner | **82.3138%** | 12.0948% | 2.5862% | **7.6267%** | 10.2129% | **98.4697 / 98.4671%** |
| Isolated-transfer second seed | 93.6639% | 5.7664% | 2.7322% | 1.3240% | 4.0562% | 98.4388 / 98.4360% |
| Selected-`aux51` | **94.5219%** | 4.7304% | **2.3415%** | 0.7820% | **3.1235%** | 98.1656 / 98.3975% |
| Channel fidelity only | 94.5053% | 5.1091% | 2.4749% | 0.9227% | 3.3977% | 98.1602 / 98.1609% |
| **All-spectator guarded** | **94.6439%** | **4.6165%** | 2.4188% | **0.7282%** | 3.1471% | 98.3317 / 98.3975% |
| Frozen Gaussian reference | 93.7525% | 5.0489% | 2.7221% | 0.7594% | 3.4815% | 98.3774 / 98.3786% |

The isolated winner's higher standalone transfer probability does not imply a better complete gate. Note that its *dwell* and some phase degrees are invisible to the isolated objective, explaining much of the optimization degeneracy. The guarded objective nominally exceeds the unmodulated frozen Gaussian by **~0.89 percentage points**, but exceeds the best other full-gate arm by only ~0.12 percentage points. Do **not** claim a decisive advantage of guarded over all plausible globally optimized full-gate strategies.

The selected-`aux51` arm does reduce `aux51` relative to the guarded arm (~2.342% versus ~2.419%), while giving a slightly worse nominal gate. This is a **modest tradeoff**, not the extreme redistribution generated by the earlier discontinuous-phase pulse. It would be incorrect to present the previously reported ~92% `aux51` suppression as a result of this finite-bandwidth matched optimization.

## 5. Combined-screen 25-state check

The basis contains all 11 canonical dressed levels plus 14 extra states selected by the union of the independent Step 2 microwave- and dipolar-coupling rankings. Same pulse parameters and finite phase ramps; no retuning for 25 states.

| Protocol | 11-state nominal (64 slices) | 25-state (24 slices) | Difference |
|---|---:|---:|---:|
| Isolated-transfer winner | 82.3138% | **82.3299%** | +0.0162 pp |
| Isolated-transfer second seed | 93.6639% | **93.6441%** | −0.0198 pp |
| Selected `aux51` | 94.5219% | **94.4992%** | −0.0228 pp |
| Logical-channel only | 94.5053% | **94.4990%** | −0.0063 pp |
| **All-spectator guarded** | **94.6439%** | **94.6278%** | −0.0161 pp |
| Frozen Gaussian reference | 93.7525% | **93.7545%** | +0.0020 pp |

The complex projected logical maps were also compared at a common 24-slice Gaussian discretization. Raw 11-to-25-state Frobenius differences range from **0.0369 to 0.0390**, but fall to **0.00050–0.00086** after fitting only a physically irrelevant global phase and allowed static local-Z phases. All six complex matrices, fitted phase angles and operator differences are preserved in `basis/logical_matrices_11_25.npz` and `basis/operator_alignment_11_to_25.json`. This independently supports agreement of the *logical map*, not just the scalar fidelity; the model truncation is still not rigorously bounded.

Some differences combine basis enlargement with the smaller echo-discretization grid. There is no claimed exact 128x128 pair-space convergence. The key qualitative ranking survives the tested combined 25-state basis. All six 25-state complex logical matrices and Z corrections are supplied as machine-readable JSON.

## 6. Held-out illustrative control errors

All candidate pulses share the **same 24 random realizations per scenario**, with independent *uniform* quasi-static amplitude, microwave carrier detuning, and exchange fluctuations. All bounds are **hypothetical**. Fixed once-per-protocol virtual-Z correction; final nonlogical population contributes failure.

- *Tight:* microwave field amplitude ±0.5%; transfer carrier residual detuning ±50 Hz; dipolar coupling ±1%.
- *Broad:* microwave field amplitude ±2%; residual detuning ±500 Hz; dipolar coupling ±2%.

| Protocol | Tight mean | Tight 10th percentile | Broad mean | Broad 10th percentile |
|---|---:|---:|---:|---:|
| Isolated-transfer winner | 80.27% | 78.16% | 29.93% | 13.82% |
| Selected `aux51` | 92.13% | 90.51% | 35.28% | 15.50% |
| Logical-channel only | 92.08% | 90.83% | 34.73% | 15.07% |
| **All-spectator guarded** | **92.26%** | **90.85%** | **35.30%** | **15.65%** |
| Frozen Gaussian reference | 91.37% | 90.03% | 34.33% | 15.29% |

The guarded-minus-channel-only mean in this **specific sample** is +0.176 pp (tight) and +0.574 pp (broad). Paired draw-bootstrap 95% intervals are [+0.064,+0.286] and [+0.440,+0.717] pp, respectively. These intervals describe *sampling uncertainty conditional on the frozen Hamiltonian and chosen uniform error model*. They are **not** uncertainty in the physical molecular parameters, are based on only 24 Monte Carlo draws, and are inadequate for claiming an experimentally robust performance advantage. In particular, a broad mean near 35% is unsuitable for a useful quantum-gate implementation.

## 7. Interpretation and release limitations

- The failure of isolated transfer as a surrogate complete-gate objective is reproduced under matched finite-bandwidth controls. The extreme isolated-winner-versus-second-seed difference reflects unconstrained dwell/phase degrees **as well as** genuine interacting dynamics; do not present it as proof of unavoidable leakage.
- With phase ramps enforced, optimizing `aux51` alone no longer produces the earlier large cancellation. Optimized spectator and complete-gate objectives yield **similar** performance. The strong earlier leakage-redistribution counterexample is conditional on phase jumps and should be separately identified.
- Ranking robustness under the 25-state union does **not** validate all 128 single-molecule states, motional wavefunctions, microwave polarization calibration, or irreversible loss.
- The model's effective optical ellipticity is highly sensitive (Step 1). The present nominal results and uncertainty maps do **not** average over uncertain molecular Hamiltonians. Carrier and Rabi calibration can mitigate but not eliminate that issue.
- Pulse bandwidth is represented by a 2-µs raised-cosine *assumption*, not measured hardware; microwave signal distortions, phase noise and correlated errors are omitted.
- The 8-state training model is not quantitatively trustworthy on its own. The recorded independent 11-state search and 25-state validation are essential. Neither selected-basis size nor finite search budget certifies a global optimal solution.

## 8. Reproduce and extend

Install and run from the repository root, ideally in a clean virtual environment:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pip install -e '.[dev]'
python -m pytest -q

# Each global checkpoint is independent and skipped if already present:
python scripts/paper_step3_optimize.py --objective isolated_transfer --seed 2026 --iterations 7 --population 30
python scripts/paper_step3_optimize.py --objective selected_aux51 --seed 2026 --iterations 7 --population 30
python scripts/paper_step3_optimize.py --objective channel_fidelity --seed 2026 --iterations 7 --population 30
python scripts/paper_step3_optimize.py --objective channel_guarded --seed 2026 --iterations 7 --population 30

# Repeat with seed 2027, run optional 8-state Powell refinements, then:
python scripts/paper_step3_refine_11.py --objective channel_guarded \
  --input outputs/paper_step3 --output outputs/paper_step3_refine_11/channel_guarded.json
python scripts/paper_step3_validate.py --steps 64 --output outputs/paper_step3_final_validation
python scripts/paper_step3_basis_columns.py
python scripts/paper_step3_noise.py --protocol channel_guarded --scenario tight
```

Baseline machine-readable output is committed in `benchmarks/data/paper_step3/`. Full checkpoint data include both search seeds, per-arm 8-state and 11-state refinements, nominal model metrics, all six 25-state complex logical maps, and 240 paired held-out noise rows.

**Source integrity:** The parent Step 2 release ZIP is retained byte-for-byte with SHA-256 manifest; the frozen source and previous publication plan are not modified. This Step 3 extension contains no figures.

**Next publication step:** Step 4, literature and novelty audit. It should ask whether the *mechanism* (especially the original strong `aux51` cancellation that fails under finite ramps) provides novelty beyond existing molecular robust optimal-control studies. These data alone do **not** support a publication claim of a universally superior pulse family.
