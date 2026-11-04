# NaCs Gate Control

> ** `v0.2.0.dev3`:** Four microwave-control objectives are compared under the same peak-field, phase-ramp (2 µs), duration and candidate constraints. The isolated-transfer objective does not reliably determine complete-gate quality; the guarded channel achieves ~94.64% nominal fidelity in the corrected 11-state model, with ~94.63% in the 25-state cross-screen. None is robust against broad hypothetical control errors.

**Open research code for complete-operation errors in trapped polar-molecule qubits.**

This repository reproduces exploratory coherent-control calculations for individually trapped \(^{23}\mathrm{Na}^{133}\mathrm{Cs}\) molecules: dressed hyperfine/rotational structure, microwave access, resonant dipolar exchange, the rotational echo, and spectator-state leakage. It also preserves **negative results** from pulse shaping, phase programming, and robust optimization.

> **Status: research software / technical report, not an experimentally validated gate-design package.** The molecular spectrum uses theoretical constants and fitted effective tweezer parameters. Some experimental line frequencies require *inferred*, not independently measured, microwave dressing shifts for consistency. The control-noise envelopes are hypothetical. No experimental gate fidelity is predicted.

There are **no figures** in this initial release. Numerical CSV/JSON benchmarks and plotting-independent reproduction scripts are included.

## Why this repository exists

A molecular gate can look good when only the intended two-state interaction is optimized, yet perform badly as a complete **write into the interacting manifold → exchange → echo → retrieve** operation. Nearby molecular states are part of the dynamics even when far detuned, because their microwave dipoles can be much larger than the intended hyperfine-changing transition.

This project asks three reproducible questions:

- Does a reduced molecular basis predict the same complete logical channel after adding strongly coupled off-resonant states?
- Can a control pulse suppress one unwanted molecular state while accidentally transferring leakage into another?
- Under which **assumed** field-amplitude and residual-detuning uncertainties does an optimized gate approach its nominal coherent limit?

It does **not** claim a new high-fidelity experimental pulse sequence.

## Quick start

Python 3.10+ and NumPy/SciPy/SymPy are required.

```bash
python -m venv .venv
source .venv/bin/activate               # Windows: .venv\\Scripts\\activate
python -m pip install -e '.[dev]'
export OPENBLAS_NUM_THREADS=1          # useful for small dense propagators
python -m pytest -q
python scripts/quickstart.py
```

The last command recomputes complete-gate fidelity and leakage for three protocols in the **corrected 11-state-per-molecule model**. No optimization or model fitting is necessary for this smoke test.

To generate the independent release audit (CSV/JSON, **no plots**):

```bash
python scripts/release_audit.py --mode all --output outputs/release_audit
python scripts/uncertainty_summary.py --output outputs/uncertainty_summary.json
```

`--mode all` includes seven microscopic Hamiltonian diagonalizations and enlarged-basis gate checks. Outputs are placed in an ignored `outputs/` directory, leaving the stored benchmarks unchanged.

## Verified outcomes (conditional on this model)

| Study | Reproducible observation | Important qualification |
|---|---|---|
| Spectroscopy | Seven experimental lines compared; three hyperfine lines show substantial raw residuals | Their apparent reconciliation uses *one inferred microwave power per line*, not an independent test |
| Basis selection | The historical nine-state gate truncation omits a strongly coupled spectator, `aux51` | Earlier nine-state fidelity optimizations are **deprecated** |
| Coherent gate | A corrected 35° Gaussian-echo reference gives ~93.75% unconditional average fidelity | That antenna orientation has not been demonstrated on Ni's published apparatus |
| Echo shape | Square and truncated-Gaussian rotational echoes give very similar complete-gate results | No meaningful improvement from the Gaussian alone |
| Stage-resolved leakage | Final unwanted spectator population is ~3.48%; much arises during microwave access/retrieval | Rotational population in the *intended* interacting state is not itself leakage |
| Targeted phase control | Suppresses `aux51` ~2.72% → ~0.22% but raises other spectators ~0.76% → ~2.81% | Total gate fidelity **decreases**, ~93.75% → ~92.87% |
| Spectator-guarded control | ~94.64% nominal gate fidelity, ~4.69% total logical leakage | Local optimum; not a robust high-fidelity gate |
| Stabilization map | With assumed ±0.02% MW-amplitude and ±5 Hz detuning, locally adapted coherent mean reaches ~93.67% | Assumed ±1% exchange uncertainty; **not** a 99%-gate solution |

**Most important release warning:** The original nominal ~95% result from a nine-state model is **not** an accepted final result. It cannot be used as evidence of a physical advantage without the later omitted-spectator correction.

### What the independent release audit adds

The included [release audit](docs/REPRODUCIBILITY.md) recomputes four candidate gates in nine, eleven, thirteen and fifteen selected states. The `aux51` omission gives a large discrepancy at identical controls. The corrected eleven-to-fifteen-state comparison is close for the tested protocols, **but this is not convergence against every state in the 128-level molecule**.

A seven-scenario physical-parameter perturbation check also varies magnetic field by ±0.1 G, effective optical tensor by ±0.5%, and fitted optical ellipticity by ±0.5°. It holds the **absolute microwave carrier frequencies and electric-field amplitudes fixed**. This reveals strong model dependence, especially for optical ellipticity. These perturbations are *illustrative sensitivity probes*, **not** a published experimental uncertainty budget. See [Model validity](docs/MODEL_VALIDITY.md).

## Repository layout

```text
nacs-gate-control/
├── src/nacs/
│   ├── spectroscopy/              # 128-state Hamiltonian and microwave dipoles
│   ├── molecule/                  # molecular-model public entry points
│   ├── gate/                      # pair Hamiltonian, propagation, logical fidelity
│   ├── controls/                  # fixed polarization, Gaussian and phase programs
│   └── analysis/                  # stage-resolved leakage
├── benchmarks/
│   ├── frozen/                    # checked immutable v1.0 spectrum and matrices
│   └── data/                      # archived CSV/JSON comparisons (no plots)
├── frozen/manifests/              # SHA-256 provenance for original research ZIPs
├── scripts/                       # quickstart and release validation
├── tests/                         # regression, physical invariants, basis checks
├── reports/TECHNICAL_REPORT.md    # narrative technical report, not a journal paper
├── docs/                          # validity, assumptions and reproducibility
└── .github/workflows/tests.yml    # automated test workflow
```

`src/nacs` is a **dependency-flattened research refactor** of the prior frozen source archives, not a claim that each line is an identical copy. The original archive hashes are preserved in `frozen/manifests`; the exact frozen spectroscopy arrays are preserved and SHA-256 checked at load time. Historical experiments remain documented by their original versions, not silently promoted to final results.

## Running tests

```bash
python -m pytest -q                      # runs all tests, including ranked 15-state check
python -m pytest -q -m 'not slow'       # fast regression/invariant tests
```

Tests exercise actual Hamiltonian construction and full gate propagation; they are not just assertions that files exist. See [Reproducibility](docs/REPRODUCIBILITY.md) for hardware, random-seed, and model limitations.

## Scientific documentation

- [Technical report](reports/TECHNICAL_REPORT.md): methods, observations, negative results, what remains uncertain.
- [Model validity](docs/MODEL_VALIDITY.md): molecular calibration, state-selection, experimental control limitations and systematic error audit.
- [Reproducibility](docs/REPRODUCIBILITY.md): exact commands, fixed inputs, conventions, conditional statistical intervals, audit files.
- [Source provenance](frozen/manifests/source_archives_sha256.json): checksums of the chronological frozen research archives.

## Attribution and references

This is independent theoretical work inspired by the NaCs molecular-gate experiment; it is not an official release from Ni's laboratory nor related to it in any way. Experimental datasets and published theory remain the work of their respective authors.

- Picard *et al.*, **Nature** 637, 821–826 (2025), [DOI:10.1038/s41586-024-08177-3](https://doi.org/10.1038/s41586-024-08177-3).
- Aldegunde and Hutson, **Physical Review A** 96, 042506 (2017), [arXiv:1708.05734](https://arxiv.org/abs/1708.05734).
- Ni, Rosenband and Grimes, *Dipolar Exchange Quantum Logic Gate with Polar Molecules*, [arXiv:1805.10930](https://arxiv.org/abs/1805.10930).
- Experimental data associated with Picard *et al.*, [Harvard Dataverse DOI:10.7910/DVN/3UEBEV](https://doi.org/10.7910/DVN/3UEBEV).
- Bergonzoni, Jandura and Pupillo, **Physical Review A** 112, 032621 (2025), [DOI:10.1103/q9sd-rfp6](https://doi.org/10.1103/q9sd-rfp6).

## Release policy

- Data under `benchmarks/frozen` are **immutable**. Any recalibration must create a new version, not overwrite v1.0.
- Every reported gate comparison uses a clearly specified basis, logical fidelity definition and control set.
- A result that merely moves population between spectators is not recorded as a gate improvement.
- Sampling confidence intervals are conditional on simulated error distributions, **not** uncertainties in the physical Hamiltonian.
- Figures and a journal manuscript may be added later; neither is included in this release.

**License:** MIT, proposed for the research source code in this release. Review the author/license metadata before publishing publicly. This license does not override any licenses or rights attached to cited third-party experimental datasets or papers.

## Manuscript figures and spectroscopy table

Install plotting dependencies with `python -m pip install -e '.[figures]'`, then
run `python scripts/manuscript/reproduce.py` from the repository root.
See [Manuscript reproduction](docs/MANUSCRIPT_REPRODUCTION.md) for the five
figure inputs, output names, and the seven-line table consistency check.
