# Reproducing the manuscript figures and spectroscopy table

From the repository root, install the optional plotting dependencies and run:

```bash
python -m pip install -e '.[figures]'
python scripts/manuscript/reproduce.py
```

The command validates the seven table rows against the canonical static CSV,
then creates five vector PDFs under outputs/manuscript. The figures use stored
numerical results and the two schematics use illustrative geometry. Expected
runtime is a few seconds, depending on plotting-library startup.

| Manuscript figure | PDF filename | Source |
|---|---|---|
| 1, model workflow | fig1_model_chain.pdf | Schematic |
| 2, molecular state spaces | fig5_level_schematic.pdf | Schematic |
| 3, ellipticity recalibration | fig2_ellipticity_recalibration.pdf | benchmarks/data/paper_step1/gate_comparison.csv |
| 4, basis convergence | fig3_nonmonotonic_basis.pdf | benchmarks/data/paper_step2/paper_step2_extended/extended_nominal.csv and paper_step2_union/cross_screen.csv |
| 5, objective comparison | fig4_objective_nonidentifiability.pdf | benchmarks/data/paper_step3/nominal/nominal_11.json |

PDF filenames preserve their original source names. Manuscript ordering places
the molecular schematic second. Figure 1 uses the corrected box layout and
Figure 5 uses the corrected separated labels with leader lines.

For a different output directory:

```bash
python scripts/manuscript/reproduce.py --output /absolute/path/to/figures
```

The spectroscopy table fragment is manuscript_tables/spectroscopy.tex. It
contains the seven rows from the reviewed paper. Its residual columns use the
static export and the common 5.88 kHz reference correction. The checker reads
all seven rows and detects disagreement with the canonical CSV.

To check a separate manuscript source containing the same table:

```bash
python scripts/manuscript/check_seven_line_table.py --table /path/to/manuscript.tex
```

The plotting workflow is designed for a repository checkout with its benchmark
data. Generated PDFs may vary in font embedding and metadata across Matplotlib
versions. Numerical inputs are the benchmark files distributed with the code.
