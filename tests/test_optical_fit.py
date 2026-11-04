import csv
import json

import numpy as np

from nacs.data import DATA
from nacs.spectroscopy.calibrate import fit_three_bright_lines


def test_refit_recovers_archived_optical_parameters():
    reference = json.loads((DATA / "frozen" / "baseline_manifest.json").read_text())
    fit = fit_three_bright_lines()
    assert fit.success
    # Numerical derivatives of MHz-scale eigenvalues vary between LAPACK builds.
    # Check parameters with explicit units, and independently constrain the
    # bright-line residuals to 0.1 kHz (the smallest reported line uncertainty).
    assert abs(fit.x[0] - reference["optical_tensor_khz"]) < 0.002  # kHz
    assert abs(fit.x[1] - reference["optical_ellipticity_deg"]) < 0.001  # degrees
    with (DATA / "data/spectroscopy/seven_line_accounting.csv").open() as source:
        rows = {row["line"]: row for row in csv.DictReader(source)}
    expected = [float(rows[name]["raw_residual_khz"])
                for name in ("0->e", "0->f_z", "0->f_y")]
    np.testing.assert_allclose(fit.fun, expected, rtol=0, atol=0.1)
