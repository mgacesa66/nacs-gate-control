import json

from nacs.data import DATA
from nacs.spectroscopy.calibrate import fit_three_bright_lines


def test_refit_recovers_archived_optical_parameters():
    reference = json.loads((DATA / "frozen" / "baseline_manifest.json").read_text())
    fit = fit_three_bright_lines()
    assert fit.success
    assert round(fit.x[0], 3) == reference["optical_tensor_khz"]
    assert round(fit.x[1], 3) == reference["optical_ellipticity_deg"]
