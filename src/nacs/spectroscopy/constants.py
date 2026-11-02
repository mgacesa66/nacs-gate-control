"""NaCs X(1)Sigma+, v=0 parameters (frequencies H/h in kHz).

The default *molecular* constants are the theoretical DFT values of Aldegunde
& Hutson, PRA 96, 042506 (2017), arXiv:1708.05734 Table I (23Na133Cs row).
These are not fitted to the Picard et al. (2025) gate experiment.

The optical tensor coefficient, centrifugal distortion, and rotational g
factor are NOT experimentally calibrated here; defaults omit their effects.
Never interpret the resulting spectroscopy as a fit to magic tweezers.
"""
from dataclasses import dataclass

# Correct unit: 7.622593285 MHz/T / (10^4 G/T) = 0.7622593285 kHz/G
MU_N_OVER_H_KHZ_PER_G = 0.7622593285

@dataclass(frozen=True)
class NaCsParameters:
    B_rot_khz: float = 1_740_000.0   # B_v = 1.740 GHz (theory, not spectroscopy fit)
    D_rot_khz: float = 0.0         # unspecified for this truncated N<=1 model
    eQq_na_khz: float = -97.0      # (eQq)_Na = -0.097 MHz
    eQq_cs_khz: float = +150.0     # (eQq)_Cs = +0.150 MHz
    c_na_khz: float = 0.0142       # 14.2 Hz
    c_cs_khz: float = 0.8545       # 854.5 Hz
    c3_khz: float = 0.1056         # 105.6 Hz tensor nuclear spin-spin
    c4_khz: float = 3.9418         # 3941.8 Hz scalar nuclear spin-spin
    g_na: float = 1.478
    g_cs: float = 0.738
    shielding_na: float = 639.2e-6
    shielding_cs: float = 6278.7e-6
    g_rot: float = 0.0             # NOT sourced/calibrated: omit rotational Zeeman
    field_G: float = 864.1        # Picard et al. experimental bias field
    # Effective optical tensor shift: coefficient u_t/h (kHz) in
    # H_trap/h = optical_scalar_khz * I + optical_tensor_khz * P_2(eps.n).
    # optical_tensor_khz=0 OMITs the crucial tweezer anisotropy; not an
    # accurate model of the experimental N=1 level structure.
    optical_tensor_khz: float = 0.0
    optical_scalar_khz: float = 0.0  # same scalar for all N: overall energy only
    # complex polarization vector in the laboratory xyz basis, normalized by code.
    polarization_xyz: tuple[complex, complex, complex] = (1.0+0j, 0j, 0j)
    dipole_debye: float = 4.6      # Picard et al. reported body-frame d

SOURCES = {
    "molecular_parameters": "https://arxiv.org/pdf/1708.05734",
    "gate_spectroscopy": "https://arxiv.org/pdf/2406.15345",
    "gate_nature": "https://www.nature.com/articles/s41586-024-08177-3",
}
