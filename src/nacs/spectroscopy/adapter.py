"""Integration interface for the *future* gate simulator (no time evolution here).

The module is intentionally independent of nacs_measured_spectroscopy.py.
All H/h values are CYCLIC kHz. Convert to angular rad/ms by multiplying 2*pi.

Microwave complex amplitudes are meaningful only in the SAME eigenvector
phase (gauge) basis as the exported dressed states. The diagonalizer fixes
largest bare-basis component real-positive for reproducibility.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from .validation import polarization_xyz_coefficients


def load_export(folder):
    """Read generated npz and state_index.json without constructing model."""
    root=Path(folder)
    data=np.load(root/'nacs_dressed.npz')
    with open(root/'state_index.json',encoding='utf-8') as f:states=json.load(f)
    return data,states


def single_molecule_rwa_blocks(data, state_indices, frequency_khz,
                              omega_reference_khz, reference_pair,
                              microwave_axis='x'):
    """Return (H_static/h, H_drive/h) for future Schrödinger propagation.

    state_indices: dressed eigenindices included in reduced molecular basis.
    frequency_khz: microwave carrier (positive cyclic kHz, ~3.48e6).
    omega_reference_khz: resonant full Rabi frequency Ω_ref/(2π) in kHz.
    reference_pair: (ground_index, excited_index) for Ω_ref calibration.
    H_on/h = H_static/h + H_drive/h.
    H_off/h = H_static/h; a microwave phase phi changes drive's upper
        triangular couplings by exp(i phi), and Hermitian conjugates accordingly.

    This uses a rotating-wave approximation; states with N~0 have zero
    rotating-frame excitation number, N~1 have one. Excited-state identity
    follows from the *dominant N weight* rather than energy ordering.
    """
    inds=np.asarray(state_indices,int)
    E=np.asarray(data['energies_khz'])
    vec=np.asarray(data['vectors'])
    bare=np.asarray(data['bare_basis_quantum_numbers'])
    isexc=[]
    for i in inds:
        n1=float(np.sum(np.abs(vec[bare[:,0]==1,i])**2))
        if abs(n1-0.0)>.1 and abs(n1-1.0)>.1:
            raise ValueError(f'Rotational parity is strongly mixed for state {i} (N1 fraction={n1})')
        isexc.append(int(n1>.5))
    isexc=np.asarray(isexc)
    carrier_diag=E[inds]-frequency_khz*isexc
    carrier_diag-=np.min(carrier_diag)  # harmless global phase
    hstatic=np.diag(carrier_diag).astype(complex)

    coeff=polarization_xyz_coefficients(microwave_axis)
    D=sum(coeff[q]*data[{-1:'d_q_minus1_debye',0:'d_q_zero_debye',1:'d_q_plus1_debye'}[q]]
          for q in (-1,0,1))
    ground,excited=map(int,reference_pair)
    d_ref=D[excited,ground]
    if abs(d_ref)<1e-12:raise ValueError('Reference transition has vanishing dipole')
    drive=np.zeros_like(hstatic)
    for a,ia in enumerate(inds):
        for b,ib in enumerate(inds):
            if isexc[a]==1 and isexc[b]==0:
                # +Omega/2 relative to the user-calibrated reference transition.
                drive[a,b]=(omega_reference_khz/2)*D[ia,ib]/d_ref
                drive[b,a]=np.conj(drive[a,b])
    return hstatic,drive


def make_gate_interface(data, states, output=None, microwave_axis='x'):
    """Output compact 0/1/e/e_prime subspace matrices and state IDs.

    Does NOT claim an accurate gate model: spectroscopy must be validated,
    polarization known and actual pulse sequence constructed separately.
    """
    ids={key:states[key]['index'] for key in ('0','1','e','e_prime')}
    inds=[ids[key] for key in ('0','1','e','e_prime')]
    freq=float(data['energies_khz'][ids['e']]-data['energies_khz'][ids['1']])
    H0,Hdrive=single_molecule_rwa_blocks(data,inds,freq,1.0,
                                          (ids['1'],ids['e']),microwave_axis)
    if output is not None:
        np.savez_compressed(output,indices=np.asarray(inds),H_off_khz=H0,
                            H_drive_at_1khz_khz=Hdrive,
                            microwave_carrier_khz=freq)
    return ids,H0,Hdrive
