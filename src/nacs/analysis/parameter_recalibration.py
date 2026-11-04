"""Publication Step 1: molecular-parameter perturbation and microwave recalibration.

Frequency units cyclic kHz, time ms, dipoles Debye. Uses full microscopic
N=0,1 (128-dimensional) re-diagonalization and corrected 11-state pair model.

Three physical operational scenarios:
  fixed: both absolute carriers and microwave electric-field peaks fixed;
  frequency: carrier frequencies follow the tracked bare dressed transitions;
  frequency_rabi: additionally transfer/echo E-fields are separately adjusted
                  to recover the BASELINE target Rabi magnitudes.

All scenarios use identical stored pulse durations, phases and gate controls.
No gate-level reoptimization, no experimental noise data, no fitted parameters.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
from functools import lru_cache
import numpy as np
from scipy.linalg import expm
from scipy.optimize import linear_sum_assignment
from nacs.data import frozen_spectrum,benchmark_path
from nacs.spectroscopy.constants import NaCsParameters
from nacs.spectroscopy.core import NaCsMolecule
from nacs.gate.model import FrozenSpectrum,PairModel
from nacs.gate.propagate import ECHO_RABI_KHZ,TRANSFER_RABI_KHZ,ECHO_EPS,TRANSFER_EPS
from nacs.controls.fixed_polarization import VALID_BASIS,rotated_linear_polarization,echo_rabi
from nacs.controls.spectator import (ANGLE,BASE,Design,_pair_rot_frame,measures,
                                      echo_ref_matrix,gate_output)
from nacs.controls.gaussian_echo import echo_matrix
from nacs.controls.polarization import Uncertainty,conditional_fidelity,nominal_z_calibration

@dataclass(frozen=True)
class Probe:
    name: str
    field_G: float=0.
    tensor_frac: float=0.
    ellipticity_deg: float=0.

PROBES=(Probe('nominal'),
        Probe('B_minus_0p1G',field_G=-.1),Probe('B_plus_0p1G',field_G=.1),
        Probe('tensor_minus_0p5pct',tensor_frac=-.005),Probe('tensor_plus_0p5pct',tensor_frac=.005),
        Probe('ellipse_minus_0p1deg',ellipticity_deg=-.1),Probe('ellipse_plus_0p1deg',ellipticity_deg=.1),
        Probe('ellipse_minus_0p25deg',ellipticity_deg=-.25),Probe('ellipse_plus_0p25deg',ellipticity_deg=.25),
        Probe('ellipse_minus_0p5deg',ellipticity_deg=-.5),Probe('ellipse_plus_0p5deg',ellipticity_deg=.5))
MODES=('fixed','frequency','frequency_rabi')

@lru_cache(maxsize=1)
def nominal_objects():
    import json
    from nacs.data import benchmark_path
    f=frozen_spectrum()
    f.indices.update({'aux51':51,'aux79':79,'extra78':78,'extra80':80})
    from nacs.data import ROOT
    # Access extracted original baseline through standard package data path.
    root=ROOT
    archive=np.load(root/'benchmarks/frozen/dressed_128.npz')
    manifest=json.loads((root/'benchmarks/frozen/baseline_manifest.json').read_text())
    return f,archive['vectors'].copy(),manifest


def _params(probe):
    f,V,manifest=nominal_objects()
    chi=np.deg2rad(manifest['optical_ellipticity_deg']+probe.ellipticity_deg)
    return NaCsParameters(field_G=manifest['magnetic_field_G']+probe.field_G,
         optical_tensor_khz=manifest['optical_tensor_khz']*(1+probe.tensor_frac),
         polarization_xyz=(complex(np.cos(chi)),complex(0,np.sin(chi)),0j),dipole_debye=4.6)


@lru_cache(maxsize=24)
def tracked_spectrum(probe):
    """Track labels by global max overlap and align eigenvector gauge to nominal.

    Returns diagnostics including assignment ambiguity, which is not resolved by
    a small Hilbert basis convergence test.
    """
    fn,V0,manifest=nominal_objects()
    p=_params(probe)
    molecule=NaCsMolecule(p)
    s=molecule.diagonalize()
    overlaps=V0.conj().T@s.vectors
    ia,ib=linear_sum_assignment(-np.abs(overlaps)**2)
    assign=np.empty(128,dtype=int);assign[ia]=ib
    fidelity=np.abs(overlaps[np.arange(128),assign])**2
    V=s.vectors.copy()
    for original,assigned in enumerate(assign):
        ip=overlaps[original,assigned]
        if abs(ip)>0: V[:,assigned]*=np.conjugate(ip)/abs(ip)
    bare=molecule.bare_dipole_q_debye()
    d={q:V.conj().T@bare[q]@V for q in (-1,0,1)}
    cart=( (d[-1]-d[1])/np.sqrt(2),1j*(d[-1]+d[1])/np.sqrt(2),d[0] )
    # Preserve fixed state names, permit newly selected states to switch numeric index.
    labels={name:int(assign[i]) for name,i in fn.indices.items()}
    labels.update({f'rank{i}':int(assign[i]) for i in (38,69,40,72)})
    new=FrozenSpectrum(s.energies_khz,cart,labels,'perturbed-not-frozen')
    tracked={name:float(fidelity[i]) for name,i in fn.indices.items()}
    # The prior audit also examined four coupling-ranked spectators beyond 11 states.
    ranked_indices=(38,69,40,72)
    tracked_ranked={str(i):float(fidelity[i]) for i in ranked_indices}
    return new,{'tracking':tracked,'tracking_ranked':tracked_ranked,
                'min_tracking_ranked15':min([tracked[name] for name in VALID_BASIS]+list(tracked_ranked.values())),
                'min_tracking':min(tracked.values()),
                'min_tracking_primary':min(tracked[name] for name in VALID_BASIS),
                'mapping':labels,'parameters':p}


def _dipole(f,a,b,pol):
    return abs(sum(pol[k]*f.dipoles[k][f.indices[a],f.indices[b]] for k in range(3)))


def drive_parameters(f,mode,design):
    """Physical fixed-E-field calibrations relative to immutable base state.

    Retuning by *bare* tracked resonance centers does not correct multi-level
    AC Stark shifts. It is an explicitly optimistic carrier-tracking proxy.
    """
    old,_,_=nominal_objects()
    pol=rotated_linear_polarization(ANGLE)
    tr_orig=_dipole(old,'1','e',np.asarray(TRANSFER_EPS,dtype=float)/np.linalg.norm(TRANSFER_EPS))
    ec_orig=_dipole(old,'0','e',np.asarray(ECHO_EPS,dtype=float))
    ft=TRANSFER_RABI_KHZ*_dipole(f,'1','e',pol)/tr_orig
    fe=ECHO_RABI_KHZ*_dipole(f,'0','e',pol)/ec_orig
    ft0=TRANSFER_RABI_KHZ*_dipole(old,'1','e',pol)/tr_orig
    fe0=ECHO_RABI_KHZ*_dipole(old,'0','e',pol)/ec_orig
    shift_tr=(f.energies[f.indices['e']]-f.energies[f.indices['1']])-(old.energies[old.indices['e']]-old.energies[old.indices['1']])
    shift_ec=(f.energies[f.indices['e']]-f.energies[f.indices['0']])-(old.energies[old.indices['e']]-old.energies[old.indices['0']])
    if mode=='fixed':
        det_tr=design.transfer_detuning_khz-shift_tr
        det_ec=BASE.echo_detuning_khz-shift_ec
        efield_t=efield_e=1.
    elif mode=='frequency':
        det_tr=design.transfer_detuning_khz
        det_ec=BASE.echo_detuning_khz
        efield_t=efield_e=1.
    elif mode=='frequency_rabi':
        det_tr=design.transfer_detuning_khz
        det_ec=BASE.echo_detuning_khz
        efield_t=ft0/ft
        efield_e=fe0/fe
        ft=ft0;fe=fe0
    else: raise ValueError(mode)
    return {'transfer_rabi_khz':ft,'echo_rabi_khz':fe,
            'transfer_detuning_khz':det_tr,'echo_detuning_khz':det_ec,
            'transfer_line_shift_khz':shift_tr,'echo_line_shift_khz':shift_ec,
            'electric_field_scale_transfer':efield_t,'electric_field_scale_echo':efield_e,
            'transfer_dipole_ratio':_dipole(f,'1','e',pol)/_dipole(old,'1','e',pol),
            'echo_dipole_ratio':_dipole(f,'0','e',pol)/_dipole(old,'0','e',pol)}


def full_gate_output(model,design,params,echo,Uphase=False):
    """Coherent physical-field gate; returned columns represent all logical inputs.

    Gaussian echo propagator is cached at the reference start and conjugated
    correctly to the actual start, including the recalibrated carrier delta.
    """
    pol=rotated_linear_polarization(ANGLE)
    inp=np.eye(model.size,dtype=complex)[:,model.logical_indices]
    ft=params['transfer_rabi_khz'];det=params['transfer_detuning_khz']
    transfer=_pair_rot_frame(model,'1','e',ft,det,pol,0.,design.phase_forward,
                 design.transfer_ms,design.fraction_first,0.)
    half=expm(-2j*np.pi*model.Hfree*(design.dwell_ms/2))
    effective_echo_det=params['echo_detuning_khz']
    delta_echo=model.transition_khz('0','e')+effective_echo_det-model.nu_ref
    shift=design.transfer_ms+design.dwell_ms/2-(BASE.transfer_ms+BASE.dwell_ms/2)
    phase=np.exp(-2j*np.pi*delta_echo*shift*model.num_exc_diag)
    shifted_echo=(phase[:,None]*echo)*phase.conj()[None,:]
    reverse_start=design.transfer_ms+design.dwell_ms+BASE.echo_ms
    reverse=_pair_rot_frame(model,'1','e',ft,det,pol,
                 design.reverse_phase+design.phase_reverse,design.reverse_phase,
                 design.transfer_ms,design.fraction_first,reverse_start)
    return reverse@(half@(shifted_echo@(half@(transfer@inp))))


def build_echo(model,params,steps=64):
    baseline=replace(BASE,echo_detuning_khz=params['echo_detuning_khz'])
    # The Gaussian echo wrapper normalizes the target Rabi by NEW dipoles,
    # whereas physical fixed E-field must keep its original normalization.
    ratio=params['echo_rabi_khz']/echo_rabi(model,ANGLE)
    return echo_matrix(model,ANGLE,baseline,err=Uncertainty(amp=ratio-1.),steps=steps)


def evaluate(probe,mode,designs,steps=64,reference_z=None):
    f,info=tracked_spectrum(probe)
    m=PairModel(f,VALID_BASIS)
    # Fields may change with protocol only if transfer amplitude calibration depends
    # on the per-protocol original carrier; here target field is shared across designs.
    p0=drive_parameters(f,mode,next(iter(designs.values())))
    ec=build_echo(m,p0,steps=steps)
    results=[]
    for name,d in designs.items():
        p=drive_parameters(f,mode,d)
        out=full_gate_output(m,d,p,ec)
        q=measures(m,out)
        M=q['M']
        zf=reference_z[name] if reference_z is not None else q['z']
        results.append({'probe':probe.name,'mode':mode,'protocol':name,
            'fidelity_phase_refit':float(q['fidelity']),
            'fidelity_fixed_nominal_z':float(conditional_fidelity(M,zf)),
            'leakage':float(q['leakage']),'aux51':float(q['aux51']),
            'other_spectators':float(q['other_spectators']),
            'total_spectators':float(q['total_spectators']),
            'intended_nonlogical':float(q['intended_nonlogical']),
            'z1_rad':float(q['z'][0]),'z2_rad':float(q['z'][1]),
            **p})
    return results,info
