"""Spectator-aware, piecewise-phase transfer control with guarded leakage metrics."""
from __future__ import annotations
from dataclasses import dataclass,replace
from functools import lru_cache
import numpy as np
from scipy.linalg import expm
from nacs.analysis.stage_leakage import (PairModel,frozen,VALID_BASIS,archived_controls,
 rotated_linear_polarization,target_rabi,echo_rabi,nominal_z_calibration,
 conditional_fidelity,initial_states)
from nacs.controls.gaussian_echo import echo_matrix
from nacs.gate.propagate import Controls
ANGLE=35.
BASE=archived_controls('gaussian')
SEARCH_BASIS=('0','1','e','e_prime','f_z','f_y','aux51','aux79')
EXTRA_BASIS=VALID_BASIS+('extra78','extra80')

@dataclass(frozen=True)
class Design:
    transfer_ms:float=BASE.transfer_ms
    dwell_ms:float=BASE.dwell_ms
    transfer_detuning_khz:float=BASE.transfer_detuning_khz
    phase_forward:float=0.
    phase_reverse:float=0.
    reverse_phase:float=BASE.reverse_phase
    fraction_first:float=.5

    @property
    def duration_ms(self):return 2*self.transfer_ms+self.dwell_ms+BASE.echo_ms


def ref_design():return Design()


def design_from_vector(x):
    return Design(transfer_ms=float(x[0]),dwell_ms=float(x[1]),
       transfer_detuning_khz=float(x[2]),phase_forward=float(x[3]),
       phase_reverse=float(x[4]),reverse_phase=float(x[5]))

def vector_from_design(d):
    return np.array([d.transfer_ms,d.dwell_ms,d.transfer_detuning_khz,
       d.phase_forward,d.phase_reverse,d.reverse_phase])

BOUNDS=[(.206,.288),(.45,.72),(17.8,22.2),(-np.pi,np.pi),(-np.pi,np.pi),(-.5,.5)]

def _pair_rot_frame(model, initial, target, f, det, pol, phase1,phase2,
                    duration, fraction, tstart, phase_ramp_us=0.,ramp_steps=20):
    """Exact two-piece phase step at fixed microwave carrier frequency.

    Each phase is an explicit DDS phase command, within a common carrier frame.
    No unphysical inter-segment carrier reset is assumed. The same Ppre/post
    transform applies to both segments because drive frequency is unchanged.
    """
    h1,delta,_=model.mw_hamiltonian(initial,target,f,det,pol,phase1)
    if abs(phase2-phase1)<1e-14:
        urot=expm(-2j*np.pi*h1*duration)
    else:
        h2,delta2,_=model.mw_hamiltonian(initial,target,f,det,pol,phase2)
        assert abs(delta-delta2)<1e-10
        ramp=phase_ramp_us*1e-3
        if ramp <= 0.:
            urot=expm(-2j*np.pi*h2*duration*(1-fraction)) @ expm(-2j*np.pi*h1*duration*fraction)
        else:
            # IQ phase interpolated in time with cosine endpoints, at fixed
            # amplitude and unchanged carrier frequency. 20 midpoint slices.
            ta=duration*fraction-ramp/2.
            tb=duration*(1-fraction)-ramp/2.
            if ta<0 or tb<0: raise ValueError('Phase ramp longer than pulse segments')
            urot=expm(-2j*np.pi*h1*ta)
            for k in range(ramp_steps):
                s=(k+.5)/ramp_steps
                phi=phase1+(phase2-phase1)*(1-np.cos(np.pi*s))/2
                h,_,_=model.mw_hamiltonian(initial,target,f,det,pol,phi)
                urot=expm(-2j*np.pi*h*ramp/ramp_steps)@urot
            urot=expm(-2j*np.pi*h2*tb)@urot
    pre=np.exp(+2j*np.pi*delta*tstart*model.num_exc_diag)
    post=np.exp(-2j*np.pi*delta*(tstart+duration)*model.num_exc_diag)
    return post[:,None]*urot*pre[None,:]


def echo_ref_matrix(model, amp_error=0.,steps=64):
    # Frozen Gaussian shape, original echo drive parameters, at its reference t0.
    from nacs.controls.polarization import Uncertainty
    return echo_matrix(model,ANGLE,BASE,err=Uncertainty(amp=amp_error),steps=steps,kind='truncated')


def gate_output(model,d,amp_error=0.,det_error_khz=0., echo_ref=None, steps=64,
                trace=False,phase_ramp_us=0.,ramp_steps=20):
    """Returns full pair amplitudes at selected checkpoints for four logical inputs.

    Optimizer controls only forward/reverse phase step, transfer duration,
    common microwave carrier and dwell, reverse phase baseline. Echo is frozen.
    """
    pol=rotated_linear_polarization(ANGLE)
    fr=target_rabi(model,ANGLE,amp_error)
    inp=initial_states(model)
    forward=_pair_rot_frame(model,'1','e',fr,d.transfer_detuning_khz+det_error_khz,
               pol,0.,d.phase_forward,d.transfer_ms,d.fraction_first,0.,phase_ramp_us,ramp_steps)
    s1=forward@inp
    dwell=expm(-2j*np.pi*model.Hfree*(d.dwell_ms/2))
    s2=dwell@s1
    if echo_ref is None:echo_ref=echo_ref_matrix(model,amp_error,steps=steps)
    # carrier-frame phase update for changed echo start time
    de=model.transition_khz('0','e')+BASE.echo_detuning_khz-model.nu_ref
    shift=(d.transfer_ms+d.dwell_ms/2)-(BASE.transfer_ms+BASE.dwell_ms/2)
    u=np.exp(-2j*np.pi*de*shift*model.num_exc_diag)
    echo=(u[:,None]*echo_ref)*u.conj()[None,:]
    s3=echo@s2
    s4=dwell@s3
    t0=d.transfer_ms+d.dwell_ms+BASE.echo_ms
    reverse=_pair_rot_frame(model,'1','e',fr,d.transfer_detuning_khz+det_error_khz,
       pol,d.reverse_phase+d.phase_reverse,d.reverse_phase,d.transfer_ms,d.fraction_first,t0,phase_ramp_us,ramp_steps)
    s5=reverse@s4
    return (s5,[inp,s1,s2,s3,s4,s5]) if trace else s5


def measures(model,out):
    n=model.n
    logical=np.array(model.logical_indices)
    core={model.pair(a,b) for a in ('0','1','e') for b in ('0','1','e')}
    spec=np.array(sorted(set(range(model.size))-core),int)
    aux=model.local.get('aux51',None)
    aux_rows=np.array([a*n+b for a in range(n) for b in range(n)
                         if a==aux or b==aux],int) if aux is not None else np.array([],int)
    probabilities=abs(out)**2
    paux=float(probabilities[aux_rows].sum()/4.) if len(aux_rows) else 0.
    pspec=float(probabilities[spec].sum()/4.)
    pother=float(probabilities[np.setdiff1d(spec,aux_rows)].sum()/4.)
    # Complement of aux rows is not identical to other spectator states in general;
    # explicitly use spectator minus aux rows for mutually exclusive partition.
    pleak=1.-float(probabilities[logical].sum()/4.)
    M=out[logical,:]
    z=nominal_z_calibration(M)
    F=float(conditional_fidelity(M,z))
    return dict(fidelity=F,leakage=pleak,aux51=paux,other_spectators=pother,
           total_spectators=pspec,intended_nonlogical=pleak-pspec,
           z=z, M=M)


def state_breakdown(model,out):
    n=model.n; P=abs(out)**2
    labels={}
    for key,i in model.local.items():
        if key in ('0','1','e'):continue
        rows=[a*n+b for a in range(n) for b in range(n) if a==i or b==i]
        labels[key]=float(P[rows].sum()/4)
    return labels
