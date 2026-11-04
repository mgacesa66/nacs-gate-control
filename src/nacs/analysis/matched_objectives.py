"""Step 3: preregistered, resource-matched control-objective comparison.

Frozen spectroscopy, fixed 35-degree antenna, common Gaussian echo, identical
peak microwave field, carrier/duration/phase bounds and finite phase ramps.
A 'selected spectator' penalty is not a complete-channel cost; check it later.
"""
from __future__ import annotations
from dataclasses import asdict
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from nacs.controls.spectator import (Design,ref_design,design_from_vector,
    vector_from_design,BOUNDS,gate_output,measures,echo_ref_matrix,ANGLE,BASE,
    PairModel,frozen,SEARCH_BASIS,VALID_BASIS,rotated_linear_polarization,target_rabi)
from nacs.controls.gaussian_echo import envelope
from nacs.controls.polarization import Uncertainty,conditional_fidelity,nominal_z_calibration

# Frozen before numerical search, so there is no post-hoc retuning of objectives.
OBJECTIVES=('isolated_transfer','selected_aux51','channel_fidelity','channel_guarded')
LAMBDA=2.0
PHASE_RAMP_US=2.0
RAMP_STEPS=6
ECHO_STEPS_TRAIN=16
ECHO_STEPS_VALID=64
MAX_GATE_MS=1.30


def single_local_hamiltonian(model,det,phase,amp=0.):
    """One-molecule Hamiltonian in the same RWA/carrier frame as PairModel."""
    pol=rotated_linear_polarization(ANGLE)
    fr=target_rabi(model,ANGLE,amp)
    dpol=sum(pol[k]*model.d[k] for k in range(3))
    a,b=model.local['1'],model.local['e']
    r=abs(dpol[a,b]); nu=model.transition_khz('1','e')+det
    delta=nu-model.nu_ref
    H=np.diag(model.energies-delta*model.n_exc).astype(complex)
    for j in range(model.n):
        for k in range(model.n):
            if model.n_exc[j]==0 and model.n_exc[k]==1:
                H[j,k]=fr*dpol[j,k]*np.exp(1j*phase)/(2*r)
                H[k,j]=H[j,k].conjugate()
    return H


def single_phase_pulse(model, duration,det,phase0,phase1,amp=0.,ramp_us=PHASE_RAMP_US,
                       ramp_steps=RAMP_STEPS):
    """Single molecule pulse, for isolated |1>->|e> and |e>->|1> populations.

    Carrier-frame boundary factors are diagonal phases and do not change the
    isolated transfer population, so they are unnecessary for this observable.
    """
    h0=single_local_hamiltonian(model,det,phase0,amp)
    h1=single_local_hamiltonian(model,det,phase1,amp)
    ramp=1e-3*ramp_us
    ta=duration/2-ramp/2
    if ta<0:raise ValueError('Ramp exceeds pulse duration')
    U=expm(-2j*np.pi*h0*ta)
    if abs(phase1-phase0)>1e-14 and ramp>0:
        for k in range(ramp_steps):
            s=(k+.5)/ramp_steps
            phase=phase0+(phase1-phase0)*(.5-.5*np.cos(np.pi*s))
            U=expm(-2j*np.pi*single_local_hamiltonian(model,det,phase,amp)*ramp/ramp_steps)@U
    else:
        U=expm(-2j*np.pi*h0*ramp)@U
    return expm(-2j*np.pi*h1*ta)@U


def isolated_transfers(model,d,amp=0.,det_error=0.):
    F=single_phase_pulse(model,d.transfer_ms,d.transfer_detuning_khz+det_error,
                         0.,d.phase_forward,amp)
    R=single_phase_pulse(model,d.transfer_ms,d.transfer_detuning_khz+det_error,
                         d.reverse_phase+d.phase_reverse,d.reverse_phase,amp)
    return (float(abs(F[model.local['e'],model.local['1']])**2),
            float(abs(R[model.local['1'],model.local['e']])**2))


def score(kind,met,pforward,preverse):
    if kind=='isolated_transfer':return 1.-.5*(pforward+preverse)
    if kind=='selected_aux51':return 1.-met['fidelity']+LAMBDA*met['aux51']
    if kind=='channel_fidelity':return 1.-met['fidelity']
    if kind=='channel_guarded':return 1.-met['fidelity']+LAMBDA*met['total_spectators']
    raise ValueError(kind)


def feasible(x):
    d=design_from_vector(x)
    return all(lo-1e-12<=v<=hi+1e-12 for v,(lo,hi) in zip(x,BOUNDS)) and d.duration_ms<=MAX_GATE_MS


def measure_candidate(model,x,echo=None,steps=ECHO_STEPS_VALID,amp=0.,det_error=0.,
                      exchange=0.):
    d=design_from_vector(x)
    if exchange:
        # Caller should pass model constructed for exchange perturbation.
        raise ValueError('Pass the exchange-scaled model instead')
    out=gate_output(model,d,amp_error=amp,det_error_khz=det_error,
           echo_ref=echo,steps=steps,phase_ramp_us=PHASE_RAMP_US,ramp_steps=RAMP_STEPS)
    met=measures(model,out)
    pf,pr=isolated_transfers(model,d,amp,det_error)
    return dict(fidelity=float(met['fidelity']),leakage=float(met['leakage']),
        aux51=float(met['aux51']),other_spectators=float(met['other_spectators']),
        total_spectators=float(met['total_spectators']),
        intended_nonlogical=float(met['intended_nonlogical']),
        transfer_forward=pf,transfer_reverse=pr,
        duration_ms=d.duration_ms,phase_ramp_us=PHASE_RAMP_US,
        virtual_z_rad=np.asarray(met['z']).tolist())


def archived_start_vectors():
    from nacs.data import benchmark_path
    data=json.loads(benchmark_path('optimization/candidates.json').read_text())
    names=('unmodulated','gate_only','aux51','all_spectators')
    arr=[]
    for name in names:
        c=data[name]['design']
        arr.append(vector_from_design(Design(**c)))
    return names,np.array(arr)


def population(seed, size=30):
    """Every arm/seed receives precisely the same initial population."""
    rng=np.random.default_rng(seed)
    z=rng.uniform(0,1,(size,6))
    bounds=np.array(BOUNDS,float)
    x=bounds[:,0]+z*(bounds[:,1]-bounds[:,0])
    _,history=archived_start_vectors()
    # Same historical warm starts for all arms; no preference or hidden extra work.
    for k,v in enumerate(history):
        x[k]=np.clip(v,bounds[:,0],bounds[:,1])
    return x


def optimize_arm(kind,seed=2026,iterations=8,pop_size=30,steps=ECHO_STEPS_TRAIN):
    from scipy.optimize import differential_evolution
    if kind not in OBJECTIVES:raise ValueError(kind)
    model=PairModel(frozen,SEARCH_BASIS)
    echo=echo_ref_matrix(model,steps=steps)
    init=population(seed,pop_size)
    neval=0
    best=[np.inf,None]
    def fn(x):
        nonlocal neval
        neval+=1
        if not feasible(x):
            return 5.+10.*max(0.,design_from_vector(x).duration_ms-MAX_GATE_MS)
        d=design_from_vector(x)
        if kind=='isolated_transfer':
            pf,pr=isolated_transfers(model,d)
            v=score(kind,None,pf,pr)
        else:
            out=gate_output(model,d,echo_ref=echo,steps=steps,
                     phase_ramp_us=PHASE_RAMP_US,ramp_steps=RAMP_STEPS)
            met=measures(model,out)
            v=score(kind,met,0.,0.)
        if v<best[0]:best[:]=[float(v),x.tolist()]
        return float(v)
    opt=differential_evolution(fn,BOUNDS,init=init,maxiter=iterations,seed=seed,
          polish=False,updating='immediate',tol=0,atol=0)
    return dict(objective=kind,seed=seed,iterations=iterations,
       population=pop_size,evaluations=neval,best_objective=float(best[0]),
       vector=best[1],search_basis_n=len(SEARCH_BASIS),echo_steps=steps,
       ramp_us=PHASE_RAMP_US,lambda_penalty=LAMBDA,success=bool(opt.success),
       # scipy's success=False for deliberately capped equal evaluation budget.
       stopped_at_equal_budget=True)
