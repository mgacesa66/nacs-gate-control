"""Coherent Gaussian-echo evolution with fixed microwave orientation."""
from dataclasses import replace
from functools import lru_cache
import numpy as np
from scipy.linalg import expm
from nacs.controls.fixed_polarization import (frozen,PairModel,VALID_BASIS,SEARCH_BASIS,echo_rabi,
 rotated_linear_polarization,gate_matrix as square_gate,nominal_z_calibration,conditional_fidelity,
 pulse_matrix,dwell_matrix,target_rabi,Uncertainty)
from nacs.gate.propagate import Controls
def envelope(s,kind='truncated',sigma_fraction=.3):
    """s in [0,1]; peak-normalized; two explicit endpoint conventions.

    'truncated': exp(-(s-.5)^2/(2sigma^2)), abruptly truncated at endpoints;
    'zero_end': subtract value at endpoints, normalize to a unit peak.
    """
    s=np.asarray(s)
    g=np.exp(-.5*((s-.5)/sigma_fraction)**2)
    if kind=='square':
        return np.ones_like(s)
    if kind=='zero_end':
        edge=np.exp(-.5*(.5/sigma_fraction)**2)
        g=(g-edge)/(1-edge)
    elif kind!='truncated':raise ValueError(kind)
    return g

def shape_area(kind='truncated',sigma_fraction=.3):
    from scipy.integrate import trapezoid
    return float(trapezoid(envelope(np.linspace(0,1,10001),kind,sigma_fraction),dx=1/10000))

def echo_matrix(model, angle, c, err=Uncertainty(),steps=64,kind='truncated',sigma_fraction=.3):
    """Return the pair echo propagator with exact free/dipolar dynamics at every step.
    All other Model A segments use the archived propagation unchanged.
    """
    theta=angle+err.angle_deg
    pol=rotated_linear_polarization(theta)
    fr=echo_rabi(model,theta,err.amp)
    h, drive_delta,_=model.mw_hamiltonian('0','e',fr,c.echo_detuning_khz,pol,c.echo_phase)
    h_base,_,_=model.mw_hamiltonian('0','e',0.,c.echo_detuning_khz,pol,c.echo_phase)
    h_drive=h-h_base
    dt=c.echo_ms/steps
    U=np.eye(model.size,dtype=complex)
    for k in range(steps):
        a=float(envelope((k+.5)/steps,kind,sigma_fraction))
        U=expm(-2j*np.pi*(h_base+a*h_drive)*dt)@U
    t0=c.transfer_ms+c.dwell_ms/2
    frame_post=np.exp(-2j*np.pi*drive_delta*(t0+c.echo_ms)*model.num_exc_diag)
    frame_pre=np.exp(+2j*np.pi*drive_delta*t0*model.num_exc_diag)
    return frame_post[:,None]*U*frame_pre[None,:]

def gate_matrix(model,c,angle,err=Uncertainty(),steps=64,kind='truncated'):
    theta=angle+err.angle_deg
    pol=rotated_linear_polarization(theta)
    fr=target_rabi(model,theta,err.amp)
    trans=pulse_matrix(model,'1','e',fr,c.transfer_ms,c.transfer_detuning_khz+err.det_khz,
                       pol,c.transfer_phase,t_start_ms=0.)
    dwell=dwell_matrix(model,c.dwell_ms/2)
    ue=echo_matrix(model,angle,c,err,steps,kind)
    retr=pulse_matrix(model,'1','e',fr,c.transfer_ms,c.transfer_detuning_khz+err.det_khz,
                      pol,c.reverse_phase,t_start_ms=c.transfer_ms+c.dwell_ms+c.echo_ms)
    logical=np.eye(model.size,dtype=complex)[:,model.logical_indices]
    out=retr@(dwell@(ue@(dwell@(trans@logical))))
    return out[model.logical_indices,:]

def isolated_echo_probability(model,angle,T,detuning,steps=64,kind='truncated'):
    theta=angle; pol=rotated_linear_polarization(theta); fr=echo_rabi(model,theta)
    h,nu,_=model.mw_hamiltonian('0','e',fr,detuning,pol,0.)
    hzero,_,_=model.mw_hamiltonian('0','e',0.,detuning,pol,0.)
    # Extract 1-molecule Hamiltonians by indexing pair (a,0),(b,0): Hdd removed.
    # H_free on single molecule uses model energies and carrier detuning.
    n=model.n
    hl=np.diag(model.energies-nu*model.n_exc).astype(complex)
    dh=np.zeros((n,n),complex)
    for a in range(n):
        for b in range(n):
            if a!=b:dh[a,b]=h[a*n,b*n]-hzero[a*n,b*n]
    v=np.eye(n,dtype=complex)[:,model.local['0']]
    dt=T/steps
    for k in range(steps):
        s=float(envelope((k+.5)/steps,kind))
        v=expm(-2j*np.pi*(hl+s*dh)*dt)@v
    return float(abs(v[model.local['e']])**2)

def gate_fidelity(model,c,angle,steps=64,kind='truncated'):
    M=gate_matrix(model,c,angle,steps=steps,kind=kind)
    z=nominal_z_calibration(M)
    f=float(conditional_fidelity(M,z))
    leak=1.-float(np.vdot(M,M).real/4.)
    return f,leak,z
