"""Wavefunction-column propagation, avoiding dense full-pair propagator products.

Hermitian eigendecomposition exactly exponentiates each constant Hamiltonian
interval; midpoint approximation is the same as frozen Gaussian echo model.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import eigh
from nacs.controls.spectator import ANGLE,BASE,rotated_linear_polarization,target_rabi,echo_rabi
from nacs.controls.gaussian_echo import envelope


def apply(h,t,v):
    w,q=eigh(h,check_finite=False,overwrite_a=True,driver='evr')
    return q@(np.exp(-2j*np.pi*t*w)[:,None]*(q.conj().T@v))


def phase_pulse(model,v,initial,target,field,det,pol,phase_a,phase_b,t,t0,ramp_us=2.,ramp_steps=6):
    h_a,delta,_=model.mw_hamiltonian(initial,target,field,det,pol,phase_a)
    h_b,db,_=model.mw_hamiltonian(initial,target,field,det,pol,phase_b)
    assert abs(db-delta)<1e-10
    ramp=ramp_us*.001
    ta=t/2-ramp/2
    if ta<0:raise ValueError('Negative constant segment')
    num=model.num_exc_diag
    u=np.exp(+2j*np.pi*delta*t0*num)[:,None]*v
    u=apply(h_a,ta,u)
    for k in range(ramp_steps):
        s=(k+.5)/ramp_steps
        phase=phase_a+(phase_b-phase_a)*(.5-.5*np.cos(np.pi*s))
        h,_,_=model.mw_hamiltonian(initial,target,field,det,pol,phase)
        u=apply(h,ramp/ramp_steps,u)
    u=apply(h_b,ta,u)
    return np.exp(-2j*np.pi*delta*(t0+t)*num)[:,None]*u


def gaussian_echo(model,v,t0,steps=24,amp_error=0.):
    pol=rotated_linear_polarization(ANGLE)
    field=echo_rabi(model,ANGLE,amp_error)
    h,delta,_=model.mw_hamiltonian('0','e',field,BASE.echo_detuning_khz,pol,BASE.echo_phase)
    base,_,_=model.mw_hamiltonian('0','e',0.,BASE.echo_detuning_khz,pol,BASE.echo_phase)
    drive=h-base
    t=BASE.echo_ms;dt=t/steps
    num=model.num_exc_diag
    u=np.exp(+2j*np.pi*delta*t0*num)[:,None]*v
    for k in range(steps):
        a=float(envelope((k+.5)/steps,'truncated'))
        u=apply(base+a*drive,dt,u)
    return np.exp(-2j*np.pi*delta*(t0+t)*num)[:,None]*u


def gate_columns(model,d,steps=24,ramp_us=2.,ramp_steps=6,amp_error=0.,det_error=0.):
    pol=rotated_linear_polarization(ANGLE)
    field=target_rabi(model,ANGLE,amp_error)
    init=np.eye(model.size,dtype=complex)[:,model.logical_indices]
    x=phase_pulse(model,init,'1','e',field,d.transfer_detuning_khz+det_error,
         pol,0.,d.phase_forward,d.transfer_ms,0.,ramp_us,ramp_steps)
    x=apply(model.Hfree,d.dwell_ms/2,x)
    x=gaussian_echo(model,x,d.transfer_ms+d.dwell_ms/2,steps,amp_error)
    x=apply(model.Hfree,d.dwell_ms/2,x)
    t0=d.transfer_ms+d.dwell_ms+BASE.echo_ms
    x=phase_pulse(model,x,'1','e',field,d.transfer_detuning_khz+det_error,
         pol,d.reverse_phase+d.phase_reverse,d.reverse_phase,d.transfer_ms,t0,ramp_us,ramp_steps)
    return x
