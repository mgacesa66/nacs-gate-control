"""Piecewise coherent propagation, with explicit carrier-frame bookkeeping."""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm

@dataclass(frozen=True)
class Controls:
    transfer_ms:float=0.307
    dwell_ms:float=0.664
    echo_ms:float=0.01283
    transfer_detuning_khz:float=0.0
    echo_detuning_khz:float=0.0
    transfer_phase:float=0.0
    echo_phase:float=0.0
    reverse_phase:float=0.0

# Mixed polarization is an *inferred illustrative hypothesis* from frozen manifest.
TRANSFER_EPS=(0.80536,0.57825,0.13047)
ECHO_EPS=(1.,0.,0.)
TRANSFER_RABI_KHZ=1/(2*0.307)
ECHO_RABI_KHZ=1/(2*0.01283)


def pulse_matrix(model, initial, target, f_rabi_khz, t_ms,
                 delta_khz=0., epsilon=(1.,0.,0.),phase=0.,
                 dd_during=True, t_start_ms=0., phase_convention="continuous"):
    h,drive_delta,_=model.mw_hamiltonian(initial,target,f_rabi_khz,
                                          delta_khz,epsilon,phase)
    if not dd_during:h=h-model.Hdd
    u_rot=expm(-2j*np.pi*h*t_ms)
    # Convert the instantaneous drive-frame propagator into the common
    # rotational reference frame. Each pulse's carrier phase is specified at
    # the segment start, not inferred from an arbitrary global time origin.
    if phase_convention not in ('continuous','reset'):
        raise ValueError('phase_convention must be continuous or reset')
    t0=t_start_ms if phase_convention=='continuous' else 0.
    frame_post=np.exp(-2j*np.pi*drive_delta*(t0+t_ms)*model.num_exc_diag)
    frame_pre=np.exp(+2j*np.pi*drive_delta*t0*model.num_exc_diag)
    return frame_post[:,None]*u_rot*frame_pre[None,:]


def dwell_matrix(model,t_ms):
    return expm(-2j*np.pi*model.Hfree*t_ms)


def gate_matrix(model,c:Controls, dd_during_pulses=True, echo=True, phase_convention="continuous",
                transfer_rabi_scale=1., echo_rabi_scale=1.):
    ttrans=pulse_matrix(model,'1','e',TRANSFER_RABI_KHZ*transfer_rabi_scale,c.transfer_ms,
        c.transfer_detuning_khz,TRANSFER_EPS,c.transfer_phase,dd_during_pulses,
        t_start_ms=0.,phase_convention=phase_convention)
    if echo:
        ue=pulse_matrix(model,'0','e',ECHO_RABI_KHZ*echo_rabi_scale,c.echo_ms,
            c.echo_detuning_khz,ECHO_EPS,c.echo_phase,dd_during_pulses,
            t_start_ms=c.transfer_ms+c.dwell_ms/2,phase_convention=phase_convention)
        ud=dwell_matrix(model,c.dwell_ms/2)
        interaction=ud@ue@ud
    else:
        interaction=dwell_matrix(model,c.dwell_ms)
    reverse=pulse_matrix(model,'1','e',TRANSFER_RABI_KHZ*transfer_rabi_scale,c.transfer_ms,
        c.transfer_detuning_khz,TRANSFER_EPS,c.reverse_phase,dd_during_pulses,
        t_start_ms=c.transfer_ms+c.dwell_ms+c.echo_ms if echo else c.transfer_ms+c.dwell_ms,
        phase_convention=phase_convention)
    v=np.eye(model.size,dtype=complex)[:,model.logical_indices]
    out=reverse@interaction@ttrans@v
    return out[model.logical_indices,:],out


def transfer_probability(model,t_ms,detuning_khz=0):
    # Full pair Hdd=0 is not needed; use an isolated molecule, but all
    # selected-level spectator transitions are retained.
    # Single-particle transfer dynamics are identical to the associated
    # two-particle sub-block with Hdd=0.
    h,delta,_=model.mw_hamiltonian('1','e',TRANSFER_RABI_KHZ,
             detuning_khz, TRANSFER_EPS)
    # Extract the full local Hamiltonian including its other off-resonant transitions.
    local=np.diag(model.energies-delta*model.n_exc).astype(complex)
    # Use diagonal of the one-body tensor-sum to reconstruct full drive matrix.
    # h = local⊗I + I⊗local + Hdd; offdiagonal ground/rot coupling is h[(a,0),(b,0)]
    for a in range(model.n):
      for b in range(model.n):
        if a!=b:
          local[a,b]=h[a*model.n, b*model.n]
    u=expm(-2j*np.pi*local*t_ms)
    return float(abs(u[model.local['e'],model.local['1']])**2)
