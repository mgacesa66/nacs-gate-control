"""Polarization-axis and correlated microwave dressing in frozen NaCs gate.

Never alter frozen molecule data. Frequencies H/h are kHz, times ms.
The microwave generator holds its **electric-field** amplitude, not target Rabi,
fixed as polarization is rotated. Therefore target Rabi and AC Stark shift
both change, including under fractional microwave amplitude errors.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from nacs.gate.model import PairModel, FrozenSpectrum, BASIS_SETS
from nacs.gate.propagate import Controls, pulse_matrix, dwell_matrix, ECHO_RABI_KHZ, ECHO_EPS, TRANSFER_RABI_KHZ,TRANSFER_EPS
from nacs.gate.metrics import U_EXPERIMENTAL

REFERENCE=np.asarray(TRANSFER_EPS,float)/np.linalg.norm(TRANSFER_EPS)


def rotated_linear_polarization(angle_deg:float):
    """Rigid rotation of real MW electric-field direction around z (quantization).

    This is a *conditional laboratory knob*: requires adjustable antenna
    orientation or a calibrated pair of coherent orthogonal channels.
    Holding z component and norm fixed avoids changing microwave power.
    """
    t=np.deg2rad(angle_deg)
    x,y,z=REFERENCE
    return np.array([x*np.cos(t)-y*np.sin(t), x*np.sin(t)+y*np.cos(t),z])


def reference_dipole(model:PairModel):
    i,j=model.local['1'],model.local['e']
    return abs(sum(REFERENCE[k]*model.d[k][i,j] for k in range(3)))


def target_rabi(model:PairModel,angle_deg:float, fractional_amplitude_error:float=0.0):
    pol=rotated_linear_polarization(angle_deg)
    i,j=model.local['1'],model.local['e']
    d=abs(sum(pol[k]*model.d[k][i,j] for k in range(3)))
    return TRANSFER_RABI_KHZ*(1.+fractional_amplitude_error)*d/reference_dipole(model)


def complete_gate(model:PairModel,controls:Controls,angle_deg:float=0.0,
                  fractional_amplitude_error:float=0., detuning_error_khz:float=0.,
                  angle_error_deg:float=0.0):
    """Full sequence and global output-phase convention from frozen gate.

    MW amplitude noise multiplies ALL transfer couplings, hence changes
    target Rabi and off-resonant Stark shift together. No fake independent
    AC-Stark random parameter. Separate detuning represents residual field or
    carrier drift; if already accounting for power-to-frequency noise via
    amplitude error it must not double-count the same mechanism.
    """
    pol=rotated_linear_polarization(angle_deg+angle_error_deg)
    fr=target_rabi(model,angle_deg+angle_error_deg,fractional_amplitude_error)
    d=controls.transfer_detuning_khz+detuning_error_khz
    ttrans=pulse_matrix(model,'1','e',fr,controls.transfer_ms,d,pol,controls.transfer_phase,
       t_start_ms=0.)
    ud=dwell_matrix(model,controls.dwell_ms/2)
    echo=pulse_matrix(model,'0','e',ECHO_RABI_KHZ,controls.echo_ms,
       controls.echo_detuning_khz,ECHO_EPS,controls.echo_phase,
       t_start_ms=controls.transfer_ms+controls.dwell_ms/2)
    reverse=pulse_matrix(model,'1','e',fr,controls.transfer_ms,d,pol,controls.reverse_phase,
       t_start_ms=controls.transfer_ms+controls.dwell_ms+controls.echo_ms)
    ins=np.eye(model.size,dtype=complex)[:,model.logical_indices]
    out=reverse@(ud@(echo@(ud@(ttrans@ins))))
    return out[model.logical_indices,:]


def single_molecule_transfer(model:PairModel,controls:Controls,angle_deg=0, amp_error=0., det_error=0.):
    pol=rotated_linear_polarization(angle_deg)
    ftarget=target_rabi(model,angle_deg,amp_error)
    H,delta,_=model.mw_hamiltonian('1','e',ftarget,controls.transfer_detuning_khz+det_error,pol)
    n=model.n
    local=np.diag(model.energies-delta*model.n_exc).astype(complex)
    for a in range(n):
        for b in range(n):
            if a!=b:local[a,b]=H[a*n,b*n]
    u=expm(-2j*np.pi*local*controls.transfer_ms)
    return float(abs(u[model.local['e'],model.local['1']])**2)


def second_order_stark(frozen:FrozenSpectrum,angle_deg, mw_amplitude_scale=1.):
    """Perturbative differential Stark shift in complete 128-level N<=1 basis.
    Exclude resonant 1<->e coupling from diagonal terms; this is *not* the
    exact driven resonance and can fail at near-degeneracies/strong driving.
    """
    e=frozen.energies; ind=frozen.indices; a,b=ind['1'],ind['e']
    pol=rotated_linear_polarization(angle_deg)
    ref=sum(REFERENCE[k]*frozen.dipoles[k][a,b] for k in range(3))
    dpol=sum(pol[k]*frozen.dipoles[k] for k in range(3))
    d=np.abs(dpol);H=TRANSFER_RABI_KHZ*mw_amplitude_scale*d/(2*abs(ref))
    exc=np.arange(128)>=32
    nu=e[b]-e[a]
    en=e-e[a]-nu*exc
    sh_a=sum(H[a,j]**2/(en[a]-en[j]) for j in range(32,128) if j!=b)
    sh_b=sum(H[b,j]**2/(en[b]-en[j]) for j in range(32) if j!=a)
    return {'angle_deg':float(angle_deg),'target_rabi_khz':float(TRANSFER_RABI_KHZ*mw_amplitude_scale*d[a,b]/abs(ref)),
      'stark_khz':float(sh_b-sh_a),'d_stark_d_ln_amp_khz':float(2*(sh_b-sh_a)),
      'bright_0e_debye':float(d[ind['0'],b]),'target_1e_debye':float(d[a,b]),
      'ratio_bright_to_target':float(d[ind['0'],b]/max(d[a,b],1e-20))}


def trace_coeff(M):
    t=U_EXPERIMENTAL
    return np.asarray([sum(np.conj(t[j,i])*M[j,i] for i in range(4)) for j in range(4)])


def conditional_fidelity(M,angles):
    a,b=angles
    z=np.exp(0.5j*np.array([-a-b,-a+b,a-b,a+b]))
    s=np.sum(abs(M)**2).real
    return float(np.clip((s+abs(trace_coeff(M)@z)**2)/20.,0.,1.))


def nominal_z_calibration(M):
    from nacs.gate.metrics import metrics
    return np.asarray(metrics(M)['local_z_angles_rad'],float)


@dataclass(frozen=True)
class Uncertainty:
    amp:float=0.
    det_khz:float=0.
    exchange:float=0.
    angle_deg:float=0.


def error_points(seed=9184,n=120, amplitude=.02, detuning=.5, exchange=.02, polarization_deg=0.):
    rng=np.random.default_rng(seed)
    errs=rng.uniform(-1.,1.,(n,4))*[amplitude,detuning,exchange,polarization_deg]
    return [Uncertainty(*x) for x in errs]


class Ensemble:
    def __init__(self,frozen,basis='six',errors=None):
        self.errors=list(errors) if errors is not None else [Uncertainty()]
        self.models={float(j):PairModel(frozen,BASIS_SETS[basis],scale_dd=1.+float(j)) for j in sorted(set(e.exchange for e in self.errors)|{0.})}
        self.nominal=self.models[0.]
    def matrices(self,c,theta):
        return [complete_gate(self.models[e.exchange],c,theta,e.amp,e.det_khz,e.angle_deg) for e in self.errors]
    def evaluate(self,c,theta,calib=None):
        nominal=complete_gate(self.nominal,c,theta)
        z=nominal_z_calibration(nominal) if calib is None else np.asarray(calib)
        f=np.array([conditional_fidelity(m,z) for m in self.matrices(c,theta)])
        return {'nominal':conditional_fidelity(nominal,z),'mean':float(f.mean()),'p10':float(np.percentile(f,10)),
           'minimum':float(f.min()),'tail25':float(np.sort(f)[:max(1,int(.25*len(f)))].mean()),
           'z':z.tolist(),'values':f}
