"""Fixed-antenna Model A: shared microwave polarization for transfer and echo.

Canonical 11-state basis includes aux51. Uses only the frozen spectroscopy export.
"""
from dataclasses import replace
import numpy as np
from scipy.linalg import expm
from nacs.gate.model import FrozenSpectrum, PairModel, BASIS_SETS
from nacs.gate.propagate import Controls, pulse_matrix, dwell_matrix, TRANSFER_RABI_KHZ, ECHO_RABI_KHZ, ECHO_EPS
from nacs.controls.polarization import rotated_linear_polarization,target_rabi,conditional_fidelity,nominal_z_calibration,Uncertainty,error_points
from nacs.data import frozen_spectrum
frozen=frozen_spectrum()
frozen.indices.update({'aux51':51, 'aux79':79, 'extra78':78, 'extra80':80})
SEARCH_BASIS=BASIS_SETS['six']+('aux51','aux79')
VALID_BASIS=BASIS_SETS['nine']+('aux51','aux79')
EXTRA_BASIS=VALID_BASIS+('extra78','extra80')

def echo_rabi(model,angle_deg,amp_error=0.):
    """At original peak E-field of 12.83us echo, but one fixed antenna axis."""
    i,j=model.local['0'],model.local['e']
    d=sum(rotated_linear_polarization(angle_deg)[k]*model.d[k][i,j] for k in range(3))
    dref=sum(float(ECHO_EPS[k])*model.d[k][i,j] for k in range(3))
    return ECHO_RABI_KHZ*abs(d/dref)*(1.+amp_error)


def local_transfer(model, initial,target,frabi,ms,delta,pol,phase=0.):
    """Compute single-molecule transition probability using exact local drive."""
    h,nu,_=model.mw_hamiltonian(initial,target,frabi,delta,pol,phase)
    n=model.n
    local=np.diag(model.energies-nu*model.n_exc).astype(complex)
    for a in range(n):
        for b in range(n):
            if a!=b:local[a,b]=h[a*n,b*n]
    u=expm(-2j*np.pi*local*ms)
    return float(abs(u[model.local[target],model.local[initial]])**2)


def gate_matrix(model,c:Controls,angle_deg,err=Uncertainty()):
    theta=angle_deg+err.angle_deg
    pol=rotated_linear_polarization(theta)
    fr=target_rabi(model,theta,err.amp)
    fe=echo_rabi(model,theta,err.amp)
    d=c.transfer_detuning_khz+err.det_khz
    trans=pulse_matrix(model,'1','e',fr,c.transfer_ms,d,pol,c.transfer_phase,t_start_ms=0.)
    half=dwell_matrix(model,c.dwell_ms/2)
    echo=pulse_matrix(model,'0','e',fe,c.echo_ms,c.echo_detuning_khz,pol,c.echo_phase,
       t_start_ms=c.transfer_ms+c.dwell_ms/2)
    retr=pulse_matrix(model,'1','e',fr,c.transfer_ms,d,pol,c.reverse_phase,
       t_start_ms=c.transfer_ms+c.dwell_ms+c.echo_ms)
    inp=np.eye(model.size,dtype=complex)[:,model.logical_indices]
    out=retr@(half@(echo@(half@(trans@inp))))
    return out[model.logical_indices,:]


def metrics(model,c,theta):
    M=gate_matrix(model,c,theta)
    z=nominal_z_calibration(M)
    f=conditional_fidelity(M,z)
    leakage=1.-float(np.vdot(M,M).real/4.)
    return {'fidelity':f,'leakage':leakage,'virtual_z':z.tolist(),
            'total_duration_ms':float(2*c.transfer_ms+c.dwell_ms+c.echo_ms)}


def eval_ensemble(c,theta,errors,basis=VALID_BASIS):
    # Unique models only once per exchange value
    models={e.exchange:PairModel(frozen,basis,scale_dd=1.+e.exchange) for e in errors}
    nom=models.get(0.,PairModel(frozen,basis))
    M0=gate_matrix(nom,c,theta)
    z=nominal_z_calibration(M0)
    fs=np.array([conditional_fidelity(gate_matrix(models[e.exchange],c,theta,e),z) for e in errors])
    return {'nominal':float(conditional_fidelity(M0,z)),'mean':float(fs.mean()),
            'p10':float(np.percentile(fs,10)),'p25':float(np.percentile(fs,25)),
            'minimum':float(fs.min()),'tail25':float(np.sort(fs)[:max(1,len(fs)//4)].mean()),
            'z':z.tolist(),'values':fs}
