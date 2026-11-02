"""Stage-resolved population analysis; logical, intended rotational and spectators."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh,expm
from nacs.data import benchmark_path
from nacs.controls.gaussian_echo import (PairModel,frozen,VALID_BASIS,
 rotated_linear_polarization,target_rabi,echo_rabi,envelope,square_gate,
 gate_matrix as gaussian_gate_matrix,nominal_z_calibration,conditional_fidelity,
 pulse_matrix,dwell_matrix,Controls,Uncertainty)
LOGICAL=('0','1')
INTENDED=('0','1','e')
KINDS=('square','truncated')
STAGE_LABELS=('transfer','dwell1','echo','dwell2','retrieve')
STATES=('00','01','10','11')

@dataclass
class Trace:
    times: np.ndarray
    stages: list[str]
    states: np.ndarray # (samples, pairs, 4 inputs)
    checkpoints: list[int]
    model: PairModel
    controls: Controls
    kind: str
    angle: float

    def projectors(self):
        m=self.model;n=m.n
        ls=set(m.logical_indices)
        core={m.pair(a,b) for a in INTENDED for b in INTENDED}
        spectator=set(range(m.size))-core
        rotational={j for j in core-ls if m.n_exc[j//n] or m.n_exc[j%n]}
        # in the three-state subspace core minus logical, all have e
        return {'logical':np.array(sorted(ls),int),
                'intended_rot':np.array(sorted(rotational),int),
                'spectator':np.array(sorted(spectator),int),
                'ground_spectator':np.array(sorted(j for j in spectator
                  if not m.n_exc[j//n] and not m.n_exc[j%n]),int),
                'rot_spectator':np.array(sorted(j for j in spectator
                  if m.n_exc[j//n] or m.n_exc[j%n]),int)}

    def populations(self):
        # (time, input), each is a subspace projection, partition by design
        P=np.abs(self.states)**2
        v={k:P[:,idx,:].sum(axis=1) if len(idx) else np.zeros((len(P),4))
           for k,idx in self.projectors().items()}
        v['total']=P.sum(axis=1)
        return v

    def spectator_by_label(self, time_index=-1):
        """P(at least one molecule occupies label), not pair occupation number.
        Contributions from multiple labels can overlap in two-spectator states.
        """
        P=np.abs(self.states[time_index])**2
        n=self.model.n
        out={}
        for name,index in self.model.local.items():
            if name in INTENDED:continue
            idx=[a*n+b for a in range(n) for b in range(n)
                 if (a==index or b==index)]
            out[name]=np.sum(P[idx],axis=0)
        return out


def initial_states(model):
    return np.eye(model.size,dtype=complex)[:,model.logical_indices]


def _square_stage(model, x, start, duration, kind, control, pol, fr, samples, dd_on=True):
    """Exact within-stage snapshots with same rotating-frame bookkeeping as frozen."""
    if kind=='dwell':
        h=model.Hfree
        delta=0.
    else:
        i,t,det,phase=control
        h,delta,_=model.mw_hamiltonian(i,t,fr,det,pol,phase)
        if not dd_on: h=h-model.Hdd
    ev,V=eigh(h)
    y=(V.conj().T @ (np.exp(+2j*np.pi*delta*start*model.num_exc_diag)[:,None]*x))
    times=np.linspace(0,duration,samples+1)[1:]
    out=[]
    for dt in times:
        u=V@(np.exp(-2j*np.pi*ev*dt)[:,None]*y)
        u=np.exp(-2j*np.pi*delta*(start+dt)*model.num_exc_diag)[:,None]*u
        out.append(u)
    return np.stack(out)


def _gaussian_stage(model,x,start,duration,control,pol,fr,samples,dd_on=True):
    """Midpoint time-ordering precisely matches frozen gaussian_echo.echo_matrix."""
    i,t,det,phase=control
    h,delta,_=model.mw_hamiltonian(i,t,fr,det,pol,phase)
    h0,_,_=model.mw_hamiltonian(i,t,0.,det,pol,phase)
    h1=h-h0
    if not dd_on:h0=h0-model.Hdd
    dt=duration/samples
    xrot=np.exp(+2j*np.pi*delta*start*model.num_exc_diag)[:,None]*x
    out=[]
    for k in range(samples):
        a=float(envelope((k+.5)/samples,'truncated'))
        xrot=expm(-2j*np.pi*(h0+a*h1)*dt)@xrot
        xr=np.exp(-2j*np.pi*delta*(start+(k+1)*dt)*model.num_exc_diag)[:,None]*xrot
        out.append(xr)
    return np.stack(out)


def trace(model,controls,angle=35.,kind='truncated',transfer_pts=22,dwell_pts=25,echo_pts=128,
          dd_off_stages=()):
    if kind not in KINDS:raise ValueError(kind)
    c=controls
    pol=rotated_linear_polarization(angle)
    fr=target_rabi(model,angle,0.)
    fe=echo_rabi(model,angle,0.)
    stages=[('transfer',c.transfer_ms,'pulse',('1','e',c.transfer_detuning_khz,c.transfer_phase),fr,transfer_pts),
            ('dwell1',c.dwell_ms/2,'dwell',None,0,dwell_pts),
            ('echo',c.echo_ms,'pulse',('0','e',c.echo_detuning_khz,c.echo_phase),fe,echo_pts),
            ('dwell2',c.dwell_ms/2,'dwell',None,0,dwell_pts),
            ('retrieve',c.transfer_ms,'pulse',('1','e',c.transfer_detuning_khz,c.reverse_phase),fr,transfer_pts)]
    x=initial_states(model)
    xs=[x.copy()]; tt=[0.]; tags=['initial']; ch=[0];t0=0.
    for name,T,tp,con,amp,count in stages:
        if name=='echo' and kind=='truncated':
            series=_gaussian_stage(model,x,t0,T,con,pol,amp,count,dd_on=name not in dd_off_stages)
        else:
            series=_square_stage(model,x,t0,T,tp,con,pol,amp,count,dd_on=name not in dd_off_stages)
        for k in range(count):
            xs.append(series[k]);tt.append(t0+(k+1)*T/count);tags.append(name)
        t0+=T;x=series[-1];ch.append(len(xs)-1)
    return Trace(np.asarray(tt),tags,np.asarray(xs),ch,model,c,kind,angle)


def archived_controls(tag='square'):
    if tag=='square':
        fn=benchmark_path('coherent/model_a_nominal_angle35.json')
        # Model A is extracted lazily by the Gaussian package import.
        return Controls(**json.loads(fn.read_text())['controls'])
    elif tag=='gaussian':
        fn=benchmark_path('coherent/gaussian_joint_35.json')
        return Controls(**json.loads(fn.read_text())['joint_echo_controls'])
    raise ValueError(tag)


def archived_final(model,control,tag):
    if tag=='square':
        m=square_gate(model,control,35.)
    else:
        m=gaussian_gate_matrix(model,control,35.,steps=128,kind='truncated')
    return m
