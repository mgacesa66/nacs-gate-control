"""Step 2: field-aware spectral screening and genuine pair-dynamics basis checks.

Spectral scores are heuristics, not mathematical omission-error bounds.
States are tracked using a global one-to-one maximal-overlap assignment;
near avoided crossings retain both partners in the selected subspace.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.optimize import linear_sum_assignment
from nacs.analysis.parameter_recalibration import (PROBES,nominal_objects,tracked_spectrum,
        drive_parameters,build_echo,full_gate_output,_params)
from nacs.spectroscopy.core import NaCsMolecule
from nacs.controls.fixed_polarization import VALID_BASIS,rotated_linear_polarization
from nacs.gate.model import PairModel
from nacs.controls.spectator import measures,Design
from nacs.gate.propagate import TRANSFER_RABI_KHZ, ECHO_RABI_KHZ, TRANSFER_EPS,ECHO_EPS
from nacs.data import benchmark_path
import json

@lru_cache(maxsize=16)
def assignments(probe):
    _,v0,_=nominal_objects()
    v=NaCsMolecule(_params(probe)).diagonalize().vectors
    W=np.abs(v0.conj().T@v)**2
    x,y=linear_sum_assignment(-W)
    perm=np.empty(128,dtype=int);perm[x]=y
    return perm,W

def _intensity_prefactors():
    base,_,_=nominal_objects()
    eps=np.asarray(TRANSFER_EPS,dtype=complex);eps/=np.linalg.norm(eps)
    epse=np.asarray(ECHO_EPS,dtype=complex)
    def dip(pol,a,b):
        return abs(sum(pol[k]*base.dipoles[k][base.indices[a],base.indices[b]] for k in range(3)))
    return (TRANSFER_RABI_KHZ/dip(eps,'1','e'),ECHO_RABI_KHZ/dip(epse,'0','e'))

def screen_one(probe, detuning_floor_khz=10.):
    """Rank omitted levels by maximal single-photon microwave mixing estimate.

    Includes N=0 and N=1 candidates, three core source states and both
    physical microwave electric fields. The floor avoids an unbounded score
    near resonance; it does not establish an error bound or guarantee
    completeness under two-photon/dipolar processes.
    """
    f,_=tracked_spectrum(probe)
    perm,W=assignments(probe)
    pol=rotated_linear_polarization(35.)
    ptr,pec=_intensity_prefactors()
    core=('0','1','e')
    nu_tr=f.energies[f.indices['e']]-f.energies[f.indices['1']]
    nu_ec=f.energies[f.indices['e']]-f.energies[f.indices['0']]
    old,_,_=nominal_objects()
    omitted=set(range(128))-set(old.indices[k] for k in VALID_BASIS)
    out=[]
    for old_index in omitted:
        index=int(perm[old_index]);score=-1.;best=None
        for a in core:
            i=f.indices[a]
            if (index>=32)==(i>=32):continue
            d=abs(sum(pol[k]*f.dipoles[k][i,index] for k in range(3)))
            for drive,nu,fper in [('transfer',nu_tr,ptr),('echo',nu_ec,pec)]:
                det=abs(f.energies[index]-f.energies[i])-nu
                # amplitude normalized as target full Rabi (H_offdiag = f/2)
                g=fper*d/2
                x=g/max(abs(det),detuning_floor_khz)
                if x>score:score=x;best=(a,drive,float(det),float(d))
        out.append(dict(probe=probe.name,old_index=old_index,new_index=index,
              score=float(score),reference_overlap=float(W[old_index,index]),
              source=best[0] if best else '',drive=best[1] if best else '',
              detuning_khz=best[2] if best else float('nan'),
              dipole_debye=best[3] if best else 0.))
    return sorted(out,key=lambda z:z['score'],reverse=True)

def combined_ranking(probes):
    records=[]
    for p in probes: records+=screen_one(p)
    agg={}
    for r in records:
        i=r['old_index']
        if i not in agg or r['score']>agg[i]['score']:agg[i]=r.copy()
    return sorted(agg.values(),key=lambda z:z['score'],reverse=True),records

def build_basis(f,probe,extra_indices):
    perm,W=assignments(probe)
    f.indices.update({f'spec{i}':int(perm[i]) for i in extra_indices})
    return VALID_BASIS+tuple(f'spec{i}' for i in extra_indices)

def control_designs():
    a=json.loads(benchmark_path('optimization/candidates.json').read_text())
    return {name:Design(**v['design']) for name,v in a.items()}

def model_gate(probe,mode,basis_ids,designs=None,steps=64):
    f,_=tracked_spectrum(probe)
    labels=build_basis(f,probe,basis_ids)
    model=PairModel(f,labels)
    designs=designs or control_designs()
    first=next(iter(designs.values()))
    p0=drive_parameters(f,mode,first)
    echo=build_echo(model,p0,steps=steps)
    rows=[]; mats={}
    for label,d in designs.items():
        param=drive_parameters(f,mode,d)
        out=full_gate_output(model,d,param,echo)
        q=measures(model,out)
        M=q.pop('M');z=q.pop('z')
        mats[label]=M
        rows.append({'probe':probe.name,'mode':mode,'basis_n':len(labels),
          'pair_dimension':model.size,'protocol':label,
          'fidelity':q['fidelity'],'logical_leakage':q['leakage'],
          'total_spectators':q['total_spectators'],'aux51':q['aux51'],
          'other_spectators':q['other_spectators'],
          'intended_nonlogical':q['intended_nonlogical'],
          'virtual_z1':float(z[0]),'virtual_z2':float(z[1]),
          'unitarity_error':float(np.max(np.abs(np.sum(np.abs(out)**2,axis=0)-1.))),
          'extra_indices':','.join(map(str,basis_ids))})
    return rows,mats


def compare_to_reference(reference_row,row,ref_M,M):
    return {'delta_fidelity_pp':100*(row['fidelity']-reference_row['fidelity']),
       'delta_leakage_pp':100*(row['logical_leakage']-reference_row['logical_leakage']),
       'delta_spectators_pp':100*(row['total_spectators']-reference_row['total_spectators']),
       'delta_aux51_pp':100*(row['aux51']-reference_row['aux51']),
       'logical_operator_frobenius_diff':float(np.linalg.norm(M-ref_M)),
       'logical_operator_max_complex_entry_diff':float(np.max(np.abs(M-ref_M)))}

def screen_dipolar_one(probe, energy_floor_khz=.1):
    """Independent pair-interaction-based candidate screen.

    For an omitted state j, rank max |<j,c|Hdd|a,b>|/|E_j+E_c-E_a-E_b|,
    where a,b,c are in the *corrected 11-state* active set and total N=1
    count is conserved. This is a first-order mixing HEURISTIC, not a bound.
    """
    from nacs.gate.model import DEBYE
    from scipy.constants import h,epsilon_0,pi
    f,_=tracked_spectrum(probe)
    perm,_=assignments(probe)
    selected=np.array([f.indices[k] for k in VALID_BASIS],dtype=int)
    omitted=set(range(128))-set(nominal_objects()[0].indices[k] for k in VALID_BASIS)
    kdd=(DEBYE**2)/(4*pi*epsilon_0*(1.9e-6)**3)/h/1e3
    # x-directed separation: tensor diag(-2,1,1)
    d=f.dipoles;energies=f.energies;nc=np.arange(128)>=32
    rows=[]
    for old in omitted:
        j=int(perm[old]);best=(-1.,None)
        for a in selected:
            for b in selected:
                for c in selected:
                    if int(nc[a])+int(nc[b])!=int(nc[j])+int(nc[c]):continue
                    J=kdd*(-2*d[0][a,j]*d[0][b,c]+d[1][a,j]*d[1][b,c]+d[2][a,j]*d[2][b,c])
                    det=float(energies[j]+energies[c]-energies[a]-energies[b])
                    score=abs(J)/max(energy_floor_khz,abs(det))
                    if score>best[0]:best=(score,(a,b,c,det,abs(J)))
        detail=best[1]
        rows.append({'probe':probe.name,'old_index':old,'score_dd':float(best[0]),
             'source_a':int(detail[0]) if detail else -1,
             'source_b':int(detail[1]) if detail else -1,
             'partner_c':int(detail[2]) if detail else -1,
             'energy_mismatch_khz':float(detail[3]) if detail else float('nan'),
             'coupling_khz':float(detail[4]) if detail else 0.})
    return sorted(rows,key=lambda x:x['score_dd'],reverse=True)
