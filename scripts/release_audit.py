"""Release-time numerical audit. No plots, refitting, or optimization.

Run: python scripts/release_audit.py --mode all --output outputs/release_audit
All calculations use frozen baseline data. Parameter perturbations use the
same microscopic constants with small, explicitly hypothetical variations.
"""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import argparse,json,csv,hashlib
import numpy as np

from nacs.data import frozen_spectrum,benchmark_path
from nacs.gate.model import PairModel, FrozenSpectrum, BASIS_SETS
from nacs.controls.fixed_polarization import VALID_BASIS,EXTRA_BASIS
from nacs.controls.spectator import Design,gate_output,measures,echo_ref_matrix,_pair_rot_frame,ANGLE,BASE
from nacs.spectroscopy.constants import NaCsParameters
from nacs.spectroscopy.core import NaCsMolecule
from scipy.optimize import linear_sum_assignment
from scipy.linalg import expm
from nacs.controls.gaussian_echo import echo_matrix
from nacs.controls.fixed_polarization import rotated_linear_polarization
from nacs.gate.propagate import ECHO_RABI_KHZ,TRANSFER_RABI_KHZ,ECHO_EPS,TRANSFER_EPS

ROOT=Path(__file__).resolve().parents[1]

def rows_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def benchmark_audit(out):
    tab=list(csv.DictReader(benchmark_path('spectroscopy/seven_line_accounting.csv').open()))
    rows=[]
    for r in tab:
        measured=float(r['experiment_khz']); raw=float(r['static_model_khz']);adjusted=float(r['adjusted_model_khz'])
        rows.append({'line':r['line'],'measured_khz':measured,'static_model_khz':raw,
          'raw_difference_khz':round(raw-measured,5),'posthoc_accounted_difference_khz':round(adjusted-measured,5),
          'is_independently_predicted':False if abs(float(r['inferred_driven_shift_khz']))>0 else None,
          'note':'finite-power amplitudes were inferred, NOT independently measured'})
    rows_csv(out/'spectroscopy_accounting.csv',rows)
    return rows

def basis_audit(out):
    # For the release, verify complete projected 4x4 channel with identical controls.
    f=frozen_spectrum(); archived=json.loads(benchmark_path('optimization/candidates.json').read_text())
    # Additional candidates chosen previously from frozen couplings; basis is deliberately not the full 128.
    extra={'rank38':38,'rank69':69,'rank40':40,'rank72':72}
    f.indices.update(extra)
    bases={'legacy9':BASIS_SETS['nine'], 'corrected11':VALID_BASIS,
           'corrected13':EXTRA_BASIS,
           'ranked15':VALID_BASIS+tuple(extra)}
    rows=[]; cached={}
    for name,basis in bases.items():
        model=PairModel(f,basis)
        echo=echo_ref_matrix(model,steps=64 if len(basis)<=13 else 48)
        for protocol in ['unmodulated','gate_only','aux51','all_spectators']:
            d=Design(**archived[protocol]['design'])
            result=measures(model,gate_output(model,d,echo_ref=echo))
            rows.append({'basis':name,'n_one_molecule':len(basis),'n_pair':len(basis)**2,'protocol':protocol,
              **{k:float(result[k]) for k in ['fidelity','leakage','total_spectators','aux51','other_spectators']}})
            cached[(name,protocol)]=rows[-1]
    rows_csv(out/'basis_check.csv',rows)
    return rows

def sensitivity_audit(out):
    # Small physical perturbations, no parameter refitting; eigenstates tracked by overlap.
    f=frozen_spectrum()
    d=np.load(ROOT/'benchmarks/frozen/dressed_128.npz')
    oldvec=d['vectors']; basem=json.loads((ROOT/'benchmarks/frozen/baseline_manifest.json').read_text())
    basefield=basem['magnetic_field_G']; tensor=basem['optical_tensor_khz']; chi0=basem['optical_ellipticity_deg']
    cases=[('nominal',0.,0.,0.),('B_minus_0p1G',-.1,0.,0.),('B_plus_0p1G',.1,0.,0.),
           ('optical_tensor_minus_0p5pct',0.,-.005,0.),('optical_tensor_plus_0p5pct',0.,.005,0.),
           ('ellipticity_minus_0p5deg',0.,0.,-.5),('ellipticity_plus_0p5deg',0.,0.,.5)]
    # Report specific intrinsic observables plus full-gate changes with *fixed absolute carriers*.
    nominal_nu={key:float(f.energies[f.indices[target]]-f.energies[f.indices[src]]) for key,src,target in [('transfer','1','e'),('echo','0','e')]}
    archived=json.loads(benchmark_path('optimization/candidates.json').read_text())
    baseline_design=Design(**archived['all_spectators']['design'])
    rows=[]
    for name,bshift,tfactor,eshift in cases:
        chi=np.deg2rad(chi0+eshift)
        p=NaCsParameters(field_G=basefield+bshift,optical_tensor_khz=tensor*(1+tfactor),
           polarization_xyz=(complex(np.cos(chi)),complex(0,np.sin(chi)),0j),dipole_debye=4.6)
        sp=NaCsMolecule(p).diagonalize()
        # Atomic-basis vector overlap; label tracking is NOT energy-index tracking.
        # Global one-to-one overlap assignment avoids relabeling collisions.
        S=oldvec.conj().T@sp.vectors
        iold,inew=linear_sum_assignment(-abs(S)**2)
        assignment=np.empty(128,dtype=int);assignment[iold]=inew
        overlap=np.abs(S[np.arange(128),assignment])**2
        selected=list(dict.fromkeys(list(f.indices.values())+[51,79,78,80,38,69,40,72]))
        copied={k:int(assignment[i]) for k,i in f.indices.items()}
        # Align all dressed-state phases to the original eigenvectors. Without
        # this, changes in the basis gauge look like physical gate phase errors.
        V=sp.vectors.copy()
        for old_index,new_index in enumerate(assignment):
            inner=S[old_index,new_index]
            if abs(inner): V[:,new_index]*=np.conj(inner)/abs(inner)
        bare=NaCsMolecule(p).bare_dipole_q_debye()
        D={q:V.conj().T@bare[q]@V for q in (-1,0,1)}
        dm,dz,dp=(D[q] for q in (-1,0,1))
        fnew=FrozenSpectrum(sp.energies_khz,((dm-dp)/np.sqrt(2),1j*(dm+dp)/np.sqrt(2),dz),copied,'recomputed_unfrozen')
        mi,mj=copied['1'],copied['e'];aux=copied['aux51']
        mw=np.array([.80536,.57825,.13047])
        dip=lambda a,b:abs(sum(mw[k]*fnew.dipoles[k][a,b] for k in range(3)))
        # With the same actual MW electric field (NOT the same re-normalized Rabi
        # frequency), evaluate the corrected 11-state gate using absolute carrier
        # frequencies held constant. No reoptimization or spectroscopy fit.
        m=PairModel(fnew,VALID_BASIS)
        pol=rotated_linear_polarization(35.)
        f0=frozen_spectrum()
        f0.indices.update({'aux51':51,'aux79':79,'extra78':78,'extra80':80})
        refi,refe=f0.indices['1'],f0.indices['e'];ref0=f0.indices['0']
        amp_nom_trans=abs(sum(np.asarray(TRANSFER_EPS)[k]*f0.dipoles[k][refi,refe] for k in range(3)))
        amp_trans=abs(sum(pol[k]*fnew.dipoles[k][mi,mj] for k in range(3)))
        ftrans=TRANSFER_RABI_KHZ*amp_trans/amp_nom_trans
        amp_nom_echo=abs(sum(np.asarray(ECHO_EPS)[k]*f0.dipoles[k][ref0,refe] for k in range(3)))
        amp_echo=abs(sum(pol[k]*fnew.dipoles[k][copied['0'],mj] for k in range(3)))
        fecho=ECHO_RABI_KHZ*amp_echo/amp_nom_echo
        dt=baseline_design.transfer_detuning_khz+nominal_nu['transfer']-(fnew.energies[mj]-fnew.energies[mi])
        de=BASE.echo_detuning_khz+nominal_nu['echo']-(fnew.energies[mj]-fnew.energies[copied['0']])
        trans=_pair_rot_frame(m,'1','e',ftrans,dt,pol,0.,baseline_design.phase_forward,
           baseline_design.transfer_ms,baseline_design.fraction_first,0.)
        dwell=expm(-2j*np.pi*m.Hfree*(baseline_design.dwell_ms/2))
        ec=replace(BASE,transfer_ms=baseline_design.transfer_ms,dwell_ms=baseline_design.dwell_ms,
             echo_detuning_khz=de)
        # echo_matrix ordinarily uses model-based amplitude normalization, but here
        # we need the physical fixed E-field. Correct via a model-local wrapper
        # that multiplies the microwave drive by fecho / echo_rabi(model,35).
        from nacs.controls.fixed_polarization import echo_rabi
        from nacs.controls.polarization import Uncertainty
        ratio=fecho/echo_rabi(m,35.)
        echo=echo_matrix(m,35.,ec,err=Uncertainty(amp=ratio-1.),steps=64)
        reverse=_pair_rot_frame(m,'1','e',ftrans,dt,pol,
           baseline_design.reverse_phase+baseline_design.phase_reverse,
           baseline_design.reverse_phase,baseline_design.transfer_ms,
           baseline_design.fraction_first,
           baseline_design.transfer_ms+baseline_design.dwell_ms+BASE.echo_ms)
        inputs=np.eye(m.size,dtype=complex)[:,m.logical_indices]
        q=measures(m,reverse@(dwell@(echo@(dwell@(trans@inputs)))))
        rows.append({'scenario':name,'field_G':p.field_G,'optical_tensor_khz':p.optical_tensor_khz,
           'ellipticity_deg':chi0+eshift,'ground_split_1minus0_khz':float(sp.energies_khz[copied['1']]-sp.energies_khz[copied['0']]),
           'transfer_frequency_khz':float(sp.energies_khz[mj]-sp.energies_khz[mi]),
           'aux51_detuning_from_transfer_khz':float(sp.energies_khz[aux]-sp.energies_khz[mj]),
           'd1e_debye':float(dip(mi,mj)),'d1aux51_debye':float(dip(mi,aux)),
           'min_selected_tracking_overlap':float(min(overlap[i] for i in selected)),
           'nominal_fixed_carrier_fidelity':float(q['fidelity']),
           'nominal_fixed_carrier_leakage':float(q['leakage']),
           'no_recalibration':True})
    rows_csv(out/'parameter_sensitivity.csv',rows)
    return rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['all','spectroscopy','basis','sensitivity'],default='all')
    ap.add_argument('--output',type=Path,default=ROOT/'outputs'/'release_audit')
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    output={}
    if args.mode in ['all','spectroscopy']:output['spectroscopy']=benchmark_audit(args.output)
    if args.mode in ['all','basis']:output['basis']=basis_audit(args.output)
    if args.mode in ['all','sensitivity']:output['sensitivity']=sensitivity_audit(args.output)
    (args.output/'audit.json').write_text(json.dumps(output,indent=2)+'\n')
    print('audit complete:',{k:len(v) for k,v in output.items()})
    print('outputs:',args.output)
if __name__=='__main__':main()
