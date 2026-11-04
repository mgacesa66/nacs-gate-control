"""Additional numerical diagnostic: logical-map change with and without Z/global phase alignment."""
from pathlib import Path
import json
import numpy as np
from nacs.analysis.basis_convergence import model_gate,control_designs,combined_ranking,screen_dipolar_one
from nacs.analysis.parameter_recalibration import PROBES
from nacs.controls.polarization import nominal_z_calibration
from nacs.data import ROOT


def corrected(M):
    a,b=nominal_z_calibration(M)
    f=np.exp(.5j*np.array([-a-b,-a+b,a-b,a+b]))
    return f[:,None]*M

def main():
 out=ROOT/'outputs/paper_step2_alignment';out.mkdir(exist_ok=True,parents=True)
 rank,_=combined_ranking(PROBES); seed=[38,69,40,72]
 mw=seed+[r['old_index'] for r in rank if r['old_index'] not in seed]
 dd={}
 for probe in PROBES:
  for r in screen_dipolar_one(probe):dd[r['old_index']]=max(dd.get(r['old_index'],0),r['score_dd'])
 ddrank=sorted(dd,key=dd.get,reverse=True)
 union=list(dict.fromkeys(mw[:10]+ddrank[:8])); assert len(union)==14
 ds={k:v for k,v in control_designs().items() if k in ('unmodulated','aux51','all_spectators')}
 rows=[]; store={}
 for probe in (PROBES[0],PROBES[10]):
  ref,rm=model_gate(probe,'frequency_rabi',[],ds,steps=32)
  new,nm=model_gate(probe,'frequency_rabi',union,ds,steps=32)
  for key in ds:
   A=rm[key];B=nm[key]
   A2=corrected(A);B2=corrected(B)
   phase=np.angle(np.vdot(A2,B2));B2*=np.exp(-1j*phase)
   rows.append({'probe':probe.name,'protocol':key,'basis_n':len(union)+11,
       'raw_fro':float(np.linalg.norm(B-A)),
       'z_and_global_corrected_fro':float(np.linalg.norm(B2-A2)),
       'raw_max_entry':float(np.max(np.abs(B-A))),
       'aligned_max_entry':float(np.max(np.abs(B2-A2))),
       'z_ref':nominal_z_calibration(A).tolist(),'z_new':nominal_z_calibration(B).tolist()})
   store[probe.name+'_'+key+'_n11']=A
   store[probe.name+'_'+key+'_n25']=B
   print(rows[-1],flush=True)
 (out/'operator_alignment.json').write_text(json.dumps(rows,indent=2)+'\n')
 np.savez_compressed(out/'logical_maps.npz',**store)
if __name__=='__main__':main()
