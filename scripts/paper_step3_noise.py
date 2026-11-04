"""Identical held-out error realizations, nominal Z held fixed per protocol.

Independent independent quasi-static errors; not a measured hardware noise model.
Writes one completed protocol/noise-scenario checkpoint at a time.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from nacs.analysis.matched_objectives import (PairModel,frozen,VALID_BASIS,
    echo_ref_matrix,gate_output,conditional_fidelity,design_from_vector,
    PHASE_RAMP_US,RAMP_STEPS)
from nacs.data import ROOT

SCALES={'tight':(.005,.05,.01),'broad':(.02,.5,.02)}

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--input',type=Path,default=ROOT/'outputs/paper_step3_final_validation/nominal_11.json')
 p.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step3_noise')
 p.add_argument('--protocol',required=True)
 p.add_argument('--scenario',choices=tuple(SCALES),required=True)
 p.add_argument('--samples',type=int,default=24)
 p.add_argument('--steps',type=int,default=24);p.add_argument('--seed',type=int,default=32094)
 a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=True)
 fname=a.output/f'{a.protocol}_{a.scenario}.csv'
 if fname.exists():print('PRESERVED',fname,flush=True);return
 entries={x['protocol']:x for x in json.loads(a.input.read_text())}
 ref=entries[a.protocol]
 z=np.asarray(ref['virtual_z_rad'],float)
 design=design_from_vector(ref['design'])
 rng=np.random.default_rng(a.seed)
 raw=rng.uniform(-1,1,(a.samples,3))
 errs=raw*np.asarray(SCALES[a.scenario])[None,:]
 rows=[]
 for k,(amp,det,exch) in enumerate(errs):
  model=PairModel(frozen,VALID_BASIS,scale_dd=1+exch)
  echo=echo_ref_matrix(model,amp_error=float(amp),steps=a.steps)
  out=gate_output(model,design,amp_error=float(amp),det_error_khz=float(det),
       echo_ref=echo,steps=a.steps,phase_ramp_us=PHASE_RAMP_US,ramp_steps=RAMP_STEPS)
  mat=out[np.asarray(model.logical_indices),:]
  rows.append(dict(sample=k,amp_error=float(amp),det_error_khz=float(det),
      exchange_error=float(exch),fidelity=float(conditional_fidelity(mat,z)),
      logical_survival=float(np.sum(abs(mat)**2).real/4)))
  if (k+1)%8==0:print(a.protocol,a.scenario,'samples',k+1,flush=True)
 tmp=fname.with_suffix('.pending.csv')
 with tmp.open('w',newline='') as fd:
  w=csv.DictWriter(fd,rows[0].keys());w.writeheader();w.writerows(rows)
 tmp.replace(fname)
 vals=np.array([r['fidelity'] for r in rows])
 print('SAVED',fname,'mean',float(vals.mean()),'p10',float(np.percentile(vals,10)),flush=True)
if __name__=='__main__':main()
