"""Validate selected matched-control finalists on the independent MW+dipolar 25-state basis."""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np
from nacs.analysis.basis_convergence import build_basis
from nacs.analysis.parameter_recalibration import PROBES,tracked_spectrum
from nacs.analysis.matched_objectives import (PairModel,VALID_BASIS,frozen,echo_ref_matrix,
    measure_candidate,gate_output,measures,design_from_vector,PHASE_RAMP_US,RAMP_STEPS)
from nacs.data import ROOT

def main():
 p=argparse.ArgumentParser();p.add_argument('--basis',choices=['11','15','25'],required=True)
 p.add_argument('--input',type=Path,default=ROOT/'outputs/paper_step3_final_validation/nominal_11.json')
 p.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step3_basis')
 p.add_argument('--steps',type=int,default=24);a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=True)
 entries=json.loads(a.input.read_text())
 if a.basis=='11':m=PairModel(frozen,VALID_BASIS)
 else:
  f,_=tracked_spectrum(PROBES[0])
  u=json.loads((ROOT/'benchmarks/data/paper_step2/paper_step2_union/selection.json').read_text())['union']
  extra=[38,69,40,72] if a.basis=='15' else u
  m=PairModel(f,build_basis(f,PROBES[0],extra))
 print('CREATED basis',a.basis,'states',m.size,'pair dimension',flush=True)
 t=time.monotonic(); echo=echo_ref_matrix(m,steps=a.steps)
 print('ECHO basis',a.basis,'steps',a.steps,'sec',round(time.monotonic()-t,2),flush=True)
 for item in entries:
  kind=item['protocol'];outfn=a.output/f'{kind}_{a.basis}_steps{a.steps}.json'
  if outfn.exists():print('PRESERVED',outfn,flush=True);continue
  x=item['design'];t=time.monotonic()
  out=gate_output(m,design_from_vector(x),echo_ref=echo,steps=a.steps,
     phase_ramp_us=PHASE_RAMP_US,ramp_steps=RAMP_STEPS)
  q=measures(m,out)
  z=q.pop('z');M=q.pop('M')
  values={k:float(v) for k,v in q.items()}
  values.update(protocol=kind,steps=a.steps,basis_n=m.n,pair_dimension=m.size,
      state_index_selection=m.source_indices.tolist(),elapsed_s=round(time.monotonic()-t,2),
      logical_matrix_real=M.real.tolist(),logical_matrix_imag=M.imag.tolist(),
      virtual_z_rad=z.tolist())
  tmp=outfn.with_suffix('.tmp');tmp.write_text(json.dumps(values,indent=2)+'\n');tmp.replace(outfn)
  print('SAVED',kind,'basis',a.basis,'F',round(values['fidelity']*100,5),
        'spec',round(values['total_spectators']*100,4), 'seconds',values['elapsed_s'],flush=True)
if __name__=='__main__':main()
