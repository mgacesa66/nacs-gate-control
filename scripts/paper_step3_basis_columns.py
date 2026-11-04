"""High-basis 25-state validation using exact Hermitian column propagator."""
from __future__ import annotations
import json,time
from pathlib import Path
import numpy as np
from nacs.analysis.column_propagator import gate_columns
from nacs.analysis.basis_convergence import build_basis
from nacs.analysis.parameter_recalibration import PROBES,tracked_spectrum
from nacs.analysis.matched_objectives import PairModel,design_from_vector
from nacs.controls.spectator import measures
from nacs.data import ROOT

def main():
 root=ROOT/'outputs/paper_step3_basis';root.mkdir(parents=True,exist_ok=True)
 f,_=tracked_spectrum(PROBES[0]);extra=json.loads((ROOT/'benchmarks/data/paper_step2/paper_step2_union/selection.json').read_text())['union']
 model=PairModel(f,build_basis(f,PROBES[0],extra))
 rows=json.loads((ROOT/'outputs/paper_step3_final_validation/nominal_11.json').read_text())
 for item in rows:
  name=item['protocol'];path=root/f'{name}_25_columns_steps24.json'
  if path.exists():print('PRESERVED',name,flush=True);continue
  t=time.monotonic();out=gate_columns(model,design_from_vector(item['design']),steps=24)
  metric=measures(model,out);m=metric.pop('M');z=metric.pop('z')
  report={k:float(v) for k,v in metric.items()}
  report.update(protocol=name,basis_n=model.n,steps=24,pair_dimension=model.size,
     elapsed_s=time.monotonic()-t,propagator='Hermitian eigen-column propagation',
     logical_matrix_real=m.real.tolist(),logical_matrix_imag=m.imag.tolist(),
     virtual_z_rad=z.tolist())
  old=root/f'{name}_25_steps24.json'
  if old.exists():
   prior=json.loads(old.read_text())
   prev=np.array(prior['logical_matrix_real'])+1j*np.array(prior['logical_matrix_imag'])
   report['dense_expm_max_complex_diff']=float(np.max(abs(prev-m)))
   assert report['dense_expm_max_complex_diff']<1e-8
  tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(report,indent=2)+'\n');tmp.replace(path)
  print('SAVED',name,'F',round(report['fidelity']*100,6),'spectators',round(report['total_spectators']*100,5),
      'sec',round(report['elapsed_s'],1),'dense_check',report.get('dense_expm_max_complex_diff'),flush=True)
if __name__=='__main__':main()
