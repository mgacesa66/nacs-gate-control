"""Selected enlarged-basis recalibration checks for four ellipticity probes.
Only frequency+Rabi tracked mode, all controls frozen. No plots.
"""
import argparse,csv,json
from pathlib import Path
from nacs.analysis.basis_convergence import control_designs,combined_ranking,model_gate,compare_to_reference
from nacs.analysis.parameter_recalibration import PROBES
from nacs.data import ROOT

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step2_perturbed_highbasis');p.add_argument('--steps',type=int,default=64);opt=p.parse_args()
 opt.output.mkdir(exist_ok=True,parents=True)
 rank,_=combined_ranking(PROBES)
 seed=[38,69,40,72];ordered=seed+[r['old_index'] for r in rank if r['old_index'] not in seed]
 ds={k:v for k,v in control_designs().items() if k in ('unmodulated','aux51','all_spectators')}
 rows=[]
 for probe in (PROBES[5],PROBES[6],PROBES[9],PROBES[10]):
  rr0,mm0=model_gate(probe,'frequency_rabi',[],ds,steps=opt.steps)
  base={r['protocol']:r for r in rr0}
  for n in (11,21,23):
   if n==11:rs,ms=rr0,mm0
   else:rs,ms=model_gate(probe,'frequency_rabi',ordered[:n-11],ds,steps=opt.steps)
   for r in rs:
    rows.append({**r,**compare_to_reference(base[r['protocol']],r,mm0[r['protocol']],ms[r['protocol']])})
   print(probe.name,n,[(r['protocol'],round(r['fidelity'],7)) for r in rs],flush=True)
 with (opt.output/'ellipticity_highbasis.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)
 (opt.output/'metadata.json').write_text(json.dumps({'selection':ordered[:12], 'steps':opt.steps,'method':'frozen controls, 11/21/23 pair-selected states with ideal tracking'},indent=2)+'\n')
if __name__=='__main__':main()
