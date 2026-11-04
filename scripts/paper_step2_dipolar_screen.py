"""Independent pair-dipolar screen exposing microwave-only screening blind spots."""
import argparse,csv,json
from pathlib import Path
from nacs.analysis.basis_convergence import (screen_dipolar_one,model_gate,control_designs,compare_to_reference)
from nacs.analysis.parameter_recalibration import PROBES
from nacs.data import ROOT

def writecsv(f,rows):
 with f.open('w',newline='') as fd:
  w=csv.DictWriter(fd,rows[0].keys());w.writeheader();w.writerows(rows)

def main():
 a=argparse.ArgumentParser();a.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step2_dipolar');a.add_argument('--steps',type=int,default=64);opt=a.parse_args()
 opt.output.mkdir(parents=True,exist_ok=True)
 probe_list=(PROBES[0],PROBES[5],PROBES[10]);scans=[]
 for p in probe_list:scans.extend(screen_dipolar_one(p))
 scores={}
 for s in scans:scores[s['old_index']]=max(scores.get(s['old_index'],0),s['score_dd'])
 ordered=sorted(scores,key=lambda x:scores[x],reverse=True)
 microwave_seed=[38,69,40,72]
 bases={'11_mw':[], '15_mw':microwave_seed,
        '15_dd':ordered[:4], '19_dd':ordered[:8],
        '19_union':list(dict.fromkeys(microwave_seed+ordered[:4])),
        '21_union':list(dict.fromkeys(microwave_seed+ordered[:8])),
        '23_joint':list(dict.fromkeys(microwave_seed+ordered[:10]))}
 for name in list(bases):
  if len(bases[name])!=int(name[:2])-11:print('WARNING size name',name,len(bases[name])+11,flush=True)
 ds={k:v for k,v in control_designs().items() if k in ('unmodulated','aux51','all_spectators')}
 rows=[]
 for p in probe_list:
  ref,rm=model_gate(p,'frequency_rabi',[],ds,steps=opt.steps);re={r['protocol']:r for r in ref}
  for name,ids in bases.items():
   r,m=(ref,rm) if not ids else model_gate(p,'frequency_rabi',ids,ds,steps=opt.steps)
   for item in r:
    rows.append({'probe':p.name,'basis_scheme':name,'basis_n':11+len(ids),'protocol':item['protocol'],
       **{k:v for k,v in item.items() if k not in ('probe','basis_n','protocol')},
       **compare_to_reference(re[item['protocol']],item,rm[item['protocol']],m[item['protocol']])})
   print(p.name,name,[(x['protocol'],round(x['fidelity'],7)) for x in r],flush=True)
 writecsv(opt.output/'dipolar_screen_all.csv',scans)
 writecsv(opt.output/'dipolar_basis.csv',rows)
 (opt.output/'basis_order.json').write_text(json.dumps({'scores':{str(k):v for k,v in scores.items()},'bases':bases,'steps':opt.steps},indent=2)+'\n')
if __name__=='__main__':main()
