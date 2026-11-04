"""Cross-screen 25-state check: simultaneous MW- and pair-dipolar-ranked omissions."""
import argparse,csv,json
from pathlib import Path
from nacs.analysis.basis_convergence import (model_gate,control_designs,screen_dipolar_one,combined_ranking,compare_to_reference)
from nacs.analysis.parameter_recalibration import PROBES
from nacs.data import ROOT

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step2_union');p.add_argument('--steps',type=int,default=64);opt=p.parse_args()
 opt.output.mkdir(parents=True,exist_ok=True)
 rank,_=combined_ranking(PROBES);seed=[38,69,40,72]
 mw=seed+[z['old_index'] for z in rank if z['old_index'] not in seed]
 all_dd={}
 for pr in PROBES:
  for row in screen_dipolar_one(pr):
   j=row['old_index'];all_dd[j]=max(all_dd.get(j,0),row['score_dd'])
 dd=sorted(all_dd,key=all_dd.get,reverse=True)
 union=list(dict.fromkeys(mw[:10]+dd[:8]))
 print('MW',mw[:10],'DD',dd[:8],'UNION',union,flush=True)
 ds={k:v for k,v in control_designs().items() if k in ('unmodulated','aux51','all_spectators')}
 rows=[]
 for probe in (PROBES[0],PROBES[10]):
  ref,refM=model_gate(probe,'frequency_rabi',[],ds,steps=opt.steps)
  base={x['protocol']:x for x in ref}
  for title,extra in [('11',[]),('21_MW',mw[:10]),('21_DD',dd[:10]),('union',union)]:
   r,m=(ref,refM) if not extra else model_gate(probe,'frequency_rabi',extra,ds,steps=opt.steps)
   for x in r:
    rows.append({'scheme':title,'probe':probe.name,'basis_n':len(extra)+11,'protocol':x['protocol'],
      'fidelity':x['fidelity'],'logical_leakage':x['logical_leakage'],'total_spectators':x['total_spectators'],
      'aux51':x['aux51'],'other_spectators':x['other_spectators'],
      **compare_to_reference(base[x['protocol']],x,refM[x['protocol']],m[x['protocol']])})
   print(probe.name,title,'N',11+len(extra),[(x['protocol'],round(x['fidelity'],8)) for x in r],flush=True)
 with (opt.output/'cross_screen.csv').open('w',newline='') as fp:
  w=csv.DictWriter(fp,rows[0].keys());w.writeheader();w.writerows(rows)
 (opt.output/'selection.json').write_text(json.dumps({'microwave':mw[:16],'dipolar':dd[:16],'union':union,'steps':opt.steps},indent=2)+'\n')
if __name__=='__main__':main()
