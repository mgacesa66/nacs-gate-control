"""One checkpoint per objective and seed; no destructive writes to earlier data."""
from __future__ import annotations
import argparse,json,os
from pathlib import Path
from nacs.analysis.matched_objectives import OBJECTIVES,optimize_arm
from nacs.data import ROOT

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--objective',required=True,choices=OBJECTIVES)
 p.add_argument('--seed',type=int,required=True)
 p.add_argument('--iterations',type=int,default=7)
 p.add_argument('--population',type=int,default=30)
 p.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step3')
 a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=True)
 dst=a.output/f'{a.objective}_seed{a.seed}.json'
 if dst.exists():
  obj=json.loads(dst.read_text())
  if obj['iterations']==a.iterations and obj['population']==a.population:
   print('PRESERVED',dst,flush=True);return
  raise ValueError('Existing checkpoint has different budget; choose another output directory')
 res=optimize_arm(a.objective,a.seed,a.iterations,a.population)
 temporary=dst.with_suffix('.pending.json')
 temporary.write_text(json.dumps(res,indent=2)+'\n')
 temporary.replace(dst)
 print('SAVED',dst,'value',res['best_objective'],'calls',res['evaluations'],flush=True)
if __name__=='__main__':main()
