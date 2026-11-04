"""Equal-budget local refinement starting at each arm's preregistered global winner."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from nacs.analysis.matched_objectives import (OBJECTIVES,BOUNDS,SEARCH_BASIS,PairModel,frozen,
    echo_ref_matrix,design_from_vector,gate_output,measures,isolated_transfers,score,feasible,
    PHASE_RAMP_US,RAMP_STEPS,ECHO_STEPS_TRAIN)

def main():
 p=argparse.ArgumentParser();p.add_argument('--objective',choices=OBJECTIVES,required=True)
 p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 p.add_argument('--maxfev',type=int,default=240)
 a=p.parse_args()
 if a.output.exists():print('PRESERVED',a.output,flush=True);return
 records=[json.loads(f.read_text()) for f in a.input.glob(f'{a.objective}_seed*.json')]
 init=min(records,key=lambda x:x['best_objective'])
 m=PairModel(frozen,SEARCH_BASIS);echo=echo_ref_matrix(m,steps=ECHO_STEPS_TRAIN)
 count=0;best=[float('inf'),None]
 def obj(x):
  nonlocal count
  count+=1
  if not feasible(x):return 5.+10*max(0,design_from_vector(x).duration_ms-1.3)
  d=design_from_vector(x)
  if a.objective=='isolated_transfer':
   pf,pr=isolated_transfers(m,d);v=score(a.objective,None,pf,pr)
  else:
   out=gate_output(m,d,echo_ref=echo,steps=ECHO_STEPS_TRAIN,
            phase_ramp_us=PHASE_RAMP_US,ramp_steps=RAMP_STEPS)
   met=measures(m,out);v=score(a.objective,met,0,0)
  if v<best[0]:best[:]=[float(v),list(map(float,x))]
  return v
 r=minimize(obj,init['vector'],method='Powell',bounds=BOUNDS,
           options={'maxfev':a.maxfev,'xtol':1e-4,'ftol':1e-7})
 # Store winner including starting point: local refinement can go worse.
 if init['best_objective']<best[0]:best=[init['best_objective'],init['vector']]
 ret=dict(objective=a.objective,starting_seed=init['seed'],start_objective=init['best_objective'],
   maxfev=a.maxfev,evals=count,objective_value=best[0],vector=best[1],method='Powell',
   termination_message=str(r.message))
 a.output.parent.mkdir(parents=True,exist_ok=True)
 t=a.output.with_suffix('.tmp');t.write_text(json.dumps(ret,indent=2)+'\n');t.replace(a.output)
 print('SAVED',a.output,'obj',ret['objective_value'],'evals',count,'message',r.message,flush=True)
if __name__=='__main__':main()
