"""Fixed-count, deterministic 11-state refinement of equal-budget global winners.

Each objective gets the same 181 full objective evaluations (one starting
candidate + 180 trial controls). Hybrid coordinate/random neighborhood search.
Global optimality is NOT asserted.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from nacs.analysis.matched_objectives import (OBJECTIVES,BOUNDS,PairModel,frozen,
  VALID_BASIS,echo_ref_matrix,gate_output,measures,isolated_transfers,score,
  design_from_vector,PHASE_RAMP_US,RAMP_STEPS,feasible)
from paper_step3_validate import candidate_dict

SCALE=np.array([.010,.035,.35,.18,.18,.16])

def main():
 p=argparse.ArgumentParser();p.add_argument('--objective',required=True,choices=OBJECTIVES)
 p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 p.add_argument('--budget',type=int,default=180);p.add_argument('--steps',type=int,default=24)
 p.add_argument('--seed',type=int,default=9184)
 a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
 if a.output.exists():print('PRESERVED',a.output,flush=True);return
 start=candidate_dict(a.input)[a.objective]
 m=PairModel(frozen,VALID_BASIS)
 echo=echo_ref_matrix(m,steps=a.steps)
 def f(x):
  if not feasible(x):return 10.+max(0,design_from_vector(x).duration_ms-1.30)
  d=design_from_vector(x)
  if a.objective=='isolated_transfer':
   pf,pr=isolated_transfers(m,d);return float(score(a.objective,None,pf,pr))
  out=gate_output(m,d,echo_ref=echo,steps=a.steps,
       phase_ramp_us=PHASE_RAMP_US,ramp_steps=RAMP_STEPS)
  q=measures(m,out)
  return float(score(a.objective,q,0,0))
 lo,hi=np.asarray(BOUNDS).T
 rng=np.random.default_rng(a.seed)
 best=np.asarray(start,float);scorebest=f(best)
 changes=0
 for j in range(a.budget):
  block=j//(a.budget//4 or 1)
  spread=SCALE*([1.,.5,.25,.13][min(block,3)])
  k=j%16
  direction=np.zeros(6)
  if k<12:
   direction[k//2]=(1 if k%2 else -1)
  else:
   direction=rng.normal(size=6)
   direction/=max(np.linalg.norm(direction),1e-12)
  proposed=np.clip(best+direction*spread,lo,hi)
  if design_from_vector(proposed).duration_ms>1.30:
   proposed[1]=min(proposed[1],1.30-2*proposed[0]-0.057922)
  val=f(proposed)
  if val<scorebest-1e-12:
   best,scorebest=proposed,float(val);changes+=1
 data=dict(objective=a.objective,starting_vector=start,final_vector=best.tolist(),
   final_objective=float(scorebest),n_objective_evaluations=a.budget+1,
   trial_budget=a.budget,steps=a.steps,refinement_basis_n=len(VALID_BASIS),
   phase_ramp_us=PHASE_RAMP_US,seed=a.seed,accepted_updates=changes,
   method='deterministic-matched-pattern-search',not_global_optimum=True)
 temp=a.output.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(a.output)
 print('SAVED',a.output,'objective',round(scorebest,9),'updates',changes,'evaluations',a.budget+1,flush=True)
if __name__=='__main__':main()
