"""Post-search 11-state validation; same candidate pool, no new optimization."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from nacs.analysis.matched_objectives import (OBJECTIVES,measure_candidate,PairModel,frozen,
    VALID_BASIS,echo_ref_matrix,score,archived_start_vectors,PHASE_RAMP_US)
from nacs.data import ROOT

def candidate_dict(src,include_secondary=True):
 ret={}
 for objective in OBJECTIVES:
  rs=[json.loads(f.read_text()) for f in src.glob(f'{objective}_seed*.json')]
  r=min(rs,key=lambda x:x['best_objective'])
  refine=src.parent/'paper_step3_refine'/f'{objective}.json'
  if refine.exists():
   q=json.loads(refine.read_text())
   if q['objective_value']<r['best_objective']:
    r={'vector':q['vector'],'best_objective':q['objective_value'],
       'seed':q['starting_seed'],'refined':True}
  fine=src.parent/'paper_step3_refine_11'/f'{objective}.json'
  if fine.exists():
   r={'vector':json.loads(fine.read_text())['final_vector']}
  ret[objective]=r['vector']
  if include_secondary and objective=='isolated_transfer':
   other=max(rs,key=lambda x:x['best_objective'])
   ret['isolated_second_seed']=other['vector']
 ret['historical_gaussian']=archived_start_vectors()[1][0].tolist()
 return ret

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=ROOT/'outputs/paper_step3')
 p.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step3_validation')
 p.add_argument('--steps',type=int,default=64);args=p.parse_args()
 args.output.mkdir(parents=True,exist_ok=True)
 filename=args.output/'nominal_11.json'
 if filename.exists():print('PRESERVED',filename,flush=True);return
 choices=candidate_dict(args.input)
 m=PairModel(frozen,VALID_BASIS);echo=echo_ref_matrix(m,steps=args.steps)
 rows=[]
 for kind,x in choices.items():
  met=measure_candidate(m,x,echo=echo,steps=args.steps)
  met.update(protocol=kind,design=x,basis_n=len(VALID_BASIS),echo_steps=args.steps)
  rows.append(met)
  print(kind,round(met['fidelity']*100,5),'aux',round(met['aux51']*100,4),
        'spectator',round(met['total_spectators']*100,4),'isolated',round(met['transfer_forward'],6),flush=True)
 target=filename.with_suffix('.pending.json');target.write_text(json.dumps(rows,indent=2)+'\n');target.replace(filename)
 print('SAVED',filename,flush=True)
if __name__=='__main__':main()
