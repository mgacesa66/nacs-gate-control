"""Independent large selected-basis nominal check, extending Step 2 21→23→25.

Do not use 25 selected states as evidence of entire 128-state-pair convergence.
"""
import argparse,csv,json
from pathlib import Path
from nacs.analysis.basis_convergence import model_gate,compare_to_reference,control_designs,combined_ranking
from nacs.analysis.parameter_recalibration import PROBES
from nacs.data import ROOT

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step2_extended');p.add_argument('--steps',type=int,default=64);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    rank,_=combined_ranking(PROBES)
    seed=[38,69,40,72]
    ordered=seed+[z['old_index'] for z in rank if z['old_index'] not in seed]
    select=['unmodulated','aux51','all_spectators']
    ds={k:v for k,v in control_designs().items() if k in select}
    r0,m0=model_gate(PROBES[0],'frequency_rabi',[],ds,steps=a.steps)
    baseline={r['protocol']:r for r in r0}
    rows=[]
    for n in (11,15,17,19,21,23,25):
        ids=[] if n==11 else ordered[:n-11]
        if n==11:rr,mm=r0,m0
        else:rr,mm=model_gate(PROBES[0],'frequency_rabi',ids,ds,steps=a.steps)
        for r in rr:
            info=compare_to_reference(baseline[r['protocol']],r,m0[r['protocol']],mm[r['protocol']])
            rows.append({**r,**info})
        print('nominal N',n,[(r['protocol'],round(r['fidelity'],8)) for r in rr],flush=True)
    with (a.output/'extended_nominal.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)
    (a.output/'basis_order.json').write_text(json.dumps({'ordering':ordered[:20],'steps':a.steps},indent=2)+'\n')
if __name__=='__main__':main()
