"""Summarize archived held-out errors, with *conditional* paired bootstrap CIs.

These intervals quantify sampling variability in finite simulated ensembles only.
They do NOT capture uncertainty in the molecular Hamiltonian, hardware or optimizers.
"""
from pathlib import Path
import argparse,csv,json
import numpy as np
from nacs.data import benchmark_path

def summarize(nboot=10000,seed=1101):
    rows=list(csv.DictReader(benchmark_path('optimization/heldout_noise.csv').open()))
    arr={}
    for r in rows:
        arr[(r['scenario'],int(r['sample']),r['protocol'])]=float(r['fidelity'])
    rng=np.random.default_rng(seed)
    output=[]
    for scenario in ('tight','broad'):
        for a,b in [('gate_only','all_spectators'),('unmodulated','aux51'),('unmodulated','all_spectators')]:
            indices=sorted({i for tag,i,protocol in arr if tag==scenario and protocol==a})
            dif=np.array([arr[scenario,i,b]-arr[scenario,i,a] for i in indices])
            resamp=rng.integers(0,len(dif),size=(nboot,len(dif)))
            lower,upper=np.percentile(dif[resamp].mean(axis=1),[2.5,97.5])
            output.append({'scenario':scenario,'protocol_A':a,'protocol_B':b,'n_simulated_pairs':len(dif),
              'paired_mean_B_minus_A':float(dif.mean()),'bootstrap_95_low':float(lower),
              'bootstrap_95_high':float(upper),
              'interpretation':'conditional Monte Carlo sampling only; not experimental/model uncertainty'})
    return output

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path('outputs/uncertainty_summary.json'))
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    data=summarize();a.output.write_text(json.dumps(data,indent=2)+'\n')
    for r in data:
        print(r['scenario'],r['protocol_A'],r['protocol_B'],
           f"delta={r['paired_mean_B_minus_A']:+.5f}",
           f"CI [{r['bootstrap_95_low']:+.5f},{r['bootstrap_95_high']:+.5f}]")
