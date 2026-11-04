"""Step 2 systematic basis convergence, no plots, no optimizer changes.

Set OPENBLAS_NUM_THREADS=1. This run is intentionally staged: `--profile fast`
for CI and `--profile full` for the publication numerical benchmark.
"""
import argparse,csv,json,hashlib
from pathlib import Path
from nacs.analysis.basis_convergence import *
from nacs.data import ROOT
from nacs.controls.fixed_polarization import VALID_BASIS


def write_csv(p,rows):
    if rows:
        with p.open('w',newline='') as fd:
            w=csv.DictWriter(fd,rows[0].keys());w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step2')
    p.add_argument('--profile',choices=('fast','full'),default='full')
    p.add_argument('--steps',type=int,default=48)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    allprobes=PROBES if a.profile=='full' else (PROBES[0],PROBES[5],PROBES[10])
    rank,screen=combined_ranking(allprobes)
    # Keep fixed legacy rank15 state set as mandatory first extension: comparability
    # with v0.1 release checks. Then expand using a fresh independently derived rank.
    seed=[38,69,40,72]
    ordered=seed+[r['old_index'] for r in rank if r['old_index'] not in seed]
    # If the lowest-overlap eigenstate in the extended selected basis mixes with
    # a nonincluded original level, force the nearest mixing partner into the basis.
    partners=[]
    for probe in allprobes:
        permutation,W=assignments(probe)
        for i in seed+ordered[:12]:
            if W[i,permutation[i]]<.85:
                partners.append({'probe':probe.name,'old_index':i,
                    'best_overlap':float(W[i,permutation[i]]),
                    'second_old_index':int(np.argsort(W[:,permutation[i]])[-2]),
                    'second_overlap':float(np.sort(W[:,permutation[i]])[-2]),
                    'partner_from_old_vector':int(np.argsort(W[i,:])[-2])})
    sizes=(11,13,15,17,19,21) if a.profile=='full' else (11,15,19)
    # 13 follows *original* extra78/80; 15 follows historical ranked38/69/40/72.
    bases={11:[],13:[78,80],15:seed,
       17:ordered[:6],19:ordered[:8],21:ordered[:10]}
    probe_cases=[(PROBES[0],('frequency_rabi',),sizes,('unmodulated','gate_only','aux51','all_spectators'))]
    if a.profile=='full':
        for probe in PROBES[1:]:
            modes=('fixed','frequency','frequency_rabi') if probe.name.startswith('ellipse') else ('frequency_rabi',)
            levels=(11,15,19) if probe.name.startswith('ellipse') else (11,15)
            # Limit expensive 19-state sweeps to selected ellipticities ±0.1 and ±0.5
            if '0p25' in probe.name:levels=(11,15)
            probe_cases.append((probe,modes,levels,('unmodulated','aux51','all_spectators')))
    else:
        for probe in (PROBES[5],PROBES[10]):
            probe_cases.append((probe,('frequency_rabi',),(11,15,19),('unmodulated','aux51','all_spectators')))
    rows=[]; comparisons=[]
    for probe,modes,levels,protocols in probe_cases:
        for mode in modes:
            ref_rows,ref_M=model_gate(probe,mode,[],{k:v for k,v in control_designs().items() if k in protocols},steps=a.steps)
            ref={r['protocol']:r for r in ref_rows}
            for n in levels:
                rs,ms=(ref_rows,ref_M) if n==11 else model_gate(probe,mode,bases[n],{k:v for k,v in control_designs().items() if k in protocols},steps=a.steps)
                for r in rs:
                    rows.append(dict(r))
                    info=compare_to_reference(ref[r['protocol']],r,ref_M[r['protocol']],ms[r['protocol']])
                    comparisons.append({'probe':probe.name,'mode':mode,'basis_n':n,
                           'protocol':r['protocol'],**info})
                print(probe.name,mode,'N',n,'F',[(r['protocol'],round(r['fidelity'],6)) for r in rs],flush=True)
    write_csv(a.output/'ranked_candidates.csv',rank)
    write_csv(a.output/'screening_all_probes.csv',screen)
    write_csv(a.output/'gate_basis_results.csv',rows)
    write_csv(a.output/'operator_convergence.csv',comparisons)
    if partners:write_csv(a.output/'tracking_ambiguity.csv',partners)
    meta={'method':'global max-overlap tracked 128-state spectrum, actual driven 2-molecule coherent propagation',
       'profile':a.profile,'gaussian_echo_midpoint_steps':a.steps,'candidate_heuristic':'max g/max(abs(detuning),10kHz) over core 0,1,e and physical transfer and echo fields',
       'basis_order':{str(k):v for k,v in bases.items()},'num_gate_rows':len(rows),
       'probes':[p.name for p in allprobes],'calibration_caveat':'idealized bare carrier tracking, not experimental measured line center',
       'limits':['finite screened subspaces only, not full 128x128 pair propagation','unmodeled motion, loss and microwave bandwidth','no objective reoptimization on perturbed spectra', 'ranking heuristic not rigorous omitted-amplitude bound']}
    (a.output/'run_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('TOTAL',len(rows),'rows, worst changes', sorted([(abs(x['delta_fidelity_pp']),x['probe'],x['basis_n'],x['protocol']) for x in comparisons],reverse=True)[:8],flush=True)

if __name__=='__main__':main()
