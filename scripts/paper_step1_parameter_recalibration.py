"""Publication-plan Step 1 benchmark. No fitted parameter ranges or figures.

Run from repo root:
 OPENBLAS_NUM_THREADS=1 python scripts/paper_step1_parameter_recalibration.py \
      --output outputs/paper_step1
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,csv,json,argparse,datetime,hashlib
from pathlib import Path
from dataclasses import asdict
import numpy as np
from nacs.data import benchmark_path,ROOT
from nacs.controls.spectator import Design
from nacs.analysis.parameter_recalibration import (PROBES,MODES,nominal_objects,tracked_spectrum,
   drive_parameters,build_echo,full_gate_output,measures,PairModel,VALID_BASIS,evaluate)


def write_csv(path,records):
    if not records:return
    with path.open('w',newline='') as f:
        wr=csv.DictWriter(f,records[0].keys());wr.writeheader();wr.writerows(records)


def main():
    a=argparse.ArgumentParser()
    a.add_argument('--output',type=Path,default=ROOT/'outputs/paper_step1')
    a.add_argument('--steps',type=int,default=64)
    a.add_argument('--only-probe',choices=[p.name for p in PROBES])
    opt=a.parse_args();opt.output.mkdir(parents=True,exist_ok=True)
    archived=json.loads(benchmark_path('optimization/candidates.json').read_text())
    designs={k:Design(**v['design']) for k,v in archived.items()}
    # Nominal Z angles defined once for each original stored pulse program.
    basef,_=tracked_spectrum(PROBES[0]); base=PairModel(basef,VALID_BASIS)
    bp=drive_parameters(basef,'fixed',next(iter(designs.values())))
    echo=build_echo(base,bp,steps=opt.steps)
    reference_z={}
    for name,design in designs.items():
        q=measures(base,full_gate_output(base,design,drive_parameters(basef,'fixed',design),echo))
        reference_z[name]=q['z'].tolist()
    allrows=[];spectral=[];summary=[]
    for probe in PROBES:
        if opt.only_probe and opt.only_probe!=probe.name:continue
        f,diagnostics=tracked_spectrum(probe)
        d=diagnostics
        old,_,_=nominal_objects()
        p=drive_parameters(f,'fixed',next(iter(designs.values())))
        spectral.append({'probe':probe.name,**asdict(probe),
           'transfer_shift_khz':p['transfer_line_shift_khz'],
           'echo_shift_khz':p['echo_line_shift_khz'],
           'transfer_dipole_ratio':p['transfer_dipole_ratio'],
           'echo_dipole_ratio':p['echo_dipole_ratio'],
           'min_overlap_primary':d['min_tracking_primary'],
           'min_overlap_all_selected':d['min_tracking'],
           'min_overlap_ranked15':d['min_tracking_ranked15'],
           'overlap_rank38':d['tracking_ranked']['38'],
           'overlap_rank69':d['tracking_ranked']['69'],
           'overlap_1':d['tracking']['1'],'overlap_e':d['tracking']['e'],
           'overlap_aux51':d['tracking']['aux51'],
           'aux51_gap_from_e_khz':float((f.energies[f.indices['aux51']]-f.energies[f.indices['e']])-
                  (old.energies[old.indices['aux51']]-old.energies[old.indices['e']]))})
        for mode in MODES:
            rows,info=evaluate(probe,mode,designs,steps=opt.steps,reference_z=reference_z)
            allrows.extend(rows)
            x={z['protocol']:z for z in rows}
            summary.append({'probe':probe.name,'mode':mode,
               'fidelity_unmodulated':x['unmodulated']['fidelity_phase_refit'],
               'fidelity_aux51_target':x['aux51']['fidelity_phase_refit'],
               'fidelity_all_guarded':x['all_spectators']['fidelity_phase_refit'],
               'target_minus_unmodulated':x['aux51']['fidelity_phase_refit']-x['unmodulated']['fidelity_phase_refit'],
               'guarded_minus_unmodulated':x['all_spectators']['fidelity_phase_refit']-x['unmodulated']['fidelity_phase_refit'],
               'aux51_target_reduction_fraction':1-x['aux51']['aux51']/max(x['unmodulated']['aux51'],1e-12),
               'aux51_other_spectator_growth':x['aux51']['other_spectators']-x['unmodulated']['other_spectators'],
               'guarded_leakage':x['all_spectators']['leakage']})
        print(f'{probe.name:25s}', 'shifts (kHz) =',round(p['transfer_line_shift_khz'],2),round(p['echo_line_shift_khz'],2),
                 'fidelity guarded fixed/retuned =',*[round(next(x for x in allrows if x['probe']==probe.name and x['mode']==mode and x['protocol']=='all_spectators')['fidelity_phase_refit'],4) for mode in MODES],flush=True)
    write_csv(opt.output/'gate_comparison.csv',allrows)
    write_csv(opt.output/'spectral_tracking.csv',spectral)
    write_csv(opt.output/'hypothesis_summary.csv',summary)
    meta={'title':'arXiv step1 molecular parameter sensitivity and recalibration',
      'frozen_release_sha256':hashlib.sha256((ROOT.parent/'NaCs_GitHub_technical_report_code_v0_1.zip').read_bytes()).hexdigest() if (ROOT.parent/'NaCs_GitHub_technical_report_code_v0_1.zip').exists() else 'see provenance',
      'frozen_spectroscopy_checksum':nominal_objects()[0].checksum,
      'publication_plan':'PUBLICATION_PLAN_FROZEN.md','calibration_modes':list(MODES),
      'probe_hypothetical':True, 'echo_gaussian_midpoint_steps':opt.steps,
      'phase_calibration':'both fixed baseline virtual-Z and optimistic locally refitted Z reported',
      'controls':'stored locally optimized protocols; no scenario-by-scenario gate reoptimization',
      'spectral_labels':'global maximum-overlap assignment, phase gauge aligned',
      'caveats':['power-dependent AC Stark line-center corrections not calibrated','eigenstate near crossings can be ambiguous',
                 'no motional/loss/noise model in this audit','carrier changes idealized bare dressed transition frequencies']}
    (opt.output/'run_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('completed',len(spectral),'probes',len(allrows),'gate records')
if __name__=='__main__':main()
