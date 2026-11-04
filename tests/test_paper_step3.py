"""Release regression and physical-invariant tests for publication Step 3."""
import json
from pathlib import Path
import numpy as np
import pytest
from nacs.data import ROOT
from nacs.analysis.matched_objectives import (PairModel,frozen,SEARCH_BASIS,
  VALID_BASIS,Design,archived_start_vectors,population,isolated_transfers,
  single_phase_pulse,feasible,design_from_vector,PHASE_RAMP_US,gate_output,
  echo_ref_matrix,measures)
from nacs.controls.fixed_polarization import local_transfer,rotated_linear_polarization,target_rabi
from nacs.controls.spectator import ANGLE
from nacs.analysis.column_propagator import gate_columns

DATA=ROOT/'benchmarks/data/paper_step3'

def load_nominal():
    return {x['protocol']:x for x in json.loads((DATA/'nominal/nominal_11.json').read_text())}


def test_single_molecule_propagator_same_as_archived_square_control():
    model=PairModel(frozen,SEARCH_BASIS)
    d=Design(phase_forward=0.,phase_reverse=0.,reverse_phase=0.)
    pf,pr=isolated_transfers(model,d)
    ref=local_transfer(model,'1','e',target_rabi(model,ANGLE),d.transfer_ms,
        d.transfer_detuning_khz,rotated_linear_polarization(ANGLE))
    assert abs(pf-ref)<1e-10
    # Reciprocity of forward/reverse *probability* is not automatic in a
    # multilevel unitary: |U_e1|² need not equal |U_1e|².
    u=single_phase_pulse(model,d.transfer_ms,d.transfer_detuning_khz,0.,0.)
    assert abs(pr-abs(u[model.local['1'],model.local['e']])**2)<1e-12
    assert 0<=pf<=1 and 0<=pr<=1


def test_column_propagator_reproduces_frozen_matrix_method():
    model=PairModel(frozen,SEARCH_BASIS)
    _,history=archived_start_vectors()
    d=design_from_vector(history[3])
    a=gate_columns(model,d,steps=12,ramp_us=PHASE_RAMP_US,ramp_steps=6)
    b=gate_output(model,d,steps=12,phase_ramp_us=PHASE_RAMP_US,ramp_steps=6)
    assert np.max(abs(a-b))<2e-9
    assert np.max(abs(a.conj().T@a-b.conj().T@b))<2e-9


def test_four_arms_have_identical_starting_population_and_control_bounds():
    a=population(2026)
    b=population(2026)
    assert np.array_equal(a,b)
    assert a.shape==(30,6)
    from nacs.analysis.matched_objectives import BOUNDS
    bounds=np.array(BOUNDS)
    assert np.all((a>=bounds[:,0]) & (a<=bounds[:,1]))
    # Identical populations may include designs excluded by the shared 1.30-ms constraint.
    for arm in ('isolated_transfer','selected_aux51','channel_fidelity','channel_guarded'):
        for seed in (2026,2027):
            row=json.loads((DATA/f'global_checkpoints/{arm}_seed{seed}.json').read_text())
            assert row['evaluations']==240
            assert row['echo_steps']==16
            assert row['ramp_us']==2.
        fine=json.loads((DATA/f'refine_11/{arm}.json').read_text())
        assert fine['n_objective_evaluations']==181
        assert fine['refinement_basis_n']==11
        assert feasible(fine['final_vector'])


def test_nominal_complete_gate_regression_with_physical_partition():
    model=PairModel(frozen,VALID_BASIS)
    echo=echo_ref_matrix(model,steps=16)
    records=load_nominal()
    for name in ('isolated_transfer','selected_aux51','channel_guarded'):
        d=design_from_vector(records[name]['design'])
        vals=measures(model,gate_output(model,d,echo_ref=echo,steps=16,
                                  phase_ramp_us=PHASE_RAMP_US,ramp_steps=6))
        assert abs(vals['fidelity']-records[name]['fidelity'])<.0008
        assert abs(vals['total_spectators']-(vals['aux51']+vals['other_spectators']))<1e-10
        assert abs(vals['leakage']-(vals['total_spectators']+vals['intended_nonlogical']))<1e-10
        assert 0<=vals['fidelity']<=1


def test_25_state_screened_validation_keeps_nominal_ranking():
    def get(n):return json.loads((DATA/f'basis/{n}_25_columns_steps24.json').read_text())
    a,b,c=get('isolated_transfer'),get('selected_aux51'),get('channel_guarded')
    assert a['pair_dimension']==b['pair_dimension']==c['pair_dimension']==625
    assert c['fidelity']>b['fidelity']>a['fidelity']
    assert max(abs(get(n)['fidelity']-load_nominal()[n]['fidelity']) for n in
       ('isolated_transfer','selected_aux51','channel_guarded'))<.0003
    for name in ('isolated_transfer','isolated_second_seed'):
        assert get(name)['dense_expm_max_complex_diff']<3e-10


def test_noise_coordinates_paired_across_all_arms():
    import csv
    names=('isolated_transfer','selected_aux51','channel_fidelity',
           'channel_guarded','historical_gaussian')
    for scenario in ('tight','broad'):
        records=[]
        for name in names:
            with (DATA/f'heldout/{name}_{scenario}.csv').open() as fp:
                rows=list(csv.DictReader(fp))
            assert len(rows)==24
            err=np.array([[float(r[k]) for k in ('amp_error','det_error_khz','exchange_error')]
                           for r in rows])
            assert np.isfinite([float(r['fidelity']) for r in rows]).all()
            records.append(err)
        assert all(np.array_equal(records[0],r) for r in records[1:])


def test_11_to_25_state_complex_operator_agreement_after_allowed_phase_alignment():
    rows=json.loads((DATA/'basis/operator_alignment_11_to_25.json').read_text())
    assert len(rows)==6
    assert all(0<=r['aligned_frobenius']<.002 for r in rows)
    assert all(r['raw_frobenius']>=r['aligned_frobenius'] for r in rows)
    with np.load(DATA/'basis/logical_matrices_11_25.npz') as ar:
        assert len(ar.files)==12
        assert all(ar[k].shape==(4,4) for k in ar.files)
