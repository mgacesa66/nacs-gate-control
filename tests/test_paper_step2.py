import numpy as np
from nacs.analysis.basis_convergence import (screen_one,combined_ranking,build_basis,model_gate,
    compare_to_reference,assignments)
from nacs.analysis.parameter_recalibration import PROBES,tracked_spectrum
from nacs.controls.fixed_polarization import VALID_BASIS

def test_screen_full_labels_and_ranking():
    rank=screen_one(PROBES[0]);assert len(rank)==128-len(VALID_BASIS)
    assert all(row['score']>=0 for row in rank)
    assert len({r['old_index'] for r in rank})==len(rank)
    assert all(rank[k]['score']>=rank[k+1]['score'] for k in range(len(rank)-1))

def test_tracking_permutation_is_one_to_one():
    p,_=assignments(PROBES[6]);assert len(set(p))==128
    f,_=tracked_spectrum(PROBES[6]);basis=build_basis(f,PROBES[6],[38,69,40,72])
    assert len({f.indices[k] for k in basis})==len(basis)

def test_baseline_gate_fidelity_from_new_engine():
    r,m=model_gate(PROBES[0],'frequency_rabi',[],steps=32)
    vals={x['protocol']:x for x in r}
    assert np.isclose(vals['unmodulated']['fidelity'],.9375245330782948,atol=3e-5)
    assert np.isclose(vals['all_spectators']['fidelity'],.9463802424458866,atol=3e-5)
    assert all(x['unitarity_error']<1e-8 for x in r)
    assert all(v.shape==(4,4) for v in m.values())

def test_13_15_state_gates_are_physical():
    base,ref=model_gate(PROBES[0],'frequency_rabi',[],steps=24)
    rows,mat=model_gate(PROBES[0],'frequency_rabi',[38,69,40,72],steps=24)
    assert all(0<=x['logical_leakage']<1 for x in rows)
    diff=compare_to_reference(base[0],rows[0],ref['unmodulated'],mat['unmodulated'])
    assert abs(diff['delta_fidelity_pp']) < .5
    assert diff['logical_operator_frobenius_diff'] < .05

def test_independent_dipolar_screen_has_missing_near_resonances():
    from nacs.analysis.basis_convergence import screen_dipolar_one
    row=screen_dipolar_one(PROBES[0])
    assert len(row)==128-len(VALID_BASIS)
    assert {34,2}.issubset({x['old_index'] for x in row[:5]})
    assert all(row[i]['score_dd']>=row[i+1]['score_dd'] for i in range(len(row)-1))
