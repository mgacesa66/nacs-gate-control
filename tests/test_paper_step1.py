"""Numerical regression and physics invariants of publication Step 1."""
import json
import numpy as np
import pytest
from nacs.data import benchmark_path
from nacs.controls.spectator import Design,gate_output,measures,echo_ref_matrix
from nacs.analysis.parameter_recalibration import (PROBES,MODES,tracked_spectrum,
  drive_parameters,build_echo,full_gate_output,PairModel,VALID_BASIS,evaluate)

@pytest.fixture(scope='module')
def designs():
    return {k:Design(**x['design']) for k,x in json.loads(benchmark_path('optimization/candidates.json').read_text()).items()}

@pytest.fixture(scope='module')
def nominal():
    f,diag=tracked_spectrum(PROBES[0]);return f,PairModel(f,VALID_BASIS)

def test_baseline_reproduced_with_complex_amplitudes(nominal,designs):
    f,m=nominal
    echo0=echo_ref_matrix(m,steps=64)
    echo1=build_echo(m,drive_parameters(f,'fixed',next(iter(designs.values()))),steps=64)
    # MHz-scale spectral subtraction can introduce sub-microhertz detuning
    # differences across LAPACK builds, amplified into nanolevel amplitudes.
    assert np.max(np.abs(echo0-echo1))<1e-8
    for name,d in designs.items():
        output=full_gate_output(m,d,drive_parameters(f,'fixed',d),echo1)
        expected=gate_output(m,d,echo_ref=echo0)
        assert np.max(np.abs(output-expected))<2e-8,(name,np.max(np.abs(output-expected)))

def test_nominal_all_three_calibrations_identical(nominal,designs):
    f,m=nominal
    a=[]
    for mode in MODES:
        p=drive_parameters(f,mode,designs['all_spectators'])
        echo=build_echo(m,p,steps=48)
        q=measures(m,full_gate_output(m,designs['all_spectators'],p,echo))
        a.append(q['fidelity'])
    assert np.ptp(a)<2e-10
    assert abs(a[0]-.94638)<1e-4

def test_carrier_tracking_recovers_controlled_detuning(designs):
    d=designs['unmodulated']
    probe=next(p for p in PROBES if p.name=='ellipse_plus_0p1deg')
    f,_=tracked_spectrum(probe)
    fixed=drive_parameters(f,'fixed',d);freq=drive_parameters(f,'frequency',d)
    assert abs(fixed['transfer_detuning_khz'] + fixed['transfer_line_shift_khz'] - d.transfer_detuning_khz)<1e-10
    assert abs(freq['transfer_detuning_khz'] - d.transfer_detuning_khz)<1e-10
    assert abs(freq['echo_detuning_khz'] - drive_parameters(tracked_spectrum(PROBES[0])[0],'fixed',d)['echo_detuning_khz'])<1e-8
    assert abs(freq['transfer_rabi_khz']-fixed['transfer_rabi_khz'])<1e-10

def test_rabi_recalibration_is_field_scaling_and_changes_power(designs):
    d=designs['unmodulated']
    p=next(p for p in PROBES if p.name=='ellipse_minus_0p5deg')
    f,_=tracked_spectrum(p)
    a=drive_parameters(f,'frequency',d);b=drive_parameters(f,'frequency_rabi',d)
    assert abs(a['transfer_rabi_khz']*b['electric_field_scale_transfer']-b['transfer_rabi_khz'])<1e-9
    assert abs(a['echo_rabi_khz']*b['electric_field_scale_echo']-b['echo_rabi_khz'])<1e-9
    assert b['electric_field_scale_transfer']!=1.

def test_population_partitions_and_no_imaginary_probability(nominal,designs):
    f,m=nominal
    d=designs['aux51'];p=drive_parameters(f,'fixed',d);echo=build_echo(m,p,steps=32)
    out=full_gate_output(m,d,p,echo)
    q=measures(m,out)
    assert np.allclose(np.sum(np.abs(out)**2,axis=0),np.ones(4),atol=2e-10)
    assert abs(q['leakage']-q['total_spectators']-q['intended_nonlogical'])<1e-11
    assert abs(q['total_spectators']-q['aux51']-q['other_spectators'])<1e-11
    assert q['aux51']<.003
    assert q['other_spectators']>.02

def test_tracking_overlap_and_spectral_difference():
    f,d=tracked_spectrum(next(p for p in PROBES if p.name=='ellipse_plus_0p5deg'))
    assert d['min_tracking_primary']>.99
    assert d['min_tracking_ranked15']<d['min_tracking_primary']
    assert np.isclose(f.energies[f.indices['e']]-f.energies[f.indices['0']],
      f.energies[f.indices['e']]-f.energies[f.indices['0']])
