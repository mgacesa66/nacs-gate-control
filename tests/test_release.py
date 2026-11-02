"""Tests exercise actual propagation, not only metadata or file existence."""
import csv,json,hashlib
from pathlib import Path
import numpy as np
import pytest
from nacs.data import ROOT, frozen_spectrum, benchmark_path
from nacs.gate.model import PairModel,BASIS_SETS
from nacs.gate.propagate import Controls
from nacs.controls.fixed_polarization import VALID_BASIS,EXTRA_BASIS,rotated_linear_polarization
from nacs.controls.gaussian_echo import gate_matrix
from nacs.analysis.stage_leakage import archived_controls,trace
from nacs.controls.spectator import Design,gate_output,measures,state_breakdown,echo_ref_matrix

@pytest.fixture(scope='module')
def model():return PairModel(frozen_spectrum(),VALID_BASIS)

@pytest.fixture(scope='module')
def candidate_data():return json.loads(benchmark_path('optimization/candidates.json').read_text())

def test_frozen_checksum_and_source_labels():
    f=frozen_spectrum(); assert len(f.energies)==128
    assert hashlib.sha256((ROOT/'benchmarks/frozen/gate_ready.npz').read_bytes()).hexdigest()==f.checksum
    assert len(f.dipoles)==3
    assert all(np.allclose(d,d.conj().T,atol=1e-9) for d in f.dipoles)
    assert f.indices['1']==11 and f.indices['e']==37

def test_canonical_basis_contains_known_omission(model):
    assert 'aux51' in VALID_BASIS
    assert 'aux51' not in BASIS_SETS['nine']
    assert model.size==11**2
    assert 'aux51' in model.local

def test_hermiticity_and_unitarity(model):
    assert np.allclose(model.Hfree,model.Hfree.conj().T,atol=1e-9)
    from scipy.linalg import expm
    u=expm(-2j*np.pi*model.Hfree*.02)
    assert np.allclose(u.conj().T@u,np.eye(model.size),atol=1e-9)

def test_gaussian_regression(model):
    c=archived_controls('gaussian');M=gate_matrix(model,c,35,steps=64)
    from nacs.controls.spectator import Design
    mm=measures(model,gate_output(model,Design()))
    # Full protocol must match frozen reference to numerical precision.
    assert abs(mm['fidelity']-.9375245330782948)<1e-8
    assert M.shape==(4,4)

def test_population_partition_and_state_specific_objectives(model,candidate_data):
    echo=echo_ref_matrix(model,steps=32)
    vals={}
    for name in ['unmodulated','aux51','all_spectators']:
        d=Design(**candidate_data[name]['design'])
        q=measures(model,gate_output(model,d,echo_ref=echo))
        assert np.isclose(q['total_spectators'],q['aux51']+q['other_spectators'],atol=1e-11)
        assert np.isclose(q['leakage'],q['total_spectators']+q['intended_nonlogical'],atol=1e-11)
        vals[name]=q
    assert vals['aux51']['aux51'] < vals['unmodulated']['aux51']*.2
    assert vals['aux51']['other_spectators'] > vals['unmodulated']['other_spectators']
    assert vals['aux51']['fidelity'] < vals['unmodulated']['fidelity']
    assert vals['all_spectators']['fidelity'] > vals['unmodulated']['fidelity']

def test_stability_data_provenance():
    rows=list(csv.DictReader(benchmark_path('stability/map_49_cells.csv').open()))
    assert len(rows)==49
    assert len(list(csv.DictReader(benchmark_path('spectroscopy/seven_line_accounting.csv').open())))==7

def test_source_archive_manifest_not_misrepresented():
    p=ROOT/'frozen/manifests/source_archives_sha256.json'
    entries=json.loads(p.read_text()); assert len(entries)==7
    assert all(len(v)==64 for v in entries.values())

@pytest.mark.slow
def test_ranked_basis_convergence(candidate_data):
    f=frozen_spectrum()
    # New objects for tests, no mutations to frozen npz.
    for k,i in {'rank38':38,'rank69':69,'rank40':40,'rank72':72}.items():f.indices[k]=i
    designs=['unmodulated','aux51','all_spectators']
    scores={}
    for key,basis in [('eleven',VALID_BASIS),('thirteen',EXTRA_BASIS),
                      ('ranked_fifteen',VALID_BASIS+('rank38','rank69','rank40','rank72'))]:
        m=PairModel(f,basis);e=echo_ref_matrix(m,steps=48)
        scores[key]={}
        for protocol in designs:
            d=Design(**candidate_data[protocol]['design'])
            scores[key][protocol]=measures(m,gate_output(m,d,echo_ref=e))['fidelity']
    for name in designs:
        assert abs(scores['ranked_fifteen'][name]-scores['eleven'][name])<0.01


def test_microscopic_hamiltonian_and_dipole_selection():
    from nacs.spectroscopy.core import NaCsMolecule
    from nacs.spectroscopy.constants import NaCsParameters
    from nacs.controls.gaussian_echo import shape_area
    mol=NaCsMolecule(NaCsParameters(optical_tensor_khz=-3526.155))
    H=mol.hamiltonian_khz()
    assert H.shape==(128,128)
    assert np.allclose(H,H.conj().T,atol=1e-9)
    D=mol.bare_dipole_q_debye()
    assert all(v.shape==(128,128) for v in D.values())
    assert abs(shape_area('truncated')-.6801128911235)<1e-9
