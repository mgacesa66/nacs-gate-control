"""All seven Picard et al. 2024 preprint Extended Data Table I frequencies.

Separate from legacy 4-line validator and existing gate adapter.
Energies from 128-state Hamiltonian; full frequency residuals use the
paper's 'frequency - 2B' convention. Three derived N=0 splittings are
frequency-reference-independent and provided by ground_independent.py.
"""
from __future__ import annotations
import numpy as np
from .validation import polarization_xyz_coefficients

# Printed numbers from arXiv:2406.15345v1, Extended Data Table I, p. 14.
# Frequency uncertainty convention 691.7(2)-> 0.2 kHz.
# The last Rabi ratio is printed as 189(4) (not reliably a normalized ratio),
# hence used as raw text, not a dimensionless target.
LINES=(
    ('0->f_z','0','f_z',691.7,.2,.77,'z'),
    ('0->e','0','e',.2,.1,1.,'x'),
    ('0->f_y','0','f_y',-710.3,.2,.84,'y'),
    ('0->f_weak','0','f_weak',474.6,.2,.01838,'z'),
    ('g_plus->e','g_plus','e',466.1,.3,.0053,'z'),
    ('g_minus->e','g_minus','e',-927.1,.3,.0218,'y'),
    ('1->e','1','e',-1902.2,.1,None,'x'),
)


def identify_seven_states(model,sp):
    states=model.identify_gate_states(sp).copy()
    weak=model.superposition([(-1/np.sqrt(2),1,1,1.5,1.5),
                              (+1/np.sqrt(2),1,-1,1.5,1.5)])
    for name,bare,rot in [('f_weak',weak,1),
                 ('g_plus',model.bare_state(0,0,1.5,3.5),0),
                 ('g_minus',model.bare_state(0,0,1.5,.5),0)]:
        k,w=sp.best_overlap(bare,allowed_N=rot)
        states[name]={'index':int(k),'bare_overlap':float(w),
                      'energy_khz':float(sp.energies_khz[k])}
    return states


def benchmark_seven(model,sp=None,states=None):
    if sp is None:sp=model.diagonalize()
    if states is None:states=identify_seven_states(model,sp)
    ix={label:rec['index'] for label,rec in states.items()}
    ref=sp.dressed_dipole(ix['0'],ix['e'],polarization_xyz_coefficients('x'))
    out=[]
    for name,a,b,f,sigma,strength,pol in LINES:
        ia,ib=ix[a],ix[b]
        theoretical=float(sp.energies_khz[ib]-sp.energies_khz[ia]-2*model.p.B_rot_khz)
        d=sp.dressed_dipole(ia,ib,polarization_xyz_coefficients(pol))
        out.append({'name':name,'measured_offset_khz':f,'sigma_khz':sigma,
                    'predicted_offset_khz':theoretical,'residual_khz':theoretical-f,
                    'measured_relative_rabi':strength,
                    'measured_rabi_raw':('189(4), unresolved units' if name=='1->e' else ''),
                    'assumed_mw_axis':pol,
                    'predicted_relative_rabi':float(abs(d)/abs(ref)),
                    'max_xyz_relative_rabi':float(max(abs(sp.dressed_dipole(ia,ib,polarization_xyz_coefficients(q))) for q in 'xyz')/abs(ref)),
                    'bare_state_overlap_initial':states[a]['bare_overlap'],
                    'bare_state_overlap_final':states[b]['bare_overlap'],
                    'd_real_debye':float(d.real),'d_imag_debye':float(d.imag),
                    'state_i':ia,'state_f':ib})
    return out


def centered_differences(lines):
    """Reference-free differences between lines with common initial state 0."""
    d={r['name']:r for r in lines}
    ref=d['0->e']
    result=[]
    for label in ('0->f_z','0->f_y','0->f_weak'):
        a=d[label]
        obs=a['measured_offset_khz']-ref['measured_offset_khz']
        pred=a['predicted_offset_khz']-ref['predicted_offset_khz']
        result.append({'line':label+' minus 0->e',
                       'measured_khz':obs,'predicted_khz':pred,'residual_khz':pred-obs})
    return result
