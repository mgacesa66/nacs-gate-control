"""Spectroscopic predictions / comparisons, with explicit selection conventions.

All known benchmark frequencies are taken from Picard et al., arXiv:2406.15345
Table I (v1), relative to 2 B_rot. Relative Rabi strengths are absolute
magnitudes, but experimental polarization amplitudes may differ by line.

IMPORTANT: This is an *unfitted* DFT-parameter model; frequency errors
should NOT be silently corrected by shifting individual states.
"""
import csv
import numpy as np

# NOTE that 189(4) printed for |1>-|e> is inconsistent as a plain ratio;
# do not use that entry as a measured dipole amplitude.
MEASURED = [
    {'name':'0->e','frequency_offset_khz':+0.2,'relative_strength':1.0},
    {'name':'0->f_z','frequency_offset_khz':+691.7,'relative_strength':0.77},
    {'name':'0->f_y','frequency_offset_khz':-710.3,'relative_strength':0.84},
    {'name':'1->e','frequency_offset_khz':-1902.2,'relative_strength':None},
]


def polarization_xyz_coefficients(axis='x'):
    """Spherical d_q coefficients in epsilon . d with epsilon along x,y,z.

    e_x . d = (d_-1 - d_+1)/sqrt(2)
    e_y . d = i(d_-1 + d_+1)/sqrt(2)
    e_z . d = d_0
    """
    if axis=='x':return {-1:1/np.sqrt(2),0:0,+1:-1/np.sqrt(2)}
    if axis=='y':return {-1:1j/np.sqrt(2),0:0,+1:1j/np.sqrt(2)}
    if axis=='z':return {-1:0,0:1,+1:0}
    raise ValueError(axis)


def benchmark_lines(model,sp=None,states=None):
    """Return comparison records, not an assertion of experimental agreement.

    Experimental strengths are based on their microwave polarizations,
    which the paper does not specify here as absolute Cartesian calibration.
    Comparing all with x-polarization is a diagnostic, NOT apples-to-apples
    validation of experimental relative Rabi amplitudes.
    """
    if sp is None:sp=model.diagonalize()
    if states is None:states=model.identify_gate_states(sp)
    choices={k:v['index'] for k,v in states.items()}
    ref=abs(sp.dressed_dipole(choices['0'],choices['e'],polarization_xyz_coefficients('x')))
    result=[]
    for line in MEASURED:
        a,b=line['name'].split('->')
        if a not in choices or b not in choices:continue
        ia,ib=choices[a],choices[b]
        freq=sp.energies_khz[ib]-sp.energies_khz[ia]-2*model.p.B_rot_khz
        dx=sp.dressed_dipole(ia,ib,polarization_xyz_coefficients('x'))
        dy=sp.dressed_dipole(ia,ib,polarization_xyz_coefficients('y'))
        dz=sp.dressed_dipole(ia,ib,polarization_xyz_coefficients('z'))
        preferred_axis={'0->e':'x','0->f_z':'z','0->f_y':'y','1->e':'x'}[line['name']]
        dp=sp.dressed_dipole(ia,ib,polarization_xyz_coefficients(preferred_axis))
        res={**line,'comparison_polarization':preferred_axis,'predicted_offset_khz':float(freq),
             'residual_khz':float(freq-line['frequency_offset_khz']),
             'predicted_relative_strength_for_line':float(abs(dp)/ref) if ref else float('nan'),
             'predicted_max_xyz_strength':float(max(abs(dx),abs(dy),abs(dz))/ref) if ref else float('nan'),
             'd_x_real_debye':float(np.real(dx)),'d_x_imag_debye':float(np.imag(dx)),
             'i':int(ia),'j':int(ib)}
        result.append(res)
    return result


def characterise_spectators(model,sp=None,states=None,top=10):
    """All near-by microwave-coupled N=1 dressed eigenstates from |1>.

    Rabi amplitudes are for an *assumed x-polarized microwave* (not
    necessarily the experimental horn polarization). Detunings are relative
    to the dressed |1>->|e> transition; no microscopic line is invented.
    """
    if sp is None:sp=model.diagonalize()
    if states is None:states=model.identify_gate_states(sp)
    source=states['1']['index']; target=states['e']['index']
    e1=sp.energies_khz[target]-sp.energies_khz[source]
    px=polarization_xyz_coefficients('x')
    target_d=sp.dressed_dipole(source,target,px)
    data=[]
    for j,E in enumerate(sp.energies_khz):
        if j==target:continue
        pN=sum(abs(sp.vectors[k,j])**2 for k,x in enumerate(sp.basis) if x[0]==1)
        if pN<0.5:continue
        d=sp.dressed_dipole(source,j,px)
        if abs(d)<1e-9: continue
        data.append({'state':int(j),'relative_detuning_khz':float(E-sp.energies_khz[source]-e1),
                     'dipole_x_real_debye':float(np.real(d)),'dipole_x_imag_debye':float(np.imag(d)),
                     'dipole_abs_debye':float(abs(d)),
                     'relative_to_weak_target':float(abs(d)/abs(target_d)) if abs(target_d)>0 else float('inf'),
                     'dominant_bare':str(sp.dominant_components(j,1)[0][0])})
    data.sort(key=lambda x:abs(x['relative_detuning_khz']))
    return data[:top]


def write_csv(path,records):
    if not records:return
    with open(path,'w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0].keys()))
        writer.writeheader();writer.writerows(records)
