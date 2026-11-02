"""Independent 32-state NaCs N=0 nuclear-spin Hamiltonian.

No imports from core.py, no rotation/tweezer states, no Wigner functions.
The 32 by 32 matrix is explicitly assembled by the matrix elements of
H/h = -nu_N B (g_Na*(1-sigma_Na)*I_Na,z + g_Cs*(1-sigma_Cs)*I_Cs,z)
      + c4 * I_Na . I_Cs.
The off-diagonals are generated from ladder-factor formulae and the
matrix is diagonalized block-by-block in conserved M=m_Na+m_Cs.
All energies are cyclic kHz; eigenvectors are in product-state basis.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .constants import NaCsParameters, MU_N_OVER_H_KHZ_PER_G

I_NA = 1.5
I_CS = 3.5


@dataclass
class GroundSpectrum:
    basis: tuple[tuple[float,float], ...]
    H_khz: np.ndarray
    energies_khz: np.ndarray
    vectors: np.ndarray
    M_values: np.ndarray

    def eigenstate_for_bare(self, mna, mcs):
        row=self.basis.index((float(mna),float(mcs)))
        weights=np.abs(self.vectors[row,:])**2
        index=int(np.argmax(weights))
        return index,float(weights[index])

    def energy_for_bare(self, mna, mcs):
        k,overlap=self.eigenstate_for_bare(mna,mcs)
        return float(self.energies_khz[k]),overlap


def ground_basis():
    return tuple((float(mna),float(mcs))
                 for mna in np.arange(I_NA,-I_NA-1,-1)
                 for mcs in np.arange(I_CS,-I_CS-1,-1))


def _raise(j,m):
    return float(np.sqrt(max(0.,(j-m)*(j+m+1.))))


def _lower(j,m):
    return float(np.sqrt(max(0.,(j+m)*(j-m+1.))))


def ground_hamiltonian(params:NaCsParameters=NaCsParameters()):
    """Full 32x32 independent N=0 H/h in kHz; no arbitrary constant."""
    basis=ground_basis()
    index={state:k for k,state in enumerate(basis)}
    H=np.zeros((32,32),float)
    zn=MU_N_OVER_H_KHZ_PER_G * params.field_G
    for i,(a,b) in enumerate(basis):
        H[i,i] = -zn*(params.g_na*(1-params.shielding_na)*a
                      +params.g_cs*(1-params.shielding_cs)*b) + params.c4_khz*a*b
        new=(a+1.,b-1.)
        if new in index:
            j=index[new]
            term=.5*params.c4_khz*_raise(I_NA,a)*_lower(I_CS,b)
            H[i,j]+=term
            H[j,i]+=term
    assert np.allclose(H,H.T)
    return basis,H


def diagonalize_ground(params:NaCsParameters=NaCsParameters()):
    basis,H=ground_hamiltonian(params)
    energies=[]
    vectors=[]
    M_values=[]
    for M in sorted(set(a+b for a,b in basis)):
        inds=[i for i,(a,b) in enumerate(basis) if a+b==M]
        values,vec=np.linalg.eigh(H[np.ix_(inds,inds)])
        for k,value in enumerate(values):
            full=np.zeros(32)
            full[inds]=vec[:,k]
            pivot=int(np.argmax(abs(full)))
            if full[pivot]<0:full*=-1
            energies.append(float(value));vectors.append(full);M_values.append(M)
    order=np.argsort(energies)
    return GroundSpectrum(basis,H,np.asarray(energies)[order],
                          np.column_stack(vectors)[:,order],np.asarray(M_values)[order])


def compare_to_existing_core(params:NaCsParameters=NaCsParameters()):
    """Cross-check completely independently constructed H0 against 128-state core.

    The core optical scalar acts as identity and is removed from the core
    N=0 block; all rank-2 trap contributions vanish within N=0.
    """
    from .core import NaCsMolecule
    g=diagonalize_ground(params)
    core=NaCsMolecule(params)
    H_full=core.hamiltonian_khz()
    assert tuple((a[2],a[3]) for a in core.basis[:32])==g.basis
    H_core=H_full[:32,:32]-params.optical_scalar_khz*np.eye(32)
    dif=np.max(np.abs(H_core-g.H_khz))
    nu=np.linalg.eigvalsh(H_core)
    spectral=np.max(np.abs(nu-g.energies_khz))
    cross=np.max(np.abs(H_full[:32,32:]))
    return {'matrix_max_abs_difference_khz':float(dif),
            'eigenvalue_max_abs_difference_khz':float(spectral),
            'n0_to_n1_block_max_abs_khz':float(cross),
            'pass':bool(dif<1e-8 and spectral<1e-8 and cross<1e-8)}


# Directly extracted N=0 splitting observables from Picard Table I.
# E_ground - E_0 = nu_(0->e) - nu_(ground->e).
GROUND_OBSERVED={
    '1': {'bare':(-.5,2.5),'splitting_khz':1902.4, 'sigma_khz':np.hypot(.1,.1)},
    'g_plus': {'bare':(1.5,3.5),'splitting_khz':-465.9, 'sigma_khz':np.hypot(.1,.3)},
    'g_minus': {'bare':(1.5,.5),'splitting_khz':927.3, 'sigma_khz':np.hypot(.1,.3)},
}


def ground_splitting_comparison(params=NaCsParameters()):
    sp=diagonalize_ground(params)
    E0,w0=sp.energy_for_bare(1.5,2.5)
    output=[]
    for label,obs in GROUND_OBSERVED.items():
        E,w=sp.energy_for_bare(*obs['bare'])
        predicted=E-E0
        output.append({'state':label,'measured_ground_splitting_khz':obs['splitting_khz'],
             'measured_uncertainty_khz':obs['sigma_khz'],
             'predicted_ground_splitting_khz':predicted,
             'residual_khz':predicted-obs['splitting_khz'],
             'overlap':w, 'reference_overlap':w0})
    return output


def ground_sensitivity(params=NaCsParameters()):
    """Local finite-difference response of 3 measured ground splittings.

    These derivatives are diagnostic only; no molecular parameter fitting.
    Note that effective nuclear g-factors and field are correlated.
    """
    from dataclasses import replace
    directions=(('field_G',.01,'kHz/G'),
                ('c4_khz',.01,'kHz per kHz c4'),
                ('g_na',1e-4,'kHz per unit gNa'),
                ('g_cs',1e-4,'kHz per unit gCs'))
    result=[]
    for name,step,units in directions:
        lo=ground_splitting_comparison(replace(params,**{name:getattr(params,name)-step}))
        hi=ground_splitting_comparison(replace(params,**{name:getattr(params,name)+step}))
        for left,right in zip(lo,hi):
            result.append({'parameter':name,'state':left['state'],'units':units,
                           'derivative':(right['predicted_ground_splitting_khz']-
                           left['predicted_ground_splitting_khz'])/(2*step)})
    return result
