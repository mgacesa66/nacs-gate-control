"""Frozen-spectroscopy-backed, selected-state NaCs dipolar pair Hamiltonian.

Physical conventions: all Hamiltonians H/h are cyclic kHz; time ms. Never
fit molecular energies or dipoles here. Static export is immutable.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import epsilon_0, h, pi

FROZEN_SHA256 = 'e8c99f59e0cbab339f0f48ed23b21acba5fb81a86130f9297f89891b2aecdeb5'
DEBYE = 3.33564e-30
LABELS = ('0','1','e','e_prime','f_z','f_y','f_weak','g_plus','g_minus')
BASIS_SETS = {
    'minimal': ('0','1','e'),
    'four': ('0','1','e','e_prime'),
    'six': ('0','1','e','e_prime','f_z','f_y'),
    'nine': LABELS,
}


@dataclass
class FrozenSpectrum:
    energies: np.ndarray
    dipoles: tuple
    indices: dict
    checksum: str

    @classmethod
    def load(cls, root):
        root=Path(root)
        fn=root/'gate_ready.npz'
        actual=hashlib.sha256(fn.read_bytes()).hexdigest()
        if actual != FROZEN_SHA256:
            raise ValueError(f'Frozen gate_ready.npz checksum mismatch: {actual}')
        d=np.load(fn, allow_pickle=False)
        meta=json.loads((root/'baseline_manifest.json').read_text())
        if meta['version'] != 'spectroscopy-v1.0':
            raise ValueError('Unexpected spectroscopy baseline version')
        indices={k:int(v) for k,v in meta['state_labels'].items()}
        assert list(d['indices']) == [indices[k] for k in ('0','1','e','e_prime')]
        dm,dz,dp=(d['d_'+k+'_debye'] for k in ('minus1','zero','plus1'))
        dx=(dm-dp)/np.sqrt(2.0)
        dy=1j*(dm+dp)/np.sqrt(2.0)
        for mat in (dx,dy,dz):
            if not np.allclose(mat,mat.conj().T,atol=1e-10):
                raise ValueError('Dipole matrix not Hermitian')
        return cls(d['energy_khz'], (dx,dy,dz), indices, actual)


class PairModel:
    def __init__(self, frozen:FrozenSpectrum, labels=BASIS_SETS['four'],
                 distance_um=1.9, separation_xyz=(1.,0.,0.), scale_dd=1.0):
        self.frozen=frozen
        self.labels=tuple(labels)
        if not all(v in self.labels for v in ('0','1','e')): raise ValueError('Requires 0,1,e states')
        self.n=len(self.labels); self.size=self.n**2
        self.local={k:i for i,k in enumerate(self.labels)}
        take=np.array([frozen.indices[x] for x in self.labels],int)
        self.source_indices=take
        # Physics: all N=0 are the first 32 by energy; N=1 the remaining 96.
        self.n_exc=np.array([float(i>=32) for i in take])
        self.P_exc=np.diag(self.n_exc)
        self.num_exc=np.kron(self.P_exc,np.eye(self.n))+np.kron(np.eye(self.n),self.P_exc)
        self.num_exc_diag=np.diag(self.num_exc).real
        self.E0=float(frozen.energies[frozen.indices['0']])
        self.nu_ref=float(frozen.energies[frozen.indices['e']]-self.E0)
        self.energies=frozen.energies[take]-self.E0-self.nu_ref*self.n_exc
        self.d=tuple(mat[np.ix_(take,take)] for mat in frozen.dipoles)
        self.dipolar_factor=(DEBYE**2)/(4*pi*epsilon_0*(distance_um*1e-6)**3)/h/1e3*scale_dd
        rhat=np.array(separation_xyz,dtype=float)
        if not np.linalg.norm(rhat): raise ValueError('Zero separation vector')
        rhat/=np.linalg.norm(rhat)
        self.rhat=rhat
        # Keep only resonant N_exc-conserving dipolar terms (rotating wave approximation).
        dd=sum((float(i==j)-3*rhat[i]*rhat[j])*np.kron(self.d[i],self.d[j])
               for i in range(3) for j in range(3))*self.dipolar_factor
        conserves=np.isclose(self.num_exc_diag[:,None],self.num_exc_diag[None,:])
        self.Hdd=0.5*((dd*conserves)+(dd*conserves).conj().T)
        self.Hfree=(np.kron(np.diag(self.energies),np.eye(self.n))+
                    np.kron(np.eye(self.n),np.diag(self.energies))+self.Hdd)
        self.logical_indices=[self.pair(a,b) for a,b in (('0','0'),('0','1'),('1','0'),('1','1'))]

    def pair(self,a,b): return self.local[a]*self.n+self.local[b]

    def exchange_khz(self):
        return self.Hdd[self.pair('e','0'),self.pair('0','e')]

    def transition_khz(self,initial,target):
        i=self.local[initial];j=self.local[target]
        return (self.frozen.energies[self.source_indices[j]]-
                self.frozen.energies[self.source_indices[i]])

    def mw_hamiltonian(self, initial, target, f_target_khz, detuning_khz=0.,
                       epsilon=(1.,0.,0.), phase=0.):
        """Return two-molecule drive Hamiltonian in carrier rotating frame.

        Assumes one coherent microwave carrier, its phase reset to phase at
        each segment start (a controllable phase convention). All N0-N1
        matrix elements of SAME microwave drive are retained, not just target.
        """
        pol=np.array(epsilon,dtype=complex)
        if np.linalg.norm(pol)<1e-12: raise ValueError('Zero polarization')
        pol/=np.linalg.norm(pol)
        dpol=sum(pol[k]*self.d[k] for k in range(3))
        i,j=self.local[initial],self.local[target]
        if self.n_exc[i]==self.n_exc[j]: raise ValueError('Microwave requires N change')
        dtarget=abs(dpol[i,j])
        if dtarget<1e-9:raise ValueError('Target dipole nearly dark with specified polarization')
        mw=np.zeros((self.n,self.n),complex)
        for k in range(self.n):
            for l in range(self.n):
                if self.n_exc[k]==0 and self.n_exc[l]==1:
                    mw[k,l]=f_target_khz*dpol[k,l]*np.exp(1j*phase)/(2*dtarget)
                    mw[l,k]=mw[k,l].conjugate()
        # ν is cyclic kHz for selected transition; E already in ν_ref frame.
        nu=self.transition_khz(initial,target)+detuning_khz
        delta=nu-self.nu_ref
        h_local=np.diag(self.energies-delta*self.n_exc)+mw
        h_pair=np.kron(h_local,np.eye(self.n))+np.kron(np.eye(self.n),h_local)+self.Hdd
        return h_pair,delta, dtarget
