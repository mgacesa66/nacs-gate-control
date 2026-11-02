"""Independent NaCs N=0,1 hyperfine--Zeeman spectroscopy model.

Unit convention: Hamiltonian H/h in cyclic kHz, eigenvalues in kHz;
NOT angular frequencies. Dipoles in debye. Basis: |N,mN,mNa,mCs>.
No nonlinear/empirical adjustment of unknown transitions is performed.

For trap term use P2(epsilon.n) = [3 |epsilon.n|^2 - 1]/2, projected
into N<=1. It is an *effective* tensor optical Hamiltonian. Complex
polarization is supported; no vector (circular) polarizability is assumed.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import numpy as np
from scipy.linalg import eigh
from sympy.physics.wigner import wigner_3j, clebsch_gordan
from sympy import Rational
from .constants import NaCsParameters, MU_N_OVER_H_KHZ_PER_G


def _rat(x):
    return Rational(int(round(2*x)),2)

@lru_cache(maxsize=None)
def threej(a,b,c,d,e,f):
    return float(wigner_3j(*map(_rat,(a,b,c,d,e,f))))

@lru_cache(maxsize=None)
def cg(a, ma, b, mb, c, mc):
    return float(clebsch_gordan(_rat(a),_rat(b),_rat(c),_rat(ma),_rat(mb),_rat(mc)))


def spin_matrices(j):
    m=np.arange(j,-j-1,-1,dtype=float)
    z=np.diag(m).astype(complex)
    p=np.zeros((len(m),len(m)),complex)
    for idx, mi in enumerate(m):
        if idx>0:  # |m> -> |m+1>, lower row index
            p[idx-1,idx]=np.sqrt((j-mi)*(j+mi+1))
    minus=p.conj().T
    return {"x":(p+minus)/2, "y":(p-minus)/(2j),
            "z":z, "+":p, "-":minus, "0":np.eye(len(m),dtype=complex)}


def spherical_spin(j):
    s=spin_matrices(j)
    return {-1:s["-"]/np.sqrt(2), 0:s["z"], +1:-s["+"]/np.sqrt(2)}


def spin_rank2(j):
    sp=spherical_spin(j)
    return {q:sum((cg(1,q1,1,q-q1,2,q)*sp[q1]@sp[q-q1]
                   for q1 in (-1,0,1) if q-q1 in sp),
                  np.zeros((int(2*j+1),)*2,complex)) for q in range(-2,3)}


def cross_spin_rank2(j1,j2):
    a,b=spherical_spin(j1),spherical_spin(j2)
    out={}
    for q in range(-2,3):
        x=np.zeros((int((2*j1+1)*(2*j2+1)),)*2,complex)
        for qa in (-1,0,1):
            qb=q-qa
            if qb in b:
                x+=cg(1,qa,1,qb,2,q)*np.kron(a[qa],b[qb])
        out[q]=x
    return out


@dataclass
class DressedSpectrum:
    params: NaCsParameters
    basis: tuple[tuple[float,float,float,float],...]
    energies_khz: np.ndarray
    vectors: np.ndarray    # columns are eigenstates, phases follow scipy eigh
    dipole_q_debye: dict[int,np.ndarray]  # spherical q=-1,0,+1 in dressed basis
    H_khz: np.ndarray

    def bare_index(self, N, mN, mNa, mCs):
        value=(float(N),float(mN),float(mNa),float(mCs))
        return self.basis.index(value)

    def dominant_components(self, eigen_index, how_many=5):
        v=self.vectors[:,eigen_index]
        ids=np.argsort(np.abs(v)**2)[::-1][:how_many]
        return [(self.basis[i],float(abs(v[i])**2)) for i in ids]

    def best_overlap(self, vector, allowed_N=None, exclude=()):
        """Return eigenstate maximizing |<psi_dressed|vector_bare>|^2."""
        vector=np.asarray(vector,complex)
        assert vector.shape==(len(self.basis),)
        projection=np.abs(self.vectors.conj().T@vector)**2
        for i in exclude: projection[i]=-1
        if allowed_N is not None:
            for i in range(len(self.basis)):
                pN=sum(abs(self.vectors[k,i])**2 for k,st in enumerate(self.basis) if st[0]==allowed_N)
                if pN<.5: projection[i]=-1
        best=int(np.argmax(projection))
        return best,float(projection[best])

    def dressed_dipole(self, a,b,polarization_q):
        """<a|epsilon . d|b> using polarization_q coefficients of d_q.
        Coefficients are explicit, not inferred from an experimental horn.
        This returns complex Debye for the supplied linear combination.
        """
        return sum(polarization_q.get(q,0)*self.dipole_q_debye[q][a,b] for q in (-1,0,1))

    def export(self,path):
        """Stable integration handoff; energy and dressed dipoles share basis ordering."""
        np.savez_compressed(path,
            energies_khz=self.energies_khz,vectors=self.vectors,
            d_q_minus1_debye=self.dipole_q_debye[-1],
            d_q_zero_debye=self.dipole_q_debye[0],
            d_q_plus1_debye=self.dipole_q_debye[1],
            bare_basis_quantum_numbers=np.asarray(self.basis),
            field_G=self.params.field_G,
            optical_tensor_khz=self.params.optical_tensor_khz,
            dipole_debye=self.params.dipole_debye)


class NaCsMolecule:
    """Build spectroscopy from documented constants; no two-body dynamics."""
    def __init__(self,params:NaCsParameters=NaCsParameters()):
        self.p=params
        self.rot=tuple((float(n),float(m)) for n in (0,1)
                       for m in range(n,-n-1,-1))
        self.nrot=len(self.rot)   # 4
        self.nna=4; self.ncs=8
        self.basis=tuple((n,m,na,cs) for (n,m),na,cs in product(
            self.rot,np.arange(1.5,-1.51,-1),np.arange(3.5,-3.51,-1)))
        self.eye_r=np.eye(4,dtype=complex)
        self.eye_na=np.eye(4,dtype=complex)
        self.eye_cs=np.eye(8,dtype=complex)
        self.eye_nuclear=np.eye(32,dtype=complex)
        self.eye=np.eye(128,dtype=complex)
        self.Sna=spin_matrices(1.5)
        self.Scs=spin_matrices(3.5)
        self.C={k:{q:self.rotational_C(k,q) for q in range(-k,k+1)} for k in (1,2)}
        self.N={c:self._rot_angular_momentum(c) for c in 'xyz'}
        self.Tna=spin_rank2(1.5)
        self.Tcs=spin_rank2(3.5)
        self.Tcross=cross_spin_rank2(1.5,3.5)

    def rotational_C(self,k,q):
        """Racah-normalized C^(k)_q=sqrt(4pi/(2k+1))Y_kq on rotor."""
        mat=np.zeros((4,4),complex)
        for a,(n2,m2) in enumerate(self.rot):
            for b,(n1,m1) in enumerate(self.rot):
                if abs(m2-q-m1)>1e-9:continue
                mat[a,b]=((-1)**int(m2))*np.sqrt((2*n2+1)*(2*n1+1))*\
                    threej(n2,k,n1,0,0,0)*threej(n2,k,n1,-m2,q,m1)
        return mat

    def _rot_angular_momentum(self,component):
        mat=np.zeros((4,4),complex)
        for n in (0,1):
            idx=[k for k,(nn,m) in enumerate(self.rot) if nn==n]
            block=spin_matrices(n)[component]
            mat[np.ix_(idx,idx)]=block
        return mat

    def kron3(self,r,na,cs):
        return np.kron(np.kron(r,na),cs)

    def rank2_rot_nuc(self,T):
        """sum_q (-1)^q C2_q tensor T2_-q; Hermitian."""
        A=np.zeros((128,128),complex)
        for q in range(-2,3):
            A+=(-1)**q * np.kron(self.C[2][q],T[-q])
        return A

    def trap_P2(self):
        """Projection of [3|eps.n|^2-1]/2 on N<=1 (4x4)."""
        eps=np.asarray(self.p.polarization_xyz,complex)
        if np.linalg.norm(eps)==0:raise ValueError('optical polarization cannot vanish')
        eps/=np.linalg.norm(eps)
        # Cartesian expectation e*_i e_j n_i n_j; only real symmetric part survives.
        A=np.real(np.outer(np.conj(eps),eps))
        I=np.eye(4)
        c0,c1,c_1,c2,c_2=self.C[2][0],self.C[2][1],self.C[2][-1],self.C[2][2],self.C[2][-2]
        zz=(I+2*c0)/3
        xx_plus_yy=I-zz
        xx_minus_yy=np.sqrt(2/3)*(c2+c_2)
        xx=(xx_plus_yy+xx_minus_yy)/2
        yy=(xx_plus_yy-xx_minus_yy)/2
        xy=(c2-c_2)/(1j*np.sqrt(6))
        xz=(c_1-c1)/np.sqrt(6)
        yz=(c_1+c1)/(-1j*np.sqrt(6))
        nn=np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
        op=sum(A[i,j]*nn[i,j] for i in range(3) for j in range(3))
        return (3*op-I)/2

    def terms_khz(self):
        p=self.p
        terms={}
        E=np.array([p.B_rot_khz*n*(n+1)-p.D_rot_khz*(n*(n+1))**2
                    for n,m in self.rot],float)
        terms['rotation']=np.kron(np.diag(E),self.eye_nuclear)
        # Nuclear Zeeman, including nuclear shielding. HZ = -g mu_N B Iz.
        mn=MU_N_OVER_H_KHZ_PER_G * p.field_G
        hz_na=self.kron3(self.eye_r,self.Sna['z'],self.eye_cs)
        hz_cs=self.kron3(self.eye_r,self.eye_na,self.Scs['z'])
        terms['nuclear_zeeman']=-mn*(p.g_na*(1-p.shielding_na)*hz_na +
                                    p.g_cs*(1-p.shielding_cs)*hz_cs)
        terms['rotational_zeeman']=-mn*p.g_rot*self.kron3(self.N['z'],self.eye_na,self.eye_cs)
        sr=np.zeros((128,128),complex)
        for comp in 'xyz':
            sr+=p.c_na_khz*self.kron3(self.N[comp],self.Sna[comp],self.eye_cs)
            sr+=p.c_cs_khz*self.kron3(self.N[comp],self.eye_na,self.Scs[comp])
        terms['spin_rotation']=sr
        ss4=np.zeros((128,128),complex)
        for comp in 'xyz':
            ss4+=self.kron3(self.eye_r,self.Sna[comp],self.Scs[comp])
        terms['scalar_spin_spin']=p.c4_khz*ss4
        terms['quadrupole_na']= (p.eQq_na_khz*np.sqrt(6)/(4*1.5*(2*1.5-1))
            * self.rank2_rot_nuc({q:np.kron(t,self.eye_cs) for q,t in self.Tna.items()}))
        terms['quadrupole_cs']= (p.eQq_cs_khz*np.sqrt(6)/(4*3.5*(2*3.5-1))
            * self.rank2_rot_nuc({q:np.kron(self.eye_na,t) for q,t in self.Tcs.items()}))
        # c3 I_Na . (I -3 n n) . I_Cs (anisotropic nuclear spin-spin)
        terms['tensor_spin_spin']=-p.c3_khz*np.sqrt(6)*self.rank2_rot_nuc(self.Tcross)
        terms['optical_tensor']=p.optical_tensor_khz*np.kron(self.trap_P2(),self.eye_nuclear)
        terms['optical_scalar']=p.optical_scalar_khz*self.eye
        return terms

    def hamiltonian_khz(self):
        terms=self.terms_khz()
        H=sum(terms.values(),np.zeros((128,128),complex))
        if not np.allclose(H,H.conj().T,atol=1e-8):
            raise RuntimeError('Hamiltonian non-Hermitian: check tensor conventions')
        return H

    def bare_dipole_q_debye(self):
        return {q:self.p.dipole_debye*np.kron(self.C[1][q],self.eye_nuclear)
                for q in (-1,0,1)}

    def diagonalize(self):
        H=self.hamiltonian_khz()
        E,V=eigh(H,check_finite=True)
        # Fix the otherwise arbitrary U(1) gauge of every eigenvector:
        # its largest bare-basis component is real and positive. This makes
        # exported complex dipoles reproducible between runs (away from crossings).
        pivots=np.argmax(abs(V),axis=0)
        for col,row in enumerate(pivots):
            phase=V[row,col]
            if abs(phase)>0: V[:,col]*=np.conj(phase)/abs(phase)
        D={q:V.conj().T@d@V for q,d in self.bare_dipole_q_debye().items()}
        return DressedSpectrum(self.p,self.basis,E,V,D,H)

    def bare_state(self,N,mN,mNa,mCs):
        v=np.zeros(128,complex)
        v[self.basis.index((float(N),float(mN),float(mNa),float(mCs)))]=1
        return v

    def superposition(self,terms):
        """terms: tuples (amplitude, N, mN, mNa, mCs)."""
        v=sum((a*self.bare_state(N,mN,mNa,mCs) for a,N,mN,mNa,mCs in terms),
              np.zeros(128,complex))
        return v/np.linalg.norm(v)

    def identify_gate_states(self,spectrum=None,orientation='x'):
        """Identify logical 0,1 and x/y/z-oriented excited states by overlap.

        These are *diabatic assignments* and must be checked with measured
        spectroscopy; they are not guaranteed at avoided crossings.
        """
        sp=spectrum if spectrum is not None else self.diagonalize()
        zero=self.bare_state(0,0,1.5,2.5)
        one=self.bare_state(0,0,-.5,2.5)
        # For N=1, orbitally x/y aligned eigenvectors of tensor light field.
        def n1_vector(mna,mcs,axis):
            if axis=='x':
                return self.superposition([(-1/np.sqrt(2),1,1,mna,mcs),
                                           (1/np.sqrt(2),1,-1,mna,mcs)])
            if axis=='y':
                return self.superposition([(1/np.sqrt(2),1,1,mna,mcs),
                                           (1/np.sqrt(2),1,-1,mna,mcs)])
            if axis=='z':return self.bare_state(1,0,mna,mcs)
            raise ValueError(axis)
        choices={}
        for name,v in [('0',zero),('1',one),('e',n1_vector(1.5,2.5,orientation)),
                       ('e_prime',n1_vector(-.5,2.5,orientation)),
                       ('f_z',n1_vector(1.5,2.5,'z')),
                       ('f_y',n1_vector(1.5,2.5,'y'))]:
            index,overlap=sp.best_overlap(v,allowed_N=(0 if name in ('0','1') else 1))
            choices[name]={'index':index,'bare_overlap':overlap,
                           'energy_khz':float(sp.energies_khz[index])}
        return choices
