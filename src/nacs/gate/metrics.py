"""Unconditional average gate fidelity, including leakage and permitted Z frames."""
import numpy as np
from scipy.optimize import minimize

# State order |00>, |01>, |10>, |11>.
U_EXPERIMENTAL=np.array([[0,0,0,-1],
                         [0,-1j,0,0],
                         [0,0,-1j,0],
                         [-1,0,0,0]],complex)


def metrics(M, target=U_EXPERIMENTAL, phase_correct=True):
    # Projected logical map M need NOT be trace preserving. Its missing norm is leakage.
    survival=float(np.trace(M.conj().T@M).real/4)
    leak=1-survival
    A=target.conj().T
    zz1=np.array([-1,-1,1,1],float)/2
    zz2=np.array([-1,1,-1,1],float)/2
    def trace_with_phases(x):
        phased=np.exp(1j*(x[0]*zz1+x[1]*zz2))
        return np.trace(A@(phased[:,None]*M))
    if phase_correct:
        # grid init prevents local-maxima choices depending on input seed
        seeds=[(a,b) for a in (0,np.pi/2,np.pi,3*np.pi/2)
                    for b in (0,np.pi/2,np.pi,3*np.pi/2)]
        seed=max(seeds,key=lambda s: abs(trace_with_phases(s)))
        sol=minimize(lambda x:-abs(trace_with_phases(x))**2,
                     np.asarray(seed),method='BFGS',options={'gtol':1e-8,'maxiter':80})
        phi=sol.x
    else:phi=np.zeros(2)
    overlap=trace_with_phases(phi)
    f=(np.trace(M.conj().T@M).real+abs(overlap)**2)/20
    return {'fidelity':float(np.real(f)),'infidelity':float(1-np.real(f)),
            'survival':float(survival),'leakage':float(leak),
            'survival_normalized_score':float((np.real(f)/(survival)) if survival>0 else 0),
            'local_z_angles_rad':phi.tolist(),
            'trace_overlap':float(abs(overlap))}
