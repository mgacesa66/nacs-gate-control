"""Equal-control-resource comparison: independent versus joint optimization."""
from dataclasses import replace
from scipy.optimize import differential_evolution, minimize_scalar
from .propagate import Controls,gate_matrix,transfer_probability
from .metrics import metrics

# Same bounds for independent and joint control, no change of the frozen molecule.
BOUNDS=((.24,.45),(-35.,35.),(.35,1.20)) # ms, kHz, ms


def score(model,c):
    m,_=gate_matrix(model,c)
    return metrics(m)


def nominal(model):
    c=Controls()
    return c,score(model,c)


def independent(model, maxiter=18, popsize=5, seed=12):
    # First optimize transfer fidelity on isolated molecule, *not* a two-qubit gate.
    f=lambda x:1-transfer_probability(model,float(x[0]),float(x[1]))
    op=differential_evolution(f,BOUNDS[:2],seed=seed,maxiter=maxiter,
            popsize=popsize,polish=True,tol=.0002)
    t,det=op.x
    # Then optimize dwell globally for the actual complete two-qubit target;
    # rapid coherent phase structure makes a single bounded local search unsafe.
    res=differential_evolution(lambda x:1-score(model,Controls(transfer_ms=t,
                  transfer_detuning_khz=det,dwell_ms=float(x[0])))['fidelity'],
                  [BOUNDS[2]],seed=seed+31,maxiter=maxiter*2,popsize=max(10,popsize),
                  polish=True,tol=.00005)
    c=Controls(transfer_ms=float(t),transfer_detuning_khz=float(det),
               dwell_ms=float(res.x[0]))
    return c,score(model,c),{'one_body_transfer':float(transfer_probability(model,t,det))}


def joint(model,maxiter=24,popsize=5,seed=14):
    def obj(x):
        c=Controls(transfer_ms=x[0],transfer_detuning_khz=x[1],dwell_ms=x[2])
        return 1-score(model,c)['fidelity']
    res=differential_evolution(obj,BOUNDS,seed=seed,maxiter=maxiter,
           popsize=popsize,polish=True,tol=.0002)
    c=Controls(transfer_ms=float(res.x[0]),transfer_detuning_khz=float(res.x[1]),
               dwell_ms=float(res.x[2]))
    return c,score(model,c),{'one_body_transfer':float(transfer_probability(model,c.transfer_ms,c.transfer_detuning_khz))}
