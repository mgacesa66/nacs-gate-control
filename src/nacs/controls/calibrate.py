"""Isolated-molecule calibration of echo for a static antenna orientation.

Echo calibration optimizes *measurable* carrier detuning and pulse duration;
it does not modify dipoles, power or underlying molecular constants.
"""
import numpy as np
from scipy.optimize import minimize
from .fixed_polarization import local_transfer, echo_rabi,rotated_linear_polarization

def calibrate_echo(model,angle_deg, max_echo_ms=.080):
    fr=echo_rabi(model,angle_deg)
    t0=1./(2.*fr)
    if t0>max_echo_ms:
        return None
    pol=rotated_linear_polarization(angle_deg)
    bounds=[(max(.005,t0*.72),min(max_echo_ms,t0*1.30)),(-25,25)]
    def objective(x):return -local_transfer(model,'0','e',fr,x[0],x[1],pol)
    # Multiple initial detunings protect against nearby spectator-induced shifts.
    vals=[]
    for d in [-10.,-3.,0.,3.,10.]:
        r=minimize(objective,[t0,d],bounds=bounds,method='Nelder-Mead',options={'maxiter':55,'xatol':1e-5,'fatol':1e-8})
        if np.isfinite(r.fun):vals.append(r)
    r=min(vals,key=lambda x:x.fun)
    return {'angle_deg':float(angle_deg),'echo_rabi_khz':float(fr),'ideal_pi_time_ms':float(t0),
            'echo_ms':float(r.x[0]),'echo_detuning_khz':float(r.x[1]),
            'echo_population':float(-r.fun)}
