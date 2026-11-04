"""OPTIONAL effective-tweezer calibration to 3 measured rotational lines.

This is explicitly a FIT, not a validation of the three fitted lines.
The remaining 1->e hyperfine line is withheld as a basic out-of-fit check.
The effective tensor parameters have no established experimental uniqueness.
"""
from __future__ import annotations
import argparse
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .constants import NaCsParameters
from .core import NaCsMolecule
from .validation import benchmark_lines


def make_params(u_khz,ellipticity_deg,base=NaCsParameters()):
    a=np.deg2rad(ellipticity_deg)
    return replace(base,optical_tensor_khz=float(u_khz),
                   polarization_xyz=(complex(np.cos(a)),complex(0,np.sin(a)),0j))


def fit_three_bright_lines(base=NaCsParameters()):
    """Fit 2 nuisance optical parameters to 3 line frequencies; return result.

    Fits 0->e, 0->f_z, 0->f_y, holding out hyperfine-changing 1->e.
    Requires correct mapping of x,y,z branches to the observed transitions.
    """
    def residuals(x):
        model=NaCsMolecule(make_params(*x,base))
        sp=model.diagonalize()
        lines=benchmark_lines(model,sp)
        return np.asarray([lines[k]['residual_khz'] for k in (0,1,2)])
    result=least_squares(residuals,[-3500,54.7],
                         bounds=([-6000,48],[ -1000,85]),
                         xtol=1e-9,ftol=1e-9,gtol=1e-9,max_nfev=40)
    return result


def run(argv=None):
    parser=argparse.ArgumentParser(description='Fit effective optical tensor to 3 Picard transitions (NOT validation)')
    parser.add_argument('--output',type=Path,default=Path('example_optical_fit'))
    parser.add_argument('--B-G',type=float,default=864.1)
    args=parser.parse_args(argv)
    base=NaCsParameters(field_G=args.B_G)
    res=fit_three_bright_lines(base)
    params=make_params(*res.x,base)
    model=NaCsMolecule(params)
    sp=model.diagonalize()
    lines=benchmark_lines(model,sp)
    print('TWO-PARAMETER FIT to 3 bright transition frequencies, not independent validation:')
    print(f'optical_tensor_khz={res.x[0]:.7f}, ellipticity_deg={res.x[1]:.7f}')
    for k,line in enumerate(lines):
        marker='FIT' if k<3 else 'HOLDOUT'
        print(f'{marker:7s} {line["name"]:8s} measured={line["frequency_offset_khz"]:+10.3f} '
              f'computed={line["predicted_offset_khz"]:+10.3f} residual={line["residual_khz"]:+10.3f} kHz')
    args.output.mkdir(parents=True,exist_ok=True)
    payload={'u_t_khz':float(res.x[0]),'ellipticity_deg':float(res.x[1]),
             'fit_lines':['0->e','0->f_z','0->f_y'],
             'holdout_line':'1->e','cost_khz2':float(2*res.cost),
             'scientific_status':'effective-tweezer fit only, not independent experimental calibration',
             'hyperfine_constants':'2017 Aldegunde/Hutson theoretical predictions'}
    (args.output/'optical_fit.json').write_text(json.dumps(payload,indent=2))
    print('Fit parameters saved in',args.output/'optical_fit.json')
    return res

if __name__=='__main__':run()
