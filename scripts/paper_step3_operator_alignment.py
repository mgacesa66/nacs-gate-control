"""Compare 11 and 25 logical maps, allowing only local static Z and global phase."""
import csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from nacs.analysis.matched_objectives import PairModel,frozen,VALID_BASIS,echo_ref_matrix,design_from_vector,gate_output,PHASE_RAMP_US,RAMP_STEPS
from nacs.data import ROOT

def fit_alignment(A,B):
    # Minimize || B - exp(i g) Z1(theta1)⊗Z2(theta2) A ||_F.
    d1=np.array([-1,-1,1,1],float)/2
    d2=np.array([-1,1,-1,1],float)/2
    def objective(z):
      ph=np.exp(1j*(z[0]+z[1]*d1+z[2]*d2))
      return float(np.linalg.norm(B-ph[:,None]*A)**2)
    seeds=[(0,0,0),(.1,.1,.1),(-.1,-.1,-.1)]
    result=min((minimize(objective,np.array(s),method='BFGS') for s in seeds),key=lambda x:x.fun)
    return np.sqrt(max(0,result.fun)),result.x.tolist()

def main():
    root=ROOT/'benchmarks/data/paper_step3'
    records=json.loads((root/'nominal/nominal_11.json').read_text())
    model=PairModel(frozen,VALID_BASIS)
    echo=echo_ref_matrix(model,steps=24)
    saved={};out=[]
    for item in records:
      k=item['protocol']
      x=gate_output(model,design_from_vector(item['design']),echo_ref=echo,steps=24,
                     phase_ramp_us=PHASE_RAMP_US,ramp_steps=RAMP_STEPS)
      A=x[np.array(model.logical_indices),:]
      obj=json.loads((root/f'basis/{k}_25_columns_steps24.json').read_text())
      B=np.asarray(obj['logical_matrix_real'])+1j*np.asarray(obj['logical_matrix_imag'])
      aligned,angles=fit_alignment(A,B)
      out.append({'protocol':k,'raw_frobenius':float(np.linalg.norm(A-B)),
         'aligned_frobenius':float(aligned),'global_and_local_z_rad':angles,
         'max_complex_entry_difference':float(np.max(abs(A-B)))})
      saved[k+'_11']=A;saved[k+'_25']=B
      print(k,'raw',round(out[-1]['raw_frobenius'],7),'aligned',round(aligned,7),flush=True)
    (root/'basis/operator_alignment_11_to_25.json').write_text(json.dumps(out,indent=2)+'\n')
    np.savez_compressed(root/'basis/logical_matrices_11_25.npz',**saved)
if __name__=='__main__':main()
