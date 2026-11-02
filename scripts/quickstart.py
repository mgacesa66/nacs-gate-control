"""Compute a complete gate benchmark with the canonical 11-state model.

Use after `python -m pip install -e .`. Does not modify frozen data or figures.
"""
import json
from nacs.data import frozen_spectrum,benchmark_path
from nacs.gate.model import PairModel
from nacs.controls.fixed_polarization import VALID_BASIS
from nacs.controls.spectator import Design,gate_output,measures

spectrum=frozen_spectrum()
model=PairModel(spectrum,VALID_BASIS)
candidates=json.loads(benchmark_path('optimization/candidates.json').read_text())
for label in ('unmodulated','aux51','all_spectators'):
    d=Design(**candidates[label]['design'])
    m=measures(model,gate_output(model,d))
    print(f'{label:16} F_avg={m["fidelity"]:.6f}  leakage={m["leakage"]:.6f} '
          f'aux51={m["aux51"]:.6f}  other_spectators={m["other_spectators"]:.6f}')
