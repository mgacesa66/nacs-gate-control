"""Frozen baseline and benchmark paths. No archive extraction at import time."""
from pathlib import Path
from functools import lru_cache
import hashlib
from nacs.gate.model import FrozenSpectrum

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'benchmarks'
EXPECTED='e8c99f59e0cbab339f0f48ed23b21acba5fb81a86130f9297f89891b2aecdeb5'

def benchmark_path(relative):
    return DATA/'data'/relative

@lru_cache(maxsize=1)
def frozen_spectrum():
    fn=DATA/'frozen'/'gate_ready.npz'
    found=hashlib.sha256(fn.read_bytes()).hexdigest()
    if found!=EXPECTED:
        raise ValueError(f'frozen spectroscopy data checksum mismatch {found}')
    return FrozenSpectrum.load(DATA/'frozen')
