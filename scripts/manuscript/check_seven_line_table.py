"""Check the typeset spectroscopy values against the frozen static CSV.

Prevents accidental reuse of the CSV 'adjusted_residual_khz' column, which
includes inferred line-specific microwave shifts for three hyperfine lines.
This is a consistency check, not a new spectrum calculation.
"""
import argparse
from pathlib import Path
import csv
import re

ROOT=Path(__file__).resolve().parents[2]
source={r['line']:r for r in csv.DictReader((ROOT/'benchmarks'/'data'/'spectroscopy'/'seven_line_accounting.csv').open())}
parser = argparse.ArgumentParser()
parser.add_argument('--table', type=Path, default=ROOT/'manuscript_tables/spectroscopy.tex')
args = parser.parse_args()
latex = args.table.read_text()
m=re.search(r'\\begin\{tabular\}\{llrrrrr\}(.*?)\\end\{tabular\}',latex,re.S)
assert m,'seven-line table missing'
text=m.group(1)
names={'0->e':r'0\to e','0->f_z':r'0\to f_z','0->f_y':r'0\to f_y',
       '0->f_weak':r'0\to f_w','g_plus->e':r'g_+\to e',
       'g_minus->e':r'g_-\to e','1->e':r'1\to e'}

def val(st):
    return float(st.replace('$','').strip())
for name, label in names.items():
    rows=[l for l in text.splitlines() if l.startswith('$'+label+'$ &')]
    assert len(rows)==1,(name,rows)
    cols=[c.strip().replace(r'\\','').strip() for c in rows[0].split(' & ')]
    actual=[val(c) for c in cols[2:6]]
    r=source[name]
    expected=[float(r['experiment_khz']),float(r['static_model_khz']),
              float(r['raw_residual_khz']),float(r['raw_residual_khz'])-5.88]
    assert all(abs(a-e)<0.00051 for a,e in zip(actual,expected)),(name,actual,expected)
print('PASS: all seven rows match the canonical frozen undriven static export, without inferred microwave shifts')
