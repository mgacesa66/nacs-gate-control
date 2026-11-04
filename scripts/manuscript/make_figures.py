"""Generate manuscript figures only from frozen Step 1--3 data.

This script does not modify any frozen source, archive, or numerical result.
Usage: python make_figures.py --input /path/to/source_frozen/benchmarks/data --output /path/to/figures
"""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args()
base=Path(a.input); out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':9,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,'figure.dpi':140,'savefig.bbox':'tight'})
# Fig 1: methodological workflow; schematic is not measured data.
fig, ax = plt.subplots(figsize=(10.4,3.25))
fig.subplots_adjust(left=.02,right=.98,bottom=.02,top=.98)
ax.set(xlim=(0,10.4),ylim=(0,3.25));ax.axis('off')
edge='#4b6078';arrow='#3a4a5e'
boxes=[(.12,2.12,2.16,.88,'Spectroscopy\nreference +\neffective fit'),(2.80,2.12,2.16,.88,'128-state dressed\nenergies and dipoles'),(5.48,2.12,2.16,.88,'Microwave + dipolar\nstate screening'),(8.16,2.12,2.12,.88,'Coherent pair\npropagation'),(3.66,.64,3.08,.86,'Logical map, leakage\nrecalibration checks')]
texts=[]
for x,y,w,h,s in boxes:
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.025,rounding_size=0.09',linewidth=1.1,edgecolor=edge,facecolor='#e8eef4'))
 texts.append((ax.text(x+w/2,y+h/2,s,ha='center',va='center',fontsize=11,linespacing=1.25), (x,y,w,h)))
for xa,xb in [(2.28,2.8),(4.96,5.48),(7.64,8.16)]:
 ax.annotate('',xy=(xb-.06,2.56),xytext=(xa+.06,2.56),arrowprops=dict(arrowstyle='->',lw=1.25,color=arrow))
ax.annotate('',xy=(6.82,1.07),xytext=(9.22,2.04),arrowprops=dict(arrowstyle='->',lw=1.25,color=arrow,connectionstyle='angle,angleA=-90,angleB=0,rad=15'))
ax.text(5.2,.14,'Scope: fixed geometry and coherent internal dynamics; no molecular motion or irreversible loss',ha='center',va='center',fontsize=9.5,color='#545454')
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
for t,(x,y,w,h) in texts:
 b=t.get_window_extent(renderer); l,bottom=ax.transData.transform((x,y));r,top=ax.transData.transform((x+w,y+h))
 assert b.x0>l and b.x1<r and b.y0>bottom and b.y1<top, t.get_text()
fig.savefig(out/'fig1_model_chain.pdf',dpi=220,bbox_inches='tight',pad_inches=.12)

plt.close(fig)
# Fig 2 calibrated ellipticity scan
p1=pd.read_csv(base/'paper_step1/gate_comparison.csv');m={'nominal':0.,'ellipse_minus_0p1deg':-.1,'ellipse_plus_0p1deg':.1,'ellipse_minus_0p25deg':-.25,'ellipse_plus_0p25deg':.25,'ellipse_minus_0p5deg':-.5,'ellipse_plus_0p5deg':.5}
sel=p1[(p1.protocol=='all_spectators') & p1.probe.isin(m)].copy();sel['epsilon']=sel.probe.map(m)
fig,ax=plt.subplots(figsize=(5.9,3.2)); styles=[('fixed','Fixed carriers','o','-'),('frequency','Retuned carriers','s','--'),('frequency_rabi','Carriers + target Rabi','^','-.')]
for mode,label,marker,linestyle in styles:
 d=sel[sel['mode']==mode].sort_values('epsilon');ax.plot(d.epsilon,100*d.fidelity_phase_refit,marker=marker,linestyle=linestyle,lw=1.4,ms=4,label=label)
ax.set_xlabel('Change in optical ellipticity (degrees)');ax.set_ylabel('Unconditional logical fidelity (%)');ax.set_xlim(-.53,.53);ax.set_ylim(0,100);ax.axvline(0,color='#777',lw=.8,ls=':');ax.grid(alpha=.20);ax.legend(loc='lower center',frameon=True,ncol=1);fig.savefig(out/'fig2_ellipticity_recalibration.pdf');plt.close(fig)

# Fig3 nonmonotonic basis
p2=pd.read_csv(base/'paper_step2/paper_step2_extended/extended_nominal.csv'); d=p2[(p2.probe=='nominal') & (p2.protocol=='all_spectators')].sort_values('basis_n');
punion=pd.read_csv(base/'paper_step2/paper_step2_union/cross_screen.csv');u=punion[(punion.probe=='nominal')&(punion.protocol=='all_spectators')&(punion.scheme=='union')].iloc[0]
fig,ax=plt.subplots(figsize=(5.9,3.0));ax.plot(d.basis_n,100*d.fidelity,marker='o',ms=4,lw=1.5,label='Microwave-ranked states');ax.scatter([u.basis_n],[100*u.fidelity],s=55,marker='D',label='Microwave + dipolar screen',zorder=4);ax.axhline(100*d.iloc[0].fidelity,color='#888',ls='--',lw=.9,label='11-state reference');ax.set_xlabel('Selected states per molecule');ax.set_ylabel('Unconditional logical fidelity (%)');ax.set_xticks([11,15,17,19,21,23,25]);ax.set_ylim(94.27,94.72);ax.grid(alpha=.2);ax.legend(loc='lower left',frameon=True);fig.savefig(out/'fig3_nonmonotonic_basis.pdf');plt.close(fig)

import json
import matplotlib.pyplot as plt
from itertools import combinations
plt.rcParams['pdf.fonttype']=42
f,a=plt.subplots(figsize=(7.2,3.8),layout='constrained')
rows=json.load(open(base/'paper_step3/nominal/nominal_11.json'))
labels=[('Transfer A',98.49,82.95),('Transfer B',98.50,90.7),('Selected spectator',98.22,96.6),('Logical channel',98.07,92.1),('Total-spectator guard',98.48,96.6),('Gaussian reference',98.30,89.6)]
ts=[]
for i,(r,(s,tx,ty)) in enumerate(zip(rows,labels)):
 x=50*(r['transfer_forward']+r['transfer_reverse']);y=100*r['fidelity'];a.scatter(x,y,marker='o' if i<2 else 's',s=35,zorder=4)
 a.annotate('',(x,y),(tx,ty),arrowprops=dict(arrowstyle='-',color='.48',lw=.7,shrinkA=6,shrinkB=4))
 ts.append(a.text(tx,ty,s,ha='center',va='center',fontsize=9,backgroundcolor='white'))
a.set(xlim=(97.9,98.63),ylim=(79,98),xlabel='Mean isolated forward/reverse transfer probability (%)',ylabel='Unconditional full-gate fidelity (%)');a.grid(alpha=.2)
f.canvas.draw();assert all(not x.get_window_extent().overlaps(y.get_window_extent()) for x,y in combinations(ts,2))
f.savefig(out/'fig4_objective_nonidentifiability.pdf')
