"""Deterministic, deliberately non-quantitative NaCs level/basis schematic.

No numerical spectrum is recomputed. All vertical separations are illustrative.
This figure is a conceptual visualization of *static* Hamiltonian dressing,
not a Floquet/microwave-power-dressed spectrum.
"""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
OUT = args.output / 'fig5_level_schematic.pdf'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9, 'pdf.fonttype':42})
fig,axes=plt.subplots(1,3,figsize=(10.0,3.45),gridspec_kw={'width_ratios':[1,1.25,1.16]})
fig.subplots_adjust(left=.035,right=.975,top=.86,bottom=.13,wspace=.24)

for ax in axes:
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.axis('off')

# Bare rotational basis
ax=axes[0]
ax.set_title('(a) Rotational basis',fontsize=10,pad=10,weight='bold')
for y,txt,n in [(0.25,'$N=0$',32),(.70,'$N=1$',96)]:
    ax.hlines(y,.13,.88,lw=2,color='#2b465f')
    ax.text(.51,y+.055,txt,ha='center',fontsize=10)
    ax.text(.51,y-.065,f'{n} hyperfine basis states',ha='center',fontsize=8.6)
ax.text(.51,.93,'Before static-field diagonalization',ha='center',fontsize=8.3,color='#555555')
ax.text(.51,.055,'32 + 96 = 128 states',ha='center',weight='bold')

# Static eigenstates (schematic, non-energy scale)
ax=axes[1]
ax.set_title('(b) Static-field eigenstates',fontsize=10,pad=10,weight='bold')
ax.text(.50,.92,'Hyperfine + Zeeman + optical tensor',ha='center',fontsize=8.3,color='#555555')
# illustrative sets, no implied energy differences
lo=[(.15,'$g_{-}$'),(.22,'$0$'),(.29,'$g_{+}$'),(.36,'$1$')]
hi=[(.57,'$f_y$'),(.64,'$e$'),(.71,'$f_z$'),(.78,'$f_w$')]
for y,l in lo+hi:
    ax.hlines(y,.16,.67,lw=1.4,color='#2b465f')
    ax.text(.12,y,f'${l.strip("$")}$',ha='right',va='center',fontsize=9)
ax.text(.72,.43,r'$\vdots$',ha='center',fontsize=15)
ax.text(.72,.85,r'$\vdots$',ha='center',fontsize=15)
# microwave drive arrows indicate possible couplings, not actual line positions
ax.add_patch(FancyArrowPatch((.76,.22),(.76,.64),arrowstyle='<->',mutation_scale=10,lw=1.4,color='#c36a37'))
ax.add_patch(FancyArrowPatch((.86,.36),(.86,.64),arrowstyle='<->',mutation_scale=10,lw=1.25,linestyle='--',color='#c36a37'))
ax.text(.68,.09,'Solid: bright / dashed: weak',fontsize=7.8,ha='center',color='#7b5038')
ax.text(.50,.03,'128 static dressed levels',ha='center',weight='bold')

# Selected gate models
ax=axes[2]
ax.set_title('(c) Selected gate bases',fontsize=10,pad=10,weight='bold')
ax.text(.50,.92,'For each molecule; same 128-state source',ha='center',fontsize=8.3,color='#555555')
for y,n,title,pair,clr in [(.65,11,'Corrected core',121,'#2b465f'),(.29,25,'Cross-screened endpoint',625,'#467c70')]:
    ax.add_patch(Rectangle((.06,y-.115),.88,.23,fill=False,edgecolor=clr,lw=1.6))
    ax.text(.5,y+.055,title,ha='center',fontsize=9.0)
    ax.text(.5,y-.005,f'{n} states / molecule',ha='center',weight='bold',fontsize=9.5)
    ax.text(.5,y-.07,f'{pair} two-molecule basis states',ha='center',fontsize=8.3)
ax.annotate('',xy=(.50,.41),xytext=(.50,.51),arrowprops={'arrowstyle':'->','lw':1.2})
ax.text(.5,.03,'Additional states: microwave + dipolar screens',ha='center',fontsize=7.9)

fig.suptitle('From rotational basis to static dressed states and reduced gate models',fontsize=11.2)
OUT.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(OUT,bbox_inches='tight')
plt.close(fig)
print('wrote',OUT)
