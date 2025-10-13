import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as m

plt.rc('font', family='Helvetica', size=15)
plt.rc('mathtext', fontset='dejavusans')

fig = plt.figure(figsize=(4.5, 4.5))
ax = plt.subplot()

plt.subplots_adjust(wspace=0, hspace=0, left=0.12,
                    top=0.95, right=0.95, bottom=0.12)
ax.tick_params(direction='in', top=False, right=False)
ax.tick_params(which='minor', direction='in', top=False, right=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

import math as m

def p_d(d, p, N):
    pd = (1-(1 - p**d)**m.comb(N-2, d-1))
    sums = 0
    for i in range(1, d): 
        sums += p_d(i, p, N)
    pd *= (1 - sums)
    return pd

Ns = [1e3, 1e4, 1e5, 1e6, 1e7]
Ns = [int(N) for N in Ns]
dmax = 2

# color = [(244,28,84),(255,159,0),(6,200,100),(0,170,240),(50,80,250)]
# color = color[::-1]
# color = [(r/255, g/255, b/255) for r, g, b in color]

# inferno color list 
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
cmap = cm.get_cmap("Blues")
norm = Normalize(vmin=0, vmax=len(Ns)+1)
color = [cmap(norm(i+2)) for i in range(len(Ns))]


p = np.linspace(0, .1, 10000)

# a = 0.645
# ax.plot([-1, 1], [a, a], linestyle='--', color='grey', lw=2)

for N in Ns:

    pc = N**(-(dmax-1)/dmax)

    pds = []
    for d in range(1, dmax+1):
        pd = p_d(d, p, N)
        pds.append(pd)

    nu = np.zeros(len(p))
    for j in range(len(pds)):
        nu += pds[j]
    
    ax.plot((p-pc), nu, label=f'$N = 10^{{{int(np.log10(N))}}}$', lw=2, color=color[Ns.index(N)])


ax.set_ylabel('$\\nu$', fontfamily='Times', labelpad=-15)
ax.set_xlabel('$p - p_c$', fontfamily='Times', labelpad=0)


ax.set_xlim(-0.03, 0.03)
ax.set_xticks([-0.02, 0, 0.02])
# ax.set_xticklabels(["-$10^{\\text{-}2}$", "0", "$10^{\\text{-}2}$"])
ax.set_xticklabels(["-0.02", "0", "0.02"])

ax.set_ylim(0, 1.01)
ax.set_yticks([0, 0.3, 0.7, 1])
ax.set_yticklabels(["0", "0.3", "0.7", "1"])

# ax.set_yscale('log')
# ax.set_ylim(0.1,2)
# ax.set_yticks([0.1,1])
# ax.set_yticklabels(["$10^{\\text{-}1}$", "$10^0$"])

ax.xaxis.set_tick_params(pad=5)

ax.legend(loc='upper left', fontsize=13, frameon=False)

plt.savefig(f'plots/fig6-phasetransition.pdf')
plt.close()