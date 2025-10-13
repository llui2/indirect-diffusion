import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rc('font', family='Helvetica', size=18)
plt.rc('mathtext', fontset='dejavusans')

file = 'fig4-degreedist-RR'
df = pd.read_csv(f'scripts/{file}.csv')

d1 = df['deg1'].to_numpy()
d2 = df['deg2'].to_numpy()

d2 = d1 + d2

fig, ax = plt.subplots(1, 1, figsize=(4, 4))

plt.subplots_adjust(wspace=0, hspace=0, left=0.2,
                    top=0.95, right=0.93, bottom=0.15)
ax.tick_params(direction='in', top=False, right=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

w1 = np.ones_like(d1)/len(d1)
w2 = np.ones_like(d2)/len(d2)

# ax.hist(d2, weights=w2, bins=20, color='red', alpha=1, lw=2, histtype='step', density=True)
ax.hist(d2, weights=w2, bins=6, color='#648fff', alpha=0.5, density=True, label='Simulation')

from scipy.stats import binom

import math as m
def p_d(d, p, N):
    pd = (1-(1 - p**d)**m.comb(N-2, d-1))
    sums = 0
    for i in range(1, d): 
        sums += p_d(i, p, N)
    pd *= (1 - sums)
    return pd

n = 1000
p = 30/(n-1)

data = d2
pd = p_d(2, p, n) + p
x = np.arange(min(data)-n*0.1, max(data)+n*0.1, 1)
y = binom.pmf(x, n, pd)
plt.plot(x, y, color="#648fff", ls = '-', label="$B(N, p+p^{(2)})$")

ax.set_xlabel('$k$')
ax.set_ylabel('$P(k)$', labelpad=-20)

ax.set_xlim(320, 800)

ax.set_ylim(0, 0.05)
ax.set_yticks([0, 0.05])
ax.set_yticklabels(['0', '0.05'])

leg = plt.legend(loc='upper left', fancybox=False, shadow=False, frameon=False, fontsize=13)

ax.text(-0.14, 0.85, "b", fontsize=20, ha='center', va='center', transform=ax.transAxes, color='black', fontweight='bold', fontname='DejaVu Sans')

plt.savefig(f'plots/{file}.pdf')