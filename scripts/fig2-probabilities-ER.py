import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as m

plt.rc('font', family='Helvetica', size=18)
plt.rc('mathtext', fontset='dejavusans')

fig = plt.figure(figsize=(5, 5))
ax = plt.subplot()

plt.subplots_adjust(wspace=0, hspace=0, left=0.15,
                    top=0.95, right=0.95, bottom=0.13)
ax.tick_params(direction='in', top=False, right=False)
ax.tick_params(which='minor', direction='in', top=False, right=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

N = 1000
dmax = 5

color = ["#FF0028", "#248FFF", "#FFB000", "#785EF0", "#21C215"]
color = ["#3f5aae", "#97d8c4", "#6b9ac4", "#f4ba42"]

thickness = 1.7

x = np.log(N)/N
# ax.plot([x, x], [0, 1], color='gray', linestyle='-', linewidth=thickness)
ax.fill_between([0, x], 0, 1, color='gray', alpha=0.2)


# ------------------------------------------
def p_k(k, p, N):
    pk = (1-(1 - p**k)**m.comb(N-2, k-1))
    sums = 0
    for i in range(1, k): 
        sums += p_k(i, p, N)
    pk *= (1 - sums)
    return pk

p = np.linspace(0, 1, 10000)
for k in range(2, dmax+1):
    pk = p_k(k, p, N)
    ax.plot(p, pk, color=color[k-2], linestyle='-', label=f'$d = {k}$', linewidth=thickness)

# ------------------------------------------

file = "fig2-probabilities-ER"
df = pd.read_csv(f'scripts/{file}.csv')
p1 = df['p1']

for i in range(2,dmax+1):
    ax.errorbar(p1, df[f'p{i}'], yerr=df[f'sd_p{i}'], color=color[i-2], fmt=' ', capsize=3, capthick=thickness, elinewidth=thickness)  

# ------------------------------------------

ax.set_ylabel('$p^{(d)}$', labelpad=-10)
ax.set_xlabel('$p$', labelpad=-2)

ax.set_xscale('log')

ax.set_xlim(0.005, 1)
ax.set_ylim(0, 1)

# ax.set_xticks([0, 0.1])
ax.set_yticks([0, 0.3, 0.7, 1])

ax.xaxis.set_tick_params(pad=5)
# ax.set_xticklabels(["0", "0.1"])
ax.set_yticklabels(["0", "0.3", "0.7", "1"])

ax.legend(loc=(0.55,0.1), fontsize=15, frameon=False)

ax.text(-0.14, 0.85, "a", fontsize=22, ha='center', va='center', transform=ax.transAxes, color='black', fontweight='bold', fontname='DejaVu Sans')

plt.savefig(f'plots/{file}.pdf')
plt.close()
