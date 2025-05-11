import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rc('font', family='Helvetica', size=22)
plt.rc('mathtext', fontset='dejavusans')

file = "fig3-diagram-RR"
df = pd.read_csv(f'scripts/{file}.csv')

a = df['α'].to_numpy()
k = df['k_avg'].to_numpy()
l = df['Δλ_avg'].to_numpy()

a = np.unique(a)

bin_size = 5
k = np.arange(0, 1000 + bin_size, bin_size)
data = np.zeros((len(a), len(k)))
for i in range(len(a)):
    for j in range(len(k)):
            print("i = ", i, "j = ", j, end="\r")
            mask = (df['α'] == a[i]) & (df['k_avg'] >= k[j] - bin_size) & (df['k_avg'] < k[j])
            data[i-1, j-1] = df[mask]['Δλ_avg'].mean()

fig, ax = plt.subplots(1, 1, figsize=(7, 6))

plt.subplots_adjust(wspace=0.5, hspace=2, left=0.13, bottom=0.12, right=0.89, top=0.92)
ax.tick_params(direction='out', top=True, right=True)

N = 1000
p = k / (N - 1)

# --- heatmap magma inferno
c = ax.pcolormesh(p, a, data, cmap='inferno', linewidth=0, rasterized=True, vmin=0, vmax=1)

cbar = plt.colorbar(c, fraction=0.046, pad=0.04, ticks=[0, 0.3, 0.7, 1])
cbar.ax.set_yticklabels(['0', '0.3', '0.7', '1'])

cbar.ax.set_ylabel('$\\langle \\zeta \\rangle$', rotation=0, fontsize=22, labelpad=2, y=0.55)

# --- heatmap countour
from scipy.ndimage import gaussian_filter
data2 = gaussian_filter(data, sigma=0.3) # smoothing the heatmap countour
levels = [0.1, 0.3, 0.5, 0.7, 0.9]
colors = plt.cm.magma(np.linspace(0, 1, len(levels)))
c = ax.contour(p, a, data2, levels=levels, colors="white", linestyles='-', linewidths=2)

# --- theoretical contour
levels = [0.1, 0.3, 0.5, 0.7, 0.9]

#---
from math import comb

def p_d_memoized(N):

    memo = {}
    def p_d(d, p):
        if (d, p) in memo:
            return memo[(d, p)]
        if d == 1:
            pd = p
        else:
            pd = (1 - (1 - p**d)**comb(N-2, d-1))
            sums = sum(p_d(i, p) for i in range(1, d))
            pd *= (1 - sums)
        memo[(d, p)] = pd
        return pd
    return p_d

def Δλ_theory(N, k, a, D):
    p = k / (N - 1)
    p_d_func = p_d_memoized(N)
    Δλ = sum(d**(-a) * p_d_func(d, p) for d in range(2, D))
    return Δλ
#---

k_theory = np.linspace(0, 1000, 1000)
theory_data = np.zeros((len(a), len(k_theory)))
for i in range(len(a)):
    for j in range(len(k_theory)):
        theory_data[i, j] = Δλ_theory(N, k_theory[j], a[i], D=100)
    print("step: ", 100*i/len(a), "%", end="\r")

colors = plt.cm.magma(np.linspace(0, 1, len(levels)))

p_theory = k_theory / (N - 1)

c = ax.contour(p_theory, a, theory_data, levels=levels, colors="white", linestyles='--', linewidths=2)

# countour labels
ax.text(0.14, 3.2, '$0.1$', fontsize=15, rotation=-20, ha='center', va='center', color='white')
ax.text(0.13, 1.64, '$0.3$', fontsize=15, rotation=-20, ha='center', va='center', color='white')
ax.text(0.12, 0.92, '$0.5$', fontsize=15, rotation=-20, ha='center', va='center', color='white')
ax.text(0.11, 0.47, '$0.7$', fontsize=15, rotation=-20, ha='center', va='center', color='white')
ax.text(0.1, 0.1, '$0.9$', fontsize=15, rotation=-20, ha='center', va='center', color='white')

cbar.add_lines(c)
# ---

ax.set_ylabel('$\\alpha$')
ax.set_xlabel('$p$', labelpad=-10)

xticks = [0, 0.3, 0.7, 1]
yticks = [0,1,2,3,4]

ax.set_xticks(xticks)
ax.set_yticks(yticks)

ax.set_xticklabels(["0", "0.3", "0.7", "1"])

ax.set_xlim(xticks[0], xticks[-1])
ax.set_ylim(yticks[0], yticks[-1])

ax.text(-0.14, 0.85, "b", fontsize=25, ha='center', va='center', transform=ax.transAxes, color='black', fontweight='bold', fontname='DejaVu Sans')

plt.savefig(f'plots/{file}.pdf')