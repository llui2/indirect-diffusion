import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb  # Use scipy for comb, as it's optimized
from functools import lru_cache  # For memoization

plt.rc('font', family='Times', size=20)
plt.rc('mathtext', fontset='cm')

fig, ax = plt.subplots(figsize=(5, 5))
plt.subplots_adjust(wspace=0, hspace=0, left=0.12, top=0.95, right=0.95, bottom=0.12)

ax.tick_params(direction='in', top=True, right=True)
ax.tick_params(which='minor', direction='in', top=True, right=True)

# Memoized version of p_d using LRU cache
@lru_cache(maxsize=None)
def p_d(d, p_tuple, N):
    p = np.array(p_tuple)  # Convert back to numpy array inside the function
    if d == 1:
        return 1 - (1 - p**d)**comb(N-2, d-1)
    pd = (1 - (1 - p**d)**comb(N-2, d-1))
    pd *= (1 - sum(p_d(i, p_tuple, N) for i in range(1, d)))
    return pd

Ns = [1e3, 1e4, 1e5, 1e6, 1e7]
Ns = [1e2, 1e3, 1e4, 1e5, 1e6]
Ns = [int(N) for N in Ns]

# Define color palette
colors = [(244, 28, 84), (255, 159, 0), (6, 200, 100), (0, 170, 240), (50, 80, 250)]
colors = [tuple(c / 255 for c in color) for color in colors[::-1]]

# Generate p array
p = np.linspace(0, 0.05, 10000)

for N_idx, N in enumerate(Ns):
    dmax = 5
    pc = 2*N**(-(dmax - 1) / dmax)
    
    # Convert p to a tuple for caching
    p_tuple = tuple(p)
    
    # Compute p_d values using list comprehension for better readability
    pds = [p_d(d, p_tuple, N) for d in range(1, dmax + 1)]
    
    nu = np.sum(pds, axis=0)  # Sum along the first axis
    
    ax.plot(p - pc, nu, label=f'$N = 10^{{{int(np.log10(N))}}}$', lw=2, color=colors[N_idx])

ax.set_ylabel('$\\nu$', fontfamily='Times', labelpad=-15)
ax.set_xlabel('$p - p_c$', fontfamily='Times', labelpad=-4)

ax.set_xlim(-0.03, 0.03)
ax.set_xticks([-0.02, 0, 0.02])
ax.set_xticklabels(["-0.02", "0", "0.02"])

ax.set_ylim(0, 1)
ax.set_yticks([0, 0.3, 0.7, 1])
ax.set_yticklabels(["0", "0.3", "0.7", "1"])

ax.xaxis.set_tick_params(pad=5)
ax.legend(loc='upper left', fontsize=15, frameon=False)

plt.savefig('plots/fig-phasetransition-D.pdf')
plt.close()
