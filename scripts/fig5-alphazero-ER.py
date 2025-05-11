import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rc('font', family='Helvetica', size=18)
plt.rc('mathtext', fontset='dejavusans')

def lambda2_ER(N, p, l):
    p = l * p
    gamma = 0.57721566490153286060651209008240243104215933593992
    lambda2s = p * (N - 1) - np.sqrt(2 * p * (1 - p) * (N - 1) * np.log(N)) + np.sqrt((N - 1) * p * (1 - p) / (2 * np.log(N))) * np.log(np.sqrt(2 * np.pi * np.log(N**2 / (2 * np.pi)))) - np.sqrt((N - 1) * p * (1 - p) / (2 * np.log(N))) * gamma
    return lambda2s / l

def lambda2_RR(N, k, l):
    k = p * (N-1) * l
    l_t = k - 2 * np.sqrt(k - 1)
    return l_t / l

#----
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

def p_total(d, p, N):
    p_d_func = p_d_memoized(N)
    return p + sum(p_d_func(i, p) for i in range(2, d + 1))
#----

def f(x, d, N):
    return (lambda2_ER(N, p_total(d,x,N), 1) - lambda2_ER(N, x, 1)) / N

N = 1000

fig = plt.figure(figsize=(5, 5))
ax = plt.subplot()

plt.subplots_adjust(wspace=0, hspace=0, left=0.15,
                    top=0.95, right=0.95, bottom=0.13)
ax.tick_params(direction='in', top=False, right=False)
ax.tick_params(which='minor', direction='in', top=False, right=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

file = "fig5-alphazero-ER"

df = pd.read_csv(f'scripts/{file}-2.csv')
p = df['k_avg'].to_numpy() / (N - 1)
t = df['Δλ_avg'].to_numpy()
e = df['σ_Δλ'].to_numpy()

df2 = pd.read_csv(f'scripts/{file}-D.csv')
p2 = df2['k_avg'].to_numpy() / (N - 1)
t2 = df2['Δλ_avg'].to_numpy()
e2 = df2['σ_Δλ'].to_numpy()


# Simulation
# ax.plot(p, t, color="red", marker='o', markersize=1,
#         linewidth=1)
# ax.fill_between(p, t + e, t - e, color="red", alpha=0.2)

ax.errorbar(p, t, yerr=e, color="#fe6202", fmt='', ls="", markersize=1, capsize=3, capthick=1.5, elinewidth=1.5, clip_on=False)

# ax.plot(p2, t2, color="blue", marker='o', markersize=1,
#         linewidth=1)
# ax.fill_between(p2, t2 + e2, t2 - e2, color="blue", alpha=0.2)

ax.errorbar(p2, t2, yerr=e2, color="#785ef0", fmt='', ls="", markersize=1, capsize=3, capthick=1.5, elinewidth=1.5, clip_on=False)


p_theory = np.logspace(-3, 0, num=100, endpoint=True, base=10.0)

t_theory = [f(x, 2, N) for x in p_theory]
ax.plot(p_theory, t_theory, color="#fe6202", linestyle='-', linewidth=1.3, label="$d_{\\max} = 2$")

t_theory = [(N - lambda2_ER(N, x, 1)) / N for x in p_theory]
ax.plot(p_theory, t_theory, color="#785ef0", linestyle='-', linewidth=1.3, label="$d_{\\max} = D$")

x = np.log(N)/N
# ax.plot([x, x], [-1, 2], color='gray', linestyle='-', linewidth=1)
ax.fill_between([0, x], -1, 2, color='gray', alpha=0.2)


ax.set_ylabel('$\\langle \\zeta \\rangle$', labelpad=-5)
ax.set_xlabel('$p$', labelpad=-2)

ax.set_xscale('log')
ax.set_xlim(0.005, 1.01)
ax.set_ylim(-0.03, 1.03)

ax.xaxis.set_tick_params(pad=5)

# ax.set_xticks([0, 0.3, 0.7, 1])
ax.set_yticks([0, 0.3, 0.7, 1])

# ax.set_xticklabels(["0", "0.3", "0.7", "1"])
ax.set_yticklabels(["0", "0.3", "0.7", "1"])

leg = plt.legend(loc=(0.5,0), fancybox=False, shadow=False, frameon=False, ncol=1, fontsize="15", labelspacing=.6)
plt.setp(leg.get_lines(), linewidth=1.2)

ax.text(-0.14, 0.85, "a", fontsize=22, ha='center', va='center', transform=ax.transAxes, color='black', fontweight='bold', fontname='DejaVu Sans')

plt.savefig(f'plots/{file}.pdf')
plt.close()
