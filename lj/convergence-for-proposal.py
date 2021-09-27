#!/usr/bin/python3

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import compute
import heat_capacity
import styles
import glob
from mytimer import Timer

E = np.linspace(-133.6, 0, 10000)
E = 0.5*(E[1:] + E[:-1])
dE = E[1] - E[0]
hist = None

lowest_interesting_E = -1.07
highest_interesting_E = -0.5

indices_for_err = np.array([i for i in range(
    len(E)) if lowest_interesting_E <= E[i] <= highest_interesting_E])
E_for_err = E[indices_for_err]


def normalize_S(S):
    S = S - max(S)
    total = np.sum(np.exp(S)*dE)
    return S - np.log(total)


fig, ax = plt.subplots(figsize=[5, 4], num='latest heat capacity')
# [0.005, 0.012, 25, 140])
axins = ax.inset_axes(0.5 * np.array([1, 1, 0.47/0.5, 0.47/0.5]))
# axins.set_xticklabels('')
# axins.set_yticklabels('')

Tmax = 0.6

paths = [# 'wl-tiny-0.0001+0.0001-lnw.dat',
        #  'wl-tiny-0.0001+0.01-lnw.dat',
        #  'wl-tiny-1e-05+0.001-lnw.dat',
        #  'wl-tiny-0.0001+0.001-lnw.dat',
        #  'wl-tiny-1e-05+0.0001-lnw.dat',
        #  'wl-tiny-1e-05+0.01-lnw.dat',
         'z-lj31-lnw.dat',
]

energy_boundaries, mean_e, my_lnw, _, _ = compute.read_file('z-lj31')

# Create a function for the entropy
reference_S, _, _ = compute.linear_entropy(energy_boundaries, mean_e, my_lnw)

fig_S, ax_S = plt.subplots(figsize=[7, 5], num='latest-entropy')
axins_S = ax_S.inset_axes([0.5, 0.1, 0.47, 0.47])
# plt.plot(E, normalize_S(reference_S(E)), ':', label='exact', linewidth=2)
# Cite Calvo, Doye, and Wales "Quantum partition functions from classical distributions: application fo rare-gas clusters"
ax_S.axvline(-133.58642, color='k', linestyle=':')
ax_S.axvline(-133.29382, color='k', linestyle=':')
# ax_S.axvline(-133.10462, color='k', linestyle=':')

axins_S.axvline(-133.58642, color='k', linestyle=':')
axins_S.axvline(-133.29382, color='k', linestyle=':')
# axins_S.axvline(-133.10462, color='k', linestyle=':')

heat_capacity.plot(reference_S, ax=ax, axins=axins, Tmax=Tmax)
ax.indicate_inset_zoom(axins, edgecolor="black", label=None)

ax_S.set_xlim(-133.6, -100)
axins_S.set_xlim(-133.6, -133.1)
ax_S.set_ylim(-500, -75)
axins_S.set_ylim(-500, -325)
ax_S.indicate_inset_zoom(axins_S, edgecolor="black", label=None)

minimum_moves = 1e7

for fname in paths:
    print()
    timer_fname = Timer(f'Everything for {fname}')
    base = fname[:-8]
    method = base[:base.find('-')]

    energy_boundaries, mean_e, my_lnw, my_system, p_exc = compute.read_file(
        base)

    # Create a function for the entropy
    l_function, eee, sss = compute.linear_entropy(
        energy_boundaries, mean_e, my_lnw)
    plt.figure('latest-entropy')
    ax_S.plot(E, normalize_S(l_function(E)), 
                    label=styles.pretty_label(base), marker=styles.marker(base),
                 color=styles.color(base), linestyle=styles.linestyle(base), markevery=100)
    axins_S.plot(E, normalize_S(l_function(E)), 
                    label=styles.pretty_label(base), marker=styles.marker(base),
                 color=styles.color(base), linestyle=styles.linestyle(base), markevery=100)

    plt.figure('latest heat capacity')
    if 'z-' in fname:
        heat_capacity.plot(l_function, fname=fname, ax=ax, axins=axins, Tmax=Tmax)

    try: #### REMOVE *saved.txt FILES IF DATA IS UPDATED ####
        moves = []
        errors_Cv = np.loadtxt(f'{base}-Cv-error-saved.txt')

        for frame_fname in sorted(glob.glob(f'{base}/*-lnw.dat')):
            try:
                frame_moves = int(frame_fname[len(base)+1:-8])
                if frame_moves < minimum_moves:
                    continue
                moves.append(frame_moves)
            except:
                pass
        if len(moves) != len(errors):
            print('Need to recompute, there are new moves!')
            assert(len(moves) == len(errors_Cv))
    except:
        errors_Cv = []
        moves = []
        start_conv = time.process_time()
        for frame_fname in sorted(glob.glob(f'{base}/*-lnw.dat')):
            frame_base = frame_fname[:-8]

            try:
                frame_moves = int(frame_fname[len(base)+1:-8])
                if frame_moves < minimum_moves:
                    continue
                energy_boundaries, mean_e, my_lnw, my_system, p_exc = compute.read_file(
                    frame_base)
                l_function, eee, sss = compute.linear_entropy(
                    energy_boundaries, mean_e, my_lnw)
                moves.append(frame_moves)
                errors_Cv.append(np.abs(heat_capacity.C(
                    heat_capacity.T_peak, reference_S) - heat_capacity.C(heat_capacity.T_peak, l_function)))
            except:
                pass
        print(f'Heat capacity convergence {fname} took %.2g' % (
            time.process_time() - start_conv))
        np.savetxt(f'{base}-Cv-error-saved.txt', errors_Cv)


    plt.figure('convergence_heat_capacity')
    if method in {'wl', 'itwl', 'sad'}:
        plt.loglog(moves, errors_Cv, 
                    label=styles.pretty_label(base), marker=styles.marker(
            base), color=styles.color(base), linestyle=styles.linestyle(base), markevery=2)
    elif method == 'z':
        plt.loglog(moves, errors_Cv, 
                    label=styles.pretty_label(base), color=styles.color(base),
                   linestyle=styles.linestyle(base), linewidth=3)

    del timer_fname
    print()
    # print(base[:base.find('-')]) #for debugging

plt.figure('latest-entropy')
plt.xlabel(r'$E$')
plt.ylabel(r'$S(E)$')
plt.legend()
plt.tight_layout()
plt.savefig('lj31.pdf')

plt.figure('convergence_heat_capacity')
plt.xlabel(r'# of Moves')
plt.ylabel(rf'max error in heat capacity')
plt.ylim(1e-4, 1e3)
plt.xlim(xmin=minimum_moves)
plt.legend()
# make diagonal lines for convergence
x = np.linspace(1e-10, 1e20, 2)
y = 1/np.sqrt(x)
for i in range(20):
    plt.loglog(x, y*10**(4*i/5-2), color='y', alpha=0.5)
plt.savefig('lj31-heat-capacity-convergence.svg')
plt.savefig('lj31-heat-capacity-convergence.pdf')


plt.figure('latest heat capacity')
ax.legend()
plt.xlabel('$T$')
plt.ylabel('heat capacity')
plt.savefig('lj31-heat-capacity.svg')
plt.savefig('lj31-heat-capacity.pdf')

plt.show()
