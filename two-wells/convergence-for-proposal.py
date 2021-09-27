#!/usr/bin/python3

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import system
import compute
import heat_capacity
import styles
import glob
from mytimer import Timer

T = np.linspace(0.001, 0.01, 175)

max_interesting_E = -1.442360888736597957e-01
E = np.linspace(-system.h_small, max_interesting_E, 100000)
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


plt.figure('latest-entropy')

correct_S = normalize_S(system.S(E))

correct_S_for_err = correct_S[indices_for_err]
plt.plot(E, correct_S, ':', label='exact', linewidth=2)
plt.xlim(-system.h_small*1.005, max_interesting_E)

fig, ax = plt.subplots(figsize=[5, 4], num='latest heat capacity')
# [0.005, 0.012, 25, 140])
axins = ax.inset_axes(np.array([0.27, 0.27, 0.7, 0.7]))
# axins.set_xticklabels('')
# axins.set_yticklabels('')
Tmax = 0.25
ax.set_ylim(0, 119)
ax.set_xlim(0, Tmax)

paths = ['wl-tiny-0.0001+0.001-lnw.dat',
        #  'wl-tiny-1e-05+0.001-lnw.dat',
         'itwl-tiny-0.0001+0.001-lnw.dat',
         'sad-tiny-0.0001+0.001-lnw.dat',
        #  'wl-tiny-1e-05+0.0001-lnw.dat',
        #  'wl-tiny-1e-05+0.01-lnw.dat',
         'z-i16-tiny-lnw.dat',
]

minimum_moves = 1e6

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
    #plt.plot(E, normalize_S(l_function(E)), label=base)
    if method in {'wl', 'itwl', 'sad'}:
        plt.plot(E, normalize_S(l_function(E)), 
                    label=styles.pretty_label(base), marker=styles.marker(base),
                 color=styles.color(base), linestyle=styles.linestyle(base), markevery=100)
    elif method == 'z':
        for e in energy_boundaries:
            plt.axvline(e)
        plt.plot(E, normalize_S(l_function(E)),
                    label=styles.pretty_label(base),
                 color=styles.color(base), linestyle=styles.linestyle(base))

    plt.figure('latest heat capacity')
    if 'z-' in fname:
        heat_capacity.plot(l_function, fname=fname, ax=ax, axins=axins, Tmax=Tmax)
    # correct_Cv = [heat_capacity.C(t,l_function) for t in T]
    # if method in {'wl','itwl','sad'}:
    #     plt.plot(T, correct_Cv, label=base, marker = markers[precision], color = colors[method], linestyle= linestyles[method], markevery=5)
    # elif method == 'z':
    #     plt.plot(T, correct_Cv, label=base, color = colors[method], linestyle= linestyles[method])

    start_well_histogram = time.process_time()
    plt.figure('fraction-well')
    mean_which = np.loadtxt(f'{base}-which.dat')
    plt.plot(mean_e, mean_which,
                    label=styles.pretty_label(base), )

    if os.path.exists(f'{base}-histogram.dat'):
        hist = np.loadtxt(f'{base}-histogram.dat')
        plt.figure('histogram')
        plt.plot(mean_e, hist, 
                    label=styles.pretty_label(base), )
    print(f'Plotting fraction and histogram {fname} took %.2g' % (
        time.process_time() - start_well_histogram))

    start_conv = time.process_time()
    try:
        errors = np.loadtxt(f'{base}-errors-saved.txt')
        errors_Cv = np.loadtxt(f'{base}-errors-Cv-saved.txt')
        moves = np.loadtxt(f'{base}-moves-saved.txt')
    except:
        errors = []
        errors_Cv = []
        moves = []
        for frame_fname in sorted(glob.glob(f'{base}/*-lnw.dat'))[::2]:
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
                err = np.max(np.abs(normalize_S(l_function(E))[
                            indices_for_err] - correct_S_for_err))
                # print(f'err is {err} for {frame_base}')
                errors.append(err)
                errors_Cv.append(np.abs(heat_capacity.C(
                    heat_capacity.T_peak, system.S) - heat_capacity.C(heat_capacity.T_peak, l_function)))
            except:
                pass
        np.savetxt(f'{base}-errors-saved.txt', errors)
        np.savetxt(f'{base}-errors-Cv-saved.txt', errors_Cv)
        np.savetxt(f'{base}-moves-saved.txt', moves)
    print(f'Heat capacity convergence {fname} took %.2g' % (
        time.process_time() - start_conv))

    plt.figure('convergence')
    plt.loglog(moves, errors, 
               label=styles.pretty_label(base), marker=styles.marker(base),
               color=styles.color(base), linestyle=styles.linestyle(base), markevery=2)

    plt.figure('convergence_heat_capacity')
    plt.loglog(moves, errors_Cv, 
               label=styles.pretty_label(base), marker=styles.marker(base),
               color=styles.color(base), linestyle=styles.linestyle(base), markevery=2)

    del timer_fname
    print()
    # print(base[:base.find('-')]) #for debugging

plt.figure('latest-entropy')
plt.xlabel(r'$E$')
plt.ylabel(r'$S(E)$')
plt.legend()
plt.savefig(system.system+'-entropy.pdf')

plt.figure('fraction-well')
plt.xlabel(r'E')
plt.ylabel(r'Proportion in Small Well')
plt.legend()
plt.savefig(system.system+'-which.svg')

if hist is not None:
    plt.figure('histogram')
    plt.xlabel(r'$E$')
    plt.ylabel(r'# of Visitors')
    plt.legend()
    plt.savefig(system.system+'-histogram.svg')

plt.figure('convergence')
plt.xlabel(r'# of Moves')
plt.ylabel(
    rf'max error in entropy between {lowest_interesting_E} and {highest_interesting_E}')
plt.ylim(1e-2, 1e2)
plt.legend()
# make diagonal lines for convergence
x = np.linspace(1e-10, 1e20, 2)
y = 1/np.sqrt(x)
for i in range(40):
    plt.loglog(x, y*10**(4*i/5-3), color='g', alpha=0.25)
plt.savefig(system.system+'-convergence.svg')
plt.savefig(system.system+'-convergence.pdf')

plt.figure('convergence_heat_capacity')
plt.xlabel(r'# of Moves')
plt.ylabel(rf'max error in heat capacity')
heat_capacity_convergence_minimum = 1e-3
plt.ylim(heat_capacity_convergence_minimum, 2e2)
plt.xlim(xmin=minimum_moves, xmax=1e12)
plt.legend()
# make diagonal lines for convergence
x = np.linspace(1e-10, 1e20, 2)
for i in range(1, 40, 2):
    plt.loglog(x, heat_capacity_convergence_minimum*np.sqrt(heat_capacity_convergence_minimum*10**i/x), color='xkcd:gray',
               alpha=0.25, linewidth=0.5)
plt.savefig(system.system+'-heat-capacity-convergence.svg')
plt.savefig(system.system+'-heat-capacity-convergence.pdf')


plt.figure('latest heat capacity')
heat_capacity.plot(system.S, ax=ax, axins=axins, Tmax=Tmax)

ax.legend()
plt.xlabel('$T$')
# axins.set_xlabel('$T$')
plt.ylabel('heat capacity')
plt.savefig(system.system+'-heat-capacity.svg')
plt.savefig(system.system+'-heat-capacity.pdf')

plt.show()
