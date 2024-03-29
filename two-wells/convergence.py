#!/usr/bin/python3

import os, time
import numpy as np
import matplotlib.pyplot as plt
import system, compute
import heat_capacity, styles
import glob
import time

T = np.linspace(0.001,0.01,175)

E = np.linspace(-system.h_small, 0, 10000)
E = 0.5*(E[1:] + E[:-1])
dE = E[1] - E[0]
hist = None

lowest_interesting_E = -1.07
highest_interesting_E = -0.5

indices_for_err = np.array([i for i in range(len(E)) if lowest_interesting_E <= E[i] <= highest_interesting_E])
E_for_err = E[indices_for_err]

def normalize_S(S):
    S = S - max(S)
    total = np.sum(np.exp(S)*dE)
    return S - np.log(total)

plt.figure('latest-entropy')

correct_S = normalize_S(system.S(E))

correct_S_for_err = correct_S[indices_for_err]
plt.plot(E, correct_S, ':', label='exact', linewidth=2)

fig, ax = plt.subplots(figsize=[5, 4], num='latest heat capacity')
axins = ax.inset_axes( 0.5 * np.array([1, 1, 0.47/0.5, 0.47/0.5]))#[0.005, 0.012, 25, 140])
#axins.set_xticklabels('')
#axins.set_yticklabels('')

heat_capacity.plot(system.S, ax=ax, axins=axins)
ax.indicate_inset_zoom(axins, edgecolor="black")

paths = []
for fname in sorted(glob.glob('*'+system.system+'*-lnw.dat')):
    if not ('half-barrier' in fname):
        if 'sad' in fname:
            if True:#not( '0.01+0.001' in fname):
                paths.append(fname)
        elif 'wl' in fname or 'itwl' in fname:
            if True:#'0.001+' in fname:
                paths.append(fname)
        else:
            paths.append(fname)

for fname in paths:
    print()
    start_fname = time.process_time()
    base = fname[:-8]
    method = base[:base.find('-')]

    energy_boundaries, mean_e, my_lnw, my_system, p_exc = compute.read_file(base)
    
    # Create a function for the entropy
    l_function, eee, sss = compute.linear_entropy(energy_boundaries, mean_e, my_lnw)
    plt.figure('latest-entropy')
    #plt.plot(E, normalize_S(l_function(E)), label=base)
    if method in {'wl','itwl','sad'}:
        plt.plot(E, normalize_S(l_function(E)), label=base, marker = styles.marker(base),
                 color = styles.color(base), linestyle= styles.linestyle(base), markevery=100)
    elif method == 'z':
        plt.plot(E, normalize_S(l_function(E)), label=base, color = styles.color(base), linestyle= styles.linestyle(base))

    plt.figure('latest heat capacity')
    heat_capacity.plot(l_function, fname=fname,ax=ax, axins=axins)
    # correct_Cv = [heat_capacity.C(t,l_function) for t in T]
    # if method in {'wl','itwl','sad'}:
    #     plt.plot(T, correct_Cv, label=base, marker = markers[precision], color = colors[method], linestyle= linestyles[method], markevery=5)
    # elif method == 'z':
    #     plt.plot(T, correct_Cv, label=base, color = colors[method], linestyle= linestyles[method])

    start_well_histogram = time.process_time()
    plt.figure('fraction-well')
    mean_which = np.loadtxt(f'{base}-which.dat')
    plt.plot(mean_e, mean_which, label=base)
    

    if os.path.exists(f'{base}-histogram.dat'):
        hist = np.loadtxt(f'{base}-histogram.dat')
        plt.figure('histogram')
        plt.plot(mean_e, hist, label=base)
    print(f'Plotting fraction and histogram {fname} took %.2g' % (time.process_time() - start_well_histogram))

    try: #### REMOVE *saved.txt FILES IF DATA IS UPDATED ####
        moves = []
        errors = np.loadtxt(f'{base}-S-error-saved.txt')
        errors_Cv = np.loadtxt(f'{base}-Cv-error-saved.txt')

        for frame_fname in sorted(glob.glob(f'{base}/*-lnw.dat')):
            try:
                frame_moves = int(frame_fname[len(base)+1:-8])
                if frame_moves < 1e8:
                    continue
                moves.append(frame_moves)
            except:
                pass
    except:
        errors = []
        errors_Cv = []
        moves = []
        start_conv = time.process_time()
        for frame_fname in sorted(glob.glob(f'{base}/*-lnw.dat')):
            frame_base =frame_fname[:-8]

            try:
                frame_moves = int(frame_fname[len(base)+1:-8])
                if frame_moves < 1e8:
                    continue
                energy_boundaries, mean_e, my_lnw, my_system, p_exc = compute.read_file(frame_base)
                l_function, eee, sss = compute.linear_entropy(energy_boundaries, mean_e, my_lnw)
                moves.append(frame_moves)
                err = np.max(np.abs(normalize_S(l_function(E))[indices_for_err] - correct_S_for_err))
                # print(f'err is {err} for {frame_base}')
                errors.append(err)
                errors_Cv.append(np.abs(heat_capacity.C(heat_capacity.T_peak, system.S) - heat_capacity.C(heat_capacity.T_peak, l_function)))
            except:
                pass
        print(f'Heat capacity convergence {fname} took %.2g' % (time.process_time() - start_conv))
        np.savetxt(f'{base}-S-error-saved.txt', errors)
        np.savetxt(f'{base}-Cv-error-saved.txt', errors_Cv)

    plt.figure('convergence')
    if method in {'wl','itwl','sad'}:
        plt.loglog(moves, errors, label=base, marker = styles.marker(base), color = styles.color(base), linestyle= styles.linestyle(base), markevery=2)
    elif method == 'z':
        plt.loglog(moves, errors, label=base, color = styles.color(base), linestyle= styles.linestyle(base), linewidth = 3)

    plt.figure('convergence_heat_capacity')
    if method in {'wl','itwl','sad'}:
        plt.loglog(moves, errors_Cv, label=base, marker = styles.marker(base), color = styles.color(base), linestyle= styles.linestyle(base), markevery=2)
    elif method == 'z':
        plt.loglog(moves, errors_Cv, label=base, color = styles.color(base), linestyle= styles.linestyle(base), linewidth = 3)

    print(f'Everything for {fname} took %.2g\n' % (time.process_time() - start_fname))
    #print(base[:base.find('-')]) #for debugging

plt.figure('latest-entropy')
plt.xlabel(r'$E$')
plt.ylabel(r'$S(E)$')
plt.legend()
plt.savefig(system.system+'.svg')

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
plt.ylabel(rf'max error in entropy between {lowest_interesting_E} and {highest_interesting_E}')
plt.ylim(1e-2, 1e2)
plt.legend()
#make diagonal lines for convergence
x = np.linspace(1e-10,1e20,2)
y = 1/np.sqrt(x)
for i in range(20):
    plt.loglog(x,y*10**(4*i/5-2), color = 'y',alpha=0.5)
plt.savefig(system.system+'-convergence.svg')
plt.savefig(system.system+'-convergence.pdf')

plt.figure('convergence_heat_capacity')
plt.xlabel(r'# of Moves')
plt.ylabel(rf'max error in entropy between {lowest_interesting_E} and {highest_interesting_E}')
plt.ylim(1e-2, 1e2)
plt.legend()
#make diagonal lines for convergence
x = np.linspace(1e-10,1e20,2)
y = 1/np.sqrt(x)
for i in range(20):
    plt.loglog(x,y*10**(4*i/5-2), color = 'y',alpha=0.5)
plt.savefig(system.system+'-heat-capacity-convergence.svg')
plt.savefig(system.system+'-heat-capacity-convergence.pdf')


plt.figure('latest heat capacity')
ax.legend()
plt.savefig(system.system+'-heat-capacity.svg')
plt.savefig(system.system+'-heat-capacity.pdf')

plt.show()