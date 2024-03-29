#!/usr/bin/python3

import numpy as np
import yaml, cbor, argparse, sys, os, glob
import scipy.constants as const
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import colorcet as cc

import compute

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('base', nargs='*', help = 'the yaml or cbor files')
parser.add_argument('--intensive', action='store_true')
parser.add_argument('--T-max', action='store', type=float,
                    help = "the highest energy for which we want heat capacity")

args = parser.parse_args()

prop_cycle = plt.rcParams['axes.prop_cycle']

colors = cc.glasbey_dark

def beautiful_name(base):
    name = ''
    if 'r-' == base[:2]:
        name += "Zeno's "
        base = base[2:].split('/')[0]
    elif 'z-' == base[:2]:
        name += "Zeno's "
        base = base[2:].split('/')[0]
    elif 'wl-' == base[:3]:
        name += 'WL '
        base = base[3:].split('/')[0]
    elif 'sad-' == base[:4]:
        name += 'SAD '
        base = base[4:].split('/')[0]
    elif 'itwl-' == base[:5]:
        name += r'$t^{-1}$-WL '
        base = base[5:].split('/')[0]
    if base == 'erfinv' or base == 'quadratic':
        name += ''
        base = ''
    elif base[:7] == 'erfinv-':
        base = base[7:].split('/')[0]
        name += rf' $\Delta E = {base}$'
        base = ''
    elif base[:10] == 'quadratic-':
        base = base[10:].split('/')[0]
        name += rf' $\Delta E = {base}$'
        base = ''
    return name + base

def heat_capacity(T, S_func):
    C = np.zeros_like(T)

    boltz_arg_boundaries = S_func(energy_boundaries) - energy_boundaries/T[-1]
    P_boundaries = np.exp(boltz_arg_boundaries - boltz_arg_boundaries.max())

    if (P_boundaries > 0).any():
        E = np.linspace(6*mean_e[1:-1].min() - 5*energy_boundaries[-1], energy_boundaries[P_boundaries > 0].max(), 10000)
        S = S_func(E)
        okay_indices = S==S
        S = S[okay_indices]
        E = E[okay_indices]
        for i in range(len(T)):
            boltz_arg = S - E/T[i]
            P = np.exp(boltz_arg - boltz_arg.max())
            P = P/P.sum()
            U = (E*P).sum()
            C[i] = ((E-U)**2*P).sum()/T[i]**2
    if args.intensive:
        C /= N
    return C

def canonical(T, E, X_of_E, S_func):
    X_of_T = np.zeros_like(T)
    X_of_E = 1.0*X_of_E
    badnum = np.isnan(X_of_E) | np.isinf(X_of_E)
    X_of_E[badnum] = 0

    S = S_func(E)
    for i in range(len(T)):
        boltz_arg = S - E/T[i]
        P = np.exp(boltz_arg - boltz_arg.max())
        P[badnum] = 0
        P = P/P.sum()
        X_of_T[i] = (X_of_E*P).sum()
    return X_of_T

#Read Data
moves = {}
error = {} #store the max error in each move
bases = []
print('base', args.base)

#each file has different path (including extension) so concatenating is easy
for base in args.base:
    #change base to have the cbor files. currently has the directory
    if '.cbor' in base or '.yaml' in base:
        base = base[:-5]
    bases.append(base)

for base in bases:
    print('reading', base)

    energy_boundaries = np.loadtxt(base+'-energy-boundaries.dat')
    de = abs(np.diff(energy_boundaries))
    lnw = np.loadtxt(base+'-lnw.dat')
    mean_e = np.loadtxt(base+'-mean-energy.dat')
    s_estimate = lnw[1:-1] - np.log(de)
    peak_e = mean_e[np.argmax(s_estimate)]

    if energy_boundaries[0] < energy_boundaries[-1]:
        energy_boundaries = np.flip(energy_boundaries)
        mean_e = np.flip(mean_e)
        lnw = np.flip(lnw)
    reference_function, reference_energy, reference_entropy = compute.linear_entropy(energy_boundaries, mean_e, lnw)
    unscaled_reference_function = reference_function
    if args.intensive:
        energy_boundaries, mean_e, _, system, p_exc_reference = compute.read_file(base)
        reference_N = system['N']
        N = reference_N
        reference_energy /= reference_N
        reference_entropy /= reference_N
        # energy_boundaries /= reference_N
        # lnw /= reference_N
        # peak_e /= reference_N
        # s_estimate /= reference_N
        # de /= reference_N
        # mean_e /= reference_N
        reference_function = lambda e: unscaled_reference_function(reference_N*e)/reference_N
        pressure_reference, T_reference = compute.pressure_temperature(system['density'], energy_boundaries, mean_e, p_exc_reference)
        mean_e_reference = mean_e
        print(pressure_reference)
            
        T_to_save = np.arange(0.5, 100, 0.5)
        pressure_to_save = canonical(T_to_save, mean_e, pressure_reference, unscaled_reference_function)
        np.savetxt(base+'-p-vs-T.dat',
            np.transpose(np.vstack([T_to_save, pressure_to_save])),
            fmt='%.3g')

print('done reading bases')
sigma = 1

def parse_moves(name):
    return float(os.path.basename(os.path.splitext(name)[0]))
def latex_number(x):
    if x > 10000:
        s = '%.3e' % x
        mantissa, exponent = s.split('e+')
        if mantissa == '1' or mantissa == '1.000':
            return rf'10^{{{exponent}}}'
        else:
            return rf'{mantissa}\times 10^{{{exponent}}}'
    return '%.3g' % x

plot_Cv = args.T_max is not None

E = np.linspace(mean_e[1:-1].min(), min(peak_e, energy_boundaries.max()), 10000)
# E = np.linspace(-133.59, -133.0, 10000)
# E = np.linspace(-133.59, 0, 10000)

if plot_Cv:
    T_min = 2.0/(unscaled_reference_function(mean_e[1:-1].min()) - unscaled_reference_function(mean_e[1:-1].min()-1))
    T = np.linspace(T_min, args.T_max, 1000)

    boltz_arg_boundaries = unscaled_reference_function(energy_boundaries) - energy_boundaries/T[-1]
    P_boundaries = np.exp(boltz_arg_boundaries - boltz_arg_boundaries.max())
    E = np.linspace(mean_e[1:-1].min(), energy_boundaries[P_boundaries > 0.1].max(), 1000)

if args.intensive:
    E /= reference_N

if plot_Cv:
    Cv_reference = heat_capacity(T, unscaled_reference_function)

use_inset = False

starting_moves = 1e8
for frame in range(len(list(filter(lambda f: parse_moves(f) >= starting_moves, glob.glob(os.path.join(base,'*.cbor')))))):
    which_color = 0
    plotted_something = False
    if plot_Cv:
        plt.figure('$C_V$', figsize=(8,6))
        plt.clf()
        plt.ylabel('$C_V(T)$')
        plt.xlabel('$T$')
        ax = plt.gca()
        ax.set_xlim(0, max(T))
        ax.set_ylim(0, 1.3*max(Cv_reference))
        if use_inset:
            axins = ax.inset_axes([0.1, 0.5, 0.4, 0.47])
            axins.set_xlim(T.min(), 0.01)
            # axins.set_ylim(0, 1.3*np.max(Cv_reference[T<0.01]))
        
        ax.plot(T, Cv_reference, ':', color='gray', label='reference')
        if use_inset:
            axins.plot(T, Cv_reference, ':', color='gray', label='reference')

    if args.intensive:
        plt.figure('pressure', figsize=(8,6))
        plt.clf()
        plt.ylabel('$p$')
        plt.xlabel('$T$')
        # pmax = 0
        pressure_of_T = canonical(T, mean_e_reference, pressure_reference, unscaled_reference_function)
        plt.plot(T, pressure_of_T, ':', color='gray', label='reference')
    plt.figure('S(E)', figsize=(8,6))
    plt.clf()
    plt.ylabel('$S(E)$')

    ymin, ymax = np.Infinity, -np.Infinity
    plt.plot(reference_energy, reference_entropy, ':', color='gray', label='reference')

    for base in bases:
        color = colors[which_color]
        which_color += 1

        if base not in moves:
            moves[base] = []
            error[base] = []
        frames = sorted(filter(lambda f: parse_moves(f) >= starting_moves, glob.glob(os.path.join(base,'*.cbor'))))
        if frame >= len(frames):
            continue
        f = os.path.splitext(frames[frame])[0]
        mymove = parse_moves(f)
        moves[base].append(mymove)
        print(f'working on {base} with moves {mymove} which is {f}')

        energy_boundaries, mean_e, my_lnw, system, p_exc = compute.read_file(f)
        
        # Create a function for the entropy based on this number of moves:
        l_function, eee, sss = compute.linear_entropy(energy_boundaries, mean_e, my_lnw)
        if args.intensive:
            N = system['N']
            eee /= N
            sss /= N
            # function_raw = l_function
            # l_function = lambda e: function_raw(N*e)/N

            # energy_boundaries /= N
            # mean_e /= N
            # my_lnw /= N
            plt.figure('pressure')
            density = system['density']
            pressure, T_here = compute.pressure_temperature(density, energy_boundaries, mean_e, p_exc)
            # plt.plot(mean_e/N, pressure, '.-', label=beautiful_name(f))
            # plt.xlim(E.min(), E.max())
            # if len(pressure[mean_e/N <= E.max()] > 0):
            #     pmax = max(pmax, pressure[mean_e/N <= E.max()].max())

            pressure_of_T = canonical(T, mean_e, pressure, l_function)
            plt.plot(T, pressure_of_T, '-', label=beautiful_name(f))

            plt.plot(T_here[T_here<=T.max()], pressure[T_here<=T.max()], '.-', label=beautiful_name(f))
            plt.xlim(0, T.max())

        entropy_here = l_function(E)
        offset = 0 # l_function(-133)
        plt.figure('S(E)')
        # plt.plot(E, entropy_here, label=beautiful_name(f))
        plt.plot(eee, sss - offset, '.-', label=beautiful_name(f), markersize=4)
        plt.xlim(E.min(), E.max())
        if not np.any(np.isnan(entropy_here)) and np.any(eee <= E.max()) and np.any(eee >= E.min()):
            ymin = min(ymin, sss[eee >= E.min()].min() - offset)
            ymax = max(ymax, sss[eee <= E.max()].max() - offset)
        #     plt.ylim(entropy_here.min() - offset, entropy_here.max() - offset)
        plt.title('$t=%s$' % latex_number(mymove))

        if plot_Cv:
            plt.figure('$C_V$')
            Cv= heat_capacity(T, l_function)
            ax.plot(T, Cv, '-', label=beautiful_name(f), markersize=4)
            if use_inset:
                axins.plot(T, Cv, '-', label=beautiful_name(f), markersize=4)
                ax.indicate_inset_zoom(axins, edgecolor='k')

            diff_by = 10
            # pressure, T_here = compute.pressure_temperature(density, energy_boundaries, mean_e, p_exc)
            # T_mean = 0.5*(T_here[diff_by:] + T_here[:-diff_by])
            # Cv_canonical = (mean_e[diff_by:]-mean_e[:-diff_by])/(T_here[diff_by:]-T_here[:-diff_by])
            # if args.intensive:
            #     Cv_canonical /= N
            
            # T_here = (eee[diff_by:]-eee[:-diff_by])/(sss[diff_by:]-sss[:-diff_by])
            # e_here = 0.5*(eee[diff_by:]+eee[:-diff_by])

            # Cv_canonical = (e_here[diff_by:]-e_here[:-diff_by])/(T_here[diff_by:]-T_here[:-diff_by])
            # T_mean = 0.5*(T_here[diff_by:] + T_here[:-diff_by])
            # ax.plot(T_mean, Cv_canonical, '.-', label=beautiful_name(f))
        plotted_something = True
    if not plotted_something:
        print('nothing left to plot')
        break
    plt.figure('S(E)')
    if np.isfinite(ymin) and np.isfinite(ymax):
        plt.ylim(ymin, ymax)
    plt.xlabel('E')
    if 'lj31' in bases[0]:
        # Cite Calvo, Doye, and Wales "Quantum partition functions from classical distributions: application fo rare-gas clusters"
        plt.axvline(-133.58642, color='k', linestyle=':')
        plt.axvline(-133.29382, color='k', linestyle=':')
        plt.axvline(-133.10462, color='k', linestyle=':')
    plt.legend()
    if args.intensive:
        plt.figure('pressure')
        plt.title('$t=%s$' % latex_number(mymove))
        # if np.isfinite(pmax):
        #     plt.ylim(0, pmax)
        plt.legend()
    if plot_Cv:
        plt.figure('$C_V$')
        plt.legend()
    plt.draw_if_interactive()
    plt.pause(1.01)

plt.ioff()
plt.show()
