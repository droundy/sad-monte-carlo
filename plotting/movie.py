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

args = parser.parse_args()

prop_cycle = plt.rcParams['axes.prop_cycle']

colors = cc.glasbey_dark

def beautiful_name(base):
    name = ''
    if 'r-' == base[:2]:
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

    de = abs(np.diff(np.loadtxt(base+'-energy-boundaries.dat')))
    lnw = np.loadtxt(base+'-lnw.dat')[1:-1]
    mean_e = np.loadtxt(base+'-mean-energy.dat')[:-1]
    s_estimate = lnw - np.log(de)
    peak_e = mean_e[np.argmax(s_estimate)]

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

E = np.linspace(mean_e[1:-1].min(), peak_e, 10000)
E = np.linspace(-133.3, -133.0, 10000)

starting_moves = 1e10
for frame in range(len(list(filter(lambda f: parse_moves(f) >= starting_moves, glob.glob(bases[0]+'/*.cbor'))))):
    plt.clf()
    which_color = 0
    plotted_something = False
    plt.ylabel('$S(E)$')

    for base in bases:
        color = colors[which_color]
        which_color += 1

        if base not in moves:
            moves[base] = []
            error[base] = []
        frames = sorted(filter(lambda f: parse_moves(f) >= starting_moves, glob.glob(base+'/*.cbor')))
        if frame >= len(frames):
            continue
        f = os.path.splitext(frames[frame])[0]
        mymove = parse_moves(f)
        moves[base].append(mymove)
        print(f'working on {base} with moves {mymove} which is {f}')

        energy_b = np.loadtxt(f+'-energy-boundaries.dat')
        mean_e = np.loadtxt(f+'-mean-energy.dat')
        my_lnw = np.loadtxt(f+'-lnw.dat')
        
        if energy_b.ndim == 0: #in case of a single value
            energy_b = np.array([energy_b.item()])

        if energy_b[0] < energy_b[-1]:
            energy_b = np.flip(energy_b)
            mean_e = np.flip(mean_e)
            my_lnw = np.flip(my_lnw)

        for me, be in zip(mean_e[-8:], energy_b[-8:]):
            print('energies:', me, be)
        # Create a function for the entropy based on this number of moves:
        l_function, eee, sss = compute.linear_entropy(energy_b, mean_e, my_lnw)
        # l_function, _, _ = compute.step_entropy(energy_b, mean_e, my_lnw)
        entropy_here = l_function(E)
        # plt.plot(E, entropy_here, label=beautiful_name(f))
        plt.plot(eee, sss, '.-', label=beautiful_name(f))
        plt.xlim(E.min(), E.max())
        plt.ylim(entropy_here.min(), entropy_here.max())
        plt.title('$t=%s$' % latex_number(mymove))
        plotted_something = True
    if not plotted_something:
        print('nothing left to plot')
        break
    plt.xlabel('E')
    plt.legend()
    plt.draw_if_interactive()
    plt.pause(0.01)

plt.ioff()
plt.show()
