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
        base = base[2:]
    elif 'wl-' == base[:3]:
        name += 'WL '
        base = base[3:]
    elif 'sad-' == base[:4]:
        name += 'SAD '
        base = base[4:]
    elif 'itwl-' == base[:5]:
        name += r'$t^{-1}$-WL '
        base = base[5:]
    if base == 'erfinv' or base == 'quadratic':
        name += ''
        base = ''
    elif base[:7] == 'erfinv-':
        base = base[7:]
        name += rf' $\Delta E = {base}$'
        base = ''
    elif base[:10] == 'quadratic-':
        base = base[10:]
        name += rf' $\Delta E = {base}$'
        base = ''
    return name + base

#Read Data
moves = {}
error = {} #store the max error in each move
bases = []
print('base', args.base)

energy_boundaries = {}
entropy_boundaries = {}
mean_energy = {}
lnw = {}
systems = {}
#each file has different path (including extension) so concatenating is easy
for base in args.base:
    #change base to have the cbor files. currently has the directory
    if '.cbor' in base or '.yaml' in base:
        base = base[:-5]
    bases.append(base)

for base in bases:
    print('reading', base)
    with open(base+'-system.dat') as f:
        systems[base] = yaml.safe_load(f)

    de = abs(np.diff(np.loadtxt(base+'-energy-boundaries.dat')))
    lnw = np.loadtxt(base+'-lnw.dat')[1:-1]
    mean_e = np.loadtxt(base+'-mean-energy.dat')[:-1]
    s_estimate = lnw - np.log(de)
    peak_e = mean_e[np.argmax(s_estimate)]

print('done reading bases')
sigma = 1

E = np.linspace(mean_e.min(), peak_e, 1000)

for frame in range(len(glob.glob(bases[0]+'/*-lnw.dat'))):
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
        frames = sorted(glob.glob(base+'/*.cbor'))
        if frame >= len(frames):
            continue
        f = os.path.splitext(frames[frame])[0]
        mymove = float(os.path.basename(f))
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

        # Create a function for the entropy based on this number of moves:
        l_function, _, _ = compute.linear_entropy(energy_b, mean_e, my_lnw)
        # l_function, _, _ = compute.step_entropy(energy_b, mean_e, my_lnw)
        entropy_here = l_function(E)
        plt.plot(E, entropy_here, label=beautiful_name(f))
        plotted_something = True
    if not plotted_something:
        print('nothing left to plot')
        break
    plt.xlabel('E')
    plt.legend()
    plt.draw_if_interactive()
    plt.pause(0.1)

plt.ioff()
plt.show()
