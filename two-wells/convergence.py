#!/usr/bin/python3

import os
import numpy as np
import matplotlib.pyplot as plt
import system, compute
import glob

equal_E = -.9 # -1.01

E = np.linspace(-system.h_small, 0, 10000)
E = 0.5*(E[1:] + E[:-1])
hist = None

plt.figure('latest-entropy')
plt.plot(E, system.S(E) - system.S(equal_E), ':', label='exact', linewidth=2)

for fname in sorted(glob.glob('*'+system.system+'*.cbor')):
    print(fname)
    base = fname[:-5]

    energy_boundaries, mean_e, my_lnw, my_system, p_exc = compute.read_file(base)
    
    # Create a function for the entropy
    l_function, eee, sss = compute.linear_entropy(energy_boundaries, mean_e, my_lnw)
    plt.figure('latest-entropy')
    plt.plot(E, l_function(E) - l_function(equal_E), label=base)

    plt.figure('fraction-well')
    print('base is', base)
    mean_which = np.loadtxt(f'{base}-which.dat')
    plt.plot(mean_e, mean_which, label=base)

    if os.path.exists(f'{base}-histogram.dat'):
        hist = np.loadtxt(f'{base}-histogram.dat')
        plt.figure('histogram')
        plt.plot(mean_e, hist, label=base)


plt.figure('latest-entropy')
plt.legend()
plt.savefig(system.system+'.svg')

plt.figure('fraction-well')
plt.legend()
plt.savefig(system.system+'-which.svg')

if hist is not None:
    plt.figure('histogram')
    plt.legend()
    plt.savefig(system.system+'-histogram.svg')

plt.show()