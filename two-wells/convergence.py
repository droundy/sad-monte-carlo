#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import system, compute
import glob

equal_E = -1.0

E = np.linspace(-system.h_small, 0, 10000)
E = 0.5*(E[1:] + E[:-1])
plt.plot(E, system.S(E) - system.S(equal_E), ':', label='exact', linewidth=2)

for fname in sorted(glob.glob('*'+system.system+'*.cbor')):
    print(fname)
    base = fname[:-5]

    energy_boundaries, mean_e, my_lnw, my_system, p_exc = compute.read_file(base)
    
    # Create a function for the entropy
    l_function, eee, sss = compute.linear_entropy(energy_boundaries, mean_e, my_lnw)
    plt.plot(E, l_function(E) - l_function(equal_E), label=base)

plt.legend()
plt.savefig(system.system+'.svg')
plt.show()