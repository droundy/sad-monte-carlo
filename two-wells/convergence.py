#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import system, compute
import glob

E = np.linspace(-system.h_small, 0, 10000)
E = 0.5*(E[1:] + E[:-1])
plt.plot(E, system.S(E), label='exact')

for fname in glob.glob('z-*'+system.system+'*.cbor'):
    base = fname[:-5]

    energy_boundaries, mean_e, my_lnw, system, p_exc = compute.read_file(f)
    
    # Create a function for the entropy
    l_function, eee, sss = compute.linear_entropy(energy_boundaries, mean_e, my_lnw)

plt.legend()
plt.show()