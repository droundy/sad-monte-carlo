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

correct_S = normalize_S(system.S(E))

correct_S_for_err = correct_S[indices_for_err]

np.savez(system.name(), E=E, correct_S=normalize_S(system.S(E)), correct_S_for_err=correct_S_for_err)

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
    print(fname)
    start_fname = time.process_time()
    base = fname[:-8]
    method = base[:base.find('-')]

    energy_boundaries, mean_e, my_lnw, my_system, p_exc = compute.read_file(base)
    
    # Create a function for the entropy
    l_function, eee, sss = compute.linear_entropy(energy_boundaries, mean_e, my_lnw)
    normalize_S(l_function(E))

    plt.figure('fraction-well')
    mean_which = np.loadtxt(f'{base}-which.dat')
    
    if os.path.exists(f'{base}-histogram.dat'):
        hist = np.loadtxt(f'{base}-histogram.dat')

    errors = []
    moves = []
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
            errors.append(err)
        except:
            pass
    
    np.savez(os.path.join('.',base)+'.npz',
                            E=E,
                            mean_e=mean_e,
                            mean_which=mean_which,
                            hist=hist,
                            S=normalize_S(l_function(E)),
                            moves=moves,
                            errors_S=errors)