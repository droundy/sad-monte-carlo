#!/usr/bin/python3

import numpy as np
import yaml
import cbor
import argparse
import sys
import os
import glob
import scipy.constants as const
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import colorcet as cc

import compute

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('base', nargs='*', help='the yaml or cbor files')

args = parser.parse_args()


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


# Read Data
bases = []

# each file has different path (including extension) so concatenating is easy
for base in args.base:
    # change base to have the cbor files. currently has the directory
    if '.cbor' in base or '.yaml' in base:
        base = base[:-5]
    bases.append(base)

for base in bases:
    print('reading', base)

    energy_boundaries = np.loadtxt(base+'-energy-boundaries.dat')
    de = abs(np.diff(energy_boundaries))
    lnw = np.loadtxt(base+'-lnw.dat')
    mean_e = np.loadtxt(base+'-mean-energy.dat')

    if energy_boundaries[0] < energy_boundaries[-1]:
        energy_boundaries = np.flip(energy_boundaries)
        mean_e = np.flip(mean_e)
        lnw = np.flip(lnw)
    reference_function, _, _ = compute.linear_entropy(
        energy_boundaries, mean_e, lnw)

    energy_boundaries, mean_e, _, system, p_exc_reference = compute.read_file(
        base)

    p_reference, T_reference = compute.pressure_temperature(
        system['density'], energy_boundaries, mean_e, p_exc_reference)

    small_dT = 1/16
    big_dT = 5
    T = np.array(list(np.arange(0.125, 3 + small_dT/2, small_dT)) + list(np.arange(5, 30+big_dT/2, big_dT)))
    p = canonical(T, mean_e, p_reference, reference_function)
    np.savetxt(base+'-p-vs-T.dat',
               np.transpose(np.vstack([T, p])),
               fmt='%.4g')
