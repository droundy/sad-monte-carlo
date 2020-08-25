#!/usr/bin/python3

import numpy as np
import yaml, cbor, argparse, sys, os

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('fname', nargs='*', help = 'the yaml or cbor file')

args = parser.parse_args()

for fname in args.fname:
    with open(fname,'rb') as stream:
        if 'yaml' in fname:
            try:
                data_loaded = yaml.full_load(stream)
            except IOError:
                print('An error occurred trying to read the file.')
        elif 'cbor' in fname:
            try:
                data_loaded = cbor.load(stream)
            except IOError:
                print('An error occurred trying to read the file.')
        else:
            print('What kind of file is this?!')
            exit(1)

    base = os.path.splitext(fname)[0]

    total_energy = np.array(data_loaded['total_energy'], dtype=float)
    histogram = np.array(data_loaded['histogram'])
    rel_bins = np.array(data_loaded['rel_bins'])
    bin_norm = data_loaded['bin_norm']
    max_energy = data_loaded['max_energy']
    min_energy = data_loaded['min_energy']
    lnw = np.array(data_loaded['lnw'])
    lnw -= lnw.max()
    lnw -= np.log(np.sum(np.exp(lnw))) # w = w / sum(w) to normalize

    energy_boundaries = [max_energy]
    energy_per_rel_bin = 1/bin_norm*(max_energy-min_energy)
    for b in rel_bins:
        energy_boundaries.append( energy_boundaries[-1] - b*energy_per_rel_bin)
    energy_boundaries = np.array(energy_boundaries)

    np.savetxt(base+'-energy-boundaries.dat', energy_boundaries)

    mean_energy = total_energy/histogram #includes unbounded extremes
    np.savetxt(base+'-mean-energy.dat', mean_energy)

    np.savetxt(base+'-lnw.dat', lnw) #includes unbounded extremes
