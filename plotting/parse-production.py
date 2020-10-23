#!/usr/bin/python3

import numpy as np
import yaml, cbor, argparse, sys, os, readsystem

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
    print('done reading', fname)

    base = os.path.splitext(fname)[0]

    histogram = np.array(data_loaded['histogram'])
    total_energy = np.array([float(f) for f in data_loaded['total_energy']])
    lnw = np.array(data_loaded['lnw'])
    energy_boundaries = np.array(data_loaded['energy_boundaries'])

    print('histogram', type(histogram[0]))
    print('total_energy', type(total_energy[0]))
    print('total_energy', total_energy)
    print('histogram', histogram)
    mean_energy = total_energy/histogram

    lnw += np.log(histogram+1)
    lnw -= lnw.max()

    np.savetxt(base+'-energy-boundaries.dat', energy_boundaries)
    np.savetxt(base+'-mean-energy.dat', mean_energy)
    np.savetxt(base+'-lnw.dat', lnw) #includes unbounded extremes

    system = readsystem.readsystem(data_loaded)

    with open(base+'-system.dat', 'w') as f:
        yaml.safe_dump(system, f)
