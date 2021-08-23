#!/usr/bin/python3

import numpy as np
import yaml
import cbor
import argparse
import sys
import os
import readsystem

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('fname', nargs='*', help='the yaml or cbor file')
parser.add_argument("--reparse", help="parse file even if it has already been parsed",
                    action="store_true")

args = parser.parse_args()

for fname in args.fname:
    base = os.path.splitext(fname)[0]
    if not args.reparse and os.path.exists(base+'-lnw.dat') or '*' in base:
        continue

    with open(fname, 'rb') as stream:
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

    replicas = data_loaded['replicas']

    above_total = np.array([float(r['above_total']) for r in replicas])
    above_total_squared = np.array([float(r['above_total_squared']) for r in replicas])
    above_count = np.array([r['above_count'] for r in replicas])
    below_count = np.array([r['below_count'] for r in replicas])
    count = above_count+below_count

    mean_energy = []
    mean_energy_squared = []
    lnw = []
    remaining_entopy = 0
    above_fraction = []
    for i in range(0, len(above_count)):
        lnw.append(remaining_entopy + np.log(above_count[i]/count[i]))
        remaining_entopy = remaining_entopy + np.log(below_count[i]/count[i])
        mean_energy.append(above_total[i]/above_count[i])
        mean_energy_squared.append(above_total_squared[i]/above_count[i])
        above_fraction.append([replicas[i]['cutoff_energy'], above_count[i]/count[i],
                               replicas[i]['max_energy'] - replicas[i]['cutoff_energy'],
                               above_total[i]/above_count[i], replicas[i]['unique_visitors']])

    lnw.append(remaining_entopy)
    if replicas[-1]['below_count'] > 0:
        mean_energy.append(
            replicas[-1]['below_total']/replicas[-1]['below_count'])
        mean_energy_squared.append(
            replicas[-1]['below_total_squared']/replicas[-1]['below_count'])
    else:
        mean_energy.append(np.NaN)
        mean_energy_squared.append(np.NaN)

    lnw = np.array(lnw)
    lnw -= np.log(np.sum(np.exp(lnw)))  # w = w / sum(w) to normalize
    mean_energy = np.array(mean_energy)

    energy_boundaries = []
    for i in range(0, len(above_count)):
        energy_boundaries.append(replicas[i]['cutoff_energy'])
    energy_boundaries = np.array(energy_boundaries)

    system = readsystem.readsystem(replicas[0])
    for e in replicas[0]['above_extra'].keys():
        extra = []
        for r in replicas:
            if e in r['above_extra']:
                extra.append(r['above_extra'][e][0]/r['above_extra'][e][1])
            else:
                extra.append(np.NaN)
        extra.append(np.NaN)
        np.savetxt(base+f'-{e}.dat', extra)
            

    np.savetxt(base+'-above-fraction.dat', above_fraction, fmt='%13.7g')
    np.savetxt(base+'-energy-boundaries.dat', energy_boundaries)
    np.savetxt(base+'-mean-energy.dat', mean_energy)
    np.savetxt(base+'-mean-energy-squared.dat', mean_energy_squared)
    np.savetxt(base+'-lnw.dat', lnw)  # includes unbounded extremes

    with open(base+'-system.dat', 'w') as f:
        yaml.safe_dump(system, f)
