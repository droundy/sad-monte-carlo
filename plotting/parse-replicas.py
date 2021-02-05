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

    replicas = data_loaded['replicas']

    above_total = np.array([float(r['above_total']) for r in replicas])
    above_count = np.array([r['above_count'] for r in replicas])
    below_count = np.array([r['below_count'] for r in replicas])
    count = above_count+below_count

    mean_energy = []
    lnw = []
    remaining_entopy = 0
    for i in range(0,len(above_count)):
        lnw.append(remaining_entopy + np.log(above_count[i]/count[i]))
        remaining_entopy = remaining_entopy + np.log(below_count[i]/count[i])
        mean_energy.append(above_total[i]/above_count[i])

    lnw.append(remaining_entopy)
    if replicas[-1]['below_count'] > 0:
        mean_energy.append(replicas[-1]['below_total']/replicas[-1]['below_count'])
    else:
        mean_energy.append(np.NaN)

    lnw = np.array(lnw)
    lnw -= np.log(np.sum(np.exp(lnw))) # w = w / sum(w) to normalize
    print('norm', np.sum(np.exp(lnw)))
    mean_energy = np.array(mean_energy)


    energy_boundaries = []
    for i in range(0,len(above_count)):
        energy_boundaries.append(replicas[i]['cutoff_energy'])
    energy_boundaries = np.array(energy_boundaries)

    system = readsystem.readsystem(replicas[0])

    np.savetxt(base+'-energy-boundaries.dat', energy_boundaries)
    np.savetxt(base+'-mean-energy.dat', mean_energy)
    np.savetxt(base+'-lnw.dat', lnw) #includes unbounded extremes

    with open(base+'-system.dat', 'w') as f:
        yaml.safe_dump(system, f)
