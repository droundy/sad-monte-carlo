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


rs = []
for r in data_loaded['replicas']:
    r_clean = dict()
    r_clean['accepted_swap'] = r['accepted_swap_count'] / (
                    r['rejected_swap_count'] + r['accepted_swap_count']
                )
    r_clean['accepted_moves'] = r['accepted_count'] / (
                r['rejected_count'] + r['accepted_count']
                )

    num_moves = r['accepted_count'] + r['rejected_count'] #number of energies

    r_clean['mean_E_squared'] = r['total_energy_squared']/num_moves

    r_clean['mean_E'] = r['total_energy']/num_moves

    r_clean['T'] = r['T']

    r_clean['Cv'] = (r_clean['mean_E_squared'] - r_clean['mean_E']**2) /(r_clean['T']**2)

    print(r_clean)
    rs.append(r_clean)

cvs = []
Ts = []
for r in rs:
    cvs.append(r['Cv'])
    Ts.append(r['T'])

np.savez('Cv.npz', T = Ts, Cv = cvs)


