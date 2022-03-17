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
i=0
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

    try:
        i+=1
        mean_Es = []
        var_Es = []
        cvs = []
        Ts = []
        for r in data_loaded['replicas']:
            
            #print((r, i))
            
            r_clean = dict()
            r_clean['accepted_swap'] = r['accepted_swap_count'] / (
                            r['rejected_swap_count'] + r['accepted_swap_count']
                        )
            r_clean['accepted_moves'] = r['accepted_count'] / (
                        r['rejected_count'] + r['accepted_count']
                        )

            num_moves = r['accepted_count'] + r['rejected_count'] + \
                r['accepted_swap_count'] + r['rejected_swap_count'] - r['ignored_count'] #number of samples

            swaps = 0
            steps = 0
            for k in r.keys():
                if 'swap_count' in k:
                    swaps += r[k]
                elif 'count' in k:
                    steps += r[k]
            ignored = r['ignored_count']
            T=r['T']
            #print(f'T: {T}\nswaps: {swaps}\nsteps: {steps}\nignored: {ignored}')

            r_clean['mean_E_squared'] = r['total_energy_squared']/num_moves

            r_clean['mean_E'] = r['total_energy']/num_moves

            r_clean['T'] = r['T']

            r_clean['Cv'] = (r_clean['mean_E_squared'] - r_clean['mean_E']**2) /(r_clean['T']**2)


            mean_Es.append(r_clean['mean_E'])
            var_Es.append(r_clean['mean_E_squared'] - (r_clean['mean_E'])**2)
            cvs.append(r_clean['Cv'])
            Ts.append(r_clean['T'])


        np.savez(f'{fname[:-5]}.npz', T = Ts, Cv = cvs, mean_E=mean_Es, var_E=var_Es)

    except ZeroDivisionError as err:
        print(f'Skipping file {fname} due to {err}')
        pass
    except BaseException as err:
        raise(f"Unexpected {err}, {type(err)}")


