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

print(data_loaded)