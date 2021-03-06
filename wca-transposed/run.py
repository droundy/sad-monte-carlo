#!/usr/bin/python3

import os, time

assert(not os.system('cargo build --release --bin binning --bin replicas --bin production'))

# default to run with two cores, a high upper bound on output file size, and as a restartable job
rq = 'rq run --max-output 20 -R'

time_params = "--max-iter 5e10 --movie-time '10^(1/4)' --save-time 1"
min_T = 0.1

def wca(density, N, shift_N):
    if shift_N != 0:
        Ncorrect = N + shift_N
        return f' --wca-N {Ncorrect} --wca-cell-volume "{N}/{density}"'
    else:
        return f' --wca-N {N} --wca-reduced-density {density}'

def name_wca(prefix, density, N, shift_N):
    name = f'{prefix}-wca-{density}-{N}'
    if shift_N > 0:
        name += f'+{shift_N}'
    elif shift_N < 0:
        name += f'{shift_N}'
    return name

def run_replicas(density, N, shift_N=0, job='run'):
    name = name_wca('r', density, N, shift_N)
    wca_params = wca(density, N, shift_N)
    if job == 'production':
        os.system(f'{rq} -J p-{name} ../target/release/production --base {name} --subdivide 16 --save-as p-{name}.cbor {wca_params} {time_params}')
        os.system(f'{rq} -J P-{name} ../target/release/production --base {name} --seed 2 --save-as P-{name}.cbor {wca_params} {time_params}')
    elif job == 'parse':
        os.system(f'{rq} -c all -J parse-{name} ../plotting/parse-replicas.py {name}.cbor')
    else:
        os.system(f'{rq} -c all -J {name} ../target/release/replicas --save-as {name}.cbor {wca_params} --min-T {min_T} {time_params}')

def run_sad(dE, density, N, shift_N=0, job='run'):
    name = name_wca('s', density, N, shift_N)
    wca_params = wca(density, N, shift_N)
    max_energy = N*20
    if job == 'production':
        pass
        # No point running a production run after sad, since it'll get stuck at
        # high energies!
        # os.system(f'{rq} -J p-{name} ../target/release/production --base {name} --subdivide 16 --save-as p-{name}.cbor {wca_params} {time_params}')
    elif job == 'parse':
        pass
        # os.system(f'{rq} -J parse-{name} ../plotting/parse-binning.py {name}.cbor')
    else:
        os.system(f'{rq} -J {name} ../target/release/binning --save-as {name}.cbor {wca_params} --sad-min-T {min_T} --max-allowed-energy {max_energy} --translation-scale 0.005 --histogram-bin {dE} {time_params}')

jobs = ['run', 'parse', 'production']
for density in [1.2]:
    for N in [32, 108, 256]:
        for shift_N in [0]: # [0,-1,1]:
            for job in jobs:
                run_replicas(density, N, shift_N, job=job)
                run_sad(min_T, density, N, shift_N, job=job)
