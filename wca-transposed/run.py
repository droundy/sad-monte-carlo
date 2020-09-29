#!/usr/bin/python3

import os, time

assert(not os.system('cargo build --release --bin binning --bin replicas'))

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

def run_replicas(density, N, shift_N=0):
    name = name_wca('r', density, N, shift_N)
    wca_params = wca(density, N, shift_N)
    os.system(f'{rq} -c all -J {name} ../target/release/replicas --save-as {name}.cbor {wca_params} --min-T {min_T} {time_params}')

def run_sad(dE, density, N, shift_N=0):
    name = name_wca('s', density, N, shift_N)
    wca_params = wca(density, N, shift_N)
    max_energy = N*20
    os.system(f'{rq} -J {name} ../target/release/binning --save-as {name}.cbor {wca_params} --sad-min-T {min_T} --max-allowed-energy {max_energy} --translation-scale 0.005 --histogram-bin {dE} {time_params}')

for density in [1.2]:
    for N in [32, 108, 256]:
        for shift_N in [-1,0,1]:
            run_replicas(density, N, shift_N)
            time.sleep(1)
            run_sad(min_T, density, N, shift_N)
            time.sleep(1)
