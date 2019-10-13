#!/usr/bin/python3

import sys, os

N = 38
maxrad = 3.5

dE = 0.06
minT = 0.01

n_radial = maxrad*50+1

froot = f"lj{N}-de{dE}-minT{minT}-sad"

# Need to also teach this to use different algorithms

cmd = f"rq cancel {froot} && cargo build --release && rq run -J {froot} target/release/lj-cluster --N {N} --max-allowed-energy 0 --sad-min-T {minT} --translation-scale 0.1 --energy-bin {dE} --save-as {froot}.yaml --movie-time '10^(1/4)' --save-time 1 --n-radial {n_radial} --radius {maxrad}"
print(cmd)
os.system(cmd)
