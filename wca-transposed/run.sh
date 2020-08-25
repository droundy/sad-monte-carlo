#!/bin/bash

set -ev

cargo build --release --bin wca-transposed --bin wca-binning-energy

echo default to run with two cores, a high upper bound on output file size, and as a restartable job
RQ=(rq run --max-output 20 -R)

GENERIC=(--translation-scale 0.05 --movie-time "10^(1/4)" --save-time 1)

${RQ[@]} -J t-wca -- ../target/release/wca-transposed --save-as t-wca.yaml --N 32 --reduced-density 1.2 --min-T 0.1 --f 0.5  ${GENERIC[@]}

${RQ[@]} -J sad-wca -- ../target/release/wca-binning-energy --save-as sad-wca.cbor --N 32 --reduced-density 1.2 --sad-min-T 0.1 --max-allowed-energy 640 --histogram-bin 0.1 ${GENERIC[@]}


# sleep 5

# tail -f *.out
