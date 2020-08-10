#!/bin/bash

set -ev

cargo build --release --bin wca-grand-transposed --bin lj-grand-transposed --bin erfinv-grand-transposed

echo default to run with two cores, a high upper bound on output file size, and as a restartable job
RQ=(rq run -c 2 --max-output 20 -R)

GENERIC=(--f 0.5 --translation-scale 0.05 --movie-time "10^(1/4)" --save-time 1)

# ${RQ[@]} -J wca-1 -- ../target/release/wca-grand-transposed --save-as wca-1.cbor --swap-time 1 --max-N 256 --cell-volume 100 --min-T 0.5 ${GENERIC[@]}

${RQ[@]} -J lj-1 -- ../target/release/lj-grand-transposed --save-as lj-1.cbor --swap-time 1 --radius 4 --max-N 38 --min-T 0.01 ${GENERIC[@]}

# ${RQ[@]} -J erfinv-1 -- ../target/release/erfinv-grand-transposed --save-as erfinv-1.cbor --swap-time 1 --mean-energy 0 --min-T 0.01 --max-N 10  ${GENERIC[@]}

sleep 5

tail -f *.out
