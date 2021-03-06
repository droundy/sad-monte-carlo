#!/bin/bash

set -ev

cargo build --release --bin grand-transposed

echo default to run with two cores, a high upper bound on output file size, and as a restartable job
RQ=(rq run -c 2 --max-output 20 -R)

GENERIC=(--f 0.5 --translation-scale 0.05 --movie-time "10^(1/4)" --save-time 1)

${RQ[@]} -J wca-1 -- ../target/release/grand-transposed --save-as wca-1.cbor --swap-time 1 --max-N 256 --wca-cell-volume 100 --min-T 0.5 ${GENERIC[@]}

${RQ[@]} -J lj-1 -- ../target/release/grand-transposed --save-as lj-1.cbor --swap-time 1 --lj-radius 4 --max-N 38 --min-T 0.01 ${GENERIC[@]}

${RQ[@]} -J erfinv-1 -- ../target/release/grand-transposed --save-as erfinv-1.cbor --swap-time 1 --fake-erfinv-mean-energy 0 --min-T 0.01 --max-N 10  ${GENERIC[@]}

sleep 5

tail -f *.out
