#!/bin/bash

set -ev

cargo build --release --bin transposed # --bin binning
rm -f target/release/binning

echo default to run with two cores, a high upper bound on output file size, and as a restartable job
RQ=(rq run --max-output 20 -R)

GENERIC=(--wca-N 32 --wca-reduced-density 1.2 --translation-scale 0.05 --movie-time "10^(1/4)" --save-time 1)

${RQ[@]} -J t-wca -- ../target/release/transposed --save-as t-wca.yaml --f 0.5 --min-T 0.1 ${GENERIC[@]}

# ${RQ[@]} -J sad-wca -- ../target/release/binning --save-as sad-wca.cbor --max-allowed-energy 640 --histogram-bin 0.1 --sad-min-T 0.1 ${GENERIC[@]}


sleep 5

tail -f *.out
