#!/bin/bash

set -ev

cargo build --release --bin wca-grand-transposed --bin lj-grand-transposed --bin erfinv-grand-transposed

GENERIC=(--f 0.5 --translation-scale 0.05 --movie-time "10^(1/4)" --save-time 1)

rq run -R -o wca-1 -- ../target/release/wca-grand-transposed --save-as wca-1.cbor --swap-time 1 --max-N 256 --cell-volume 100 --min-T 0.5 ${GENERIC[@]}

rq run -R -o lj-1 -- ../target/release/lj-grand-transposed --save-as lj-1.cbor --swap-time 1 --radius 4 --max-N 38 --min-T 0.01 ${GENERIC[@]}

rq run -R -o erfinv-1 -- ../target/release/erfinv-grand-transposed --save-as erfinv-1.cbor --swap-time 1 --mean-energy 0 --min-T 0.01 --max-N 10  ${GENERIC[@]}