#!/bin/sh

set -ev

cargo build --release

rq run -R -J samc-1e5-256 target/release/energy-mc --N 256 --filling-fraction 0.17 --well-width 1.5 --movie-time '10^(1/8)' --save-as samc-1e5-256.yaml --translation-scale 0.05 --save-time 0.5 --samc-t0 1e5 --max-allowed-energy=-440 --min-allowed-energy=-915 
rq run -R -J samc-1e6-256 target/release/energy-mc --N 256 --filling-fraction 0.17 --well-width 1.5 --movie-time '10^(1/8)' --save-as samc-1e6-256.yaml --translation-scale 0.05 --save-time 0.5 --samc-t0 1e6 --max-allowed-energy=-440 --min-allowed-energy=-915

rq run -R -J sad-256 target/release/energy-mc --N 256 --filling-fraction 0.17 --sad-min-T 1 --well-width 1.5 --movie-time '10^(1/8)' --save-as sad-256.yaml --translation-scale 0.05 --save-time 0.5

sleep 3
tail -f *-256.out
