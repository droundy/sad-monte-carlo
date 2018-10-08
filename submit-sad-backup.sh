#!/bin/sh

set -ev

cargo build --release

rq run -R -J samc-T13 target/release/energy-mc --N 50 --filling-fraction 0.3 --well-width 1.3 --movie-time '10^(1/8)' --save-as samc-T13.yaml --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-248 --max-allowed-energy=-114 --samc-t0 1e4

rq run -R -J wl-T13 target/release/energy-mc --N 50 --filling-fraction 0.3 --well-width 1.3 --movie-time '10^(1/8)' --save-as wl-T13.yaml --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-248 --max-allowed-energy=-114 --wl

rq run -R -J sad-T13 target/release/energy-mc --N 50 --filling-fraction 0.3 --sad-min-T '1/3' --well-width 1.3 --movie-time '10^(1/8)' --save-as sad-T13.yaml --translation-scale 0.05 --save-time 0.5

rq run -R -J sad-slow-T13 target/release/energy-mc --N 50 --filling-fraction 0.3 --sad-min-T '1/3' --well-width 1.3 --movie-time '10^(1/8)' --save-as sad-slow-T13.yaml --translation-scale 0.005 --save-time 0.5

sleep 3
tail -f *T13.out
