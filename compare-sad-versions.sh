#!/bin/sh

set -ev

cargo build --release

rq run -R -J sad-t-T13 target/release/energy-mc --N 50 --filling-fraction 0.3 --sad-min-T '1/3' --well-width 1.3 --movie-time '10^(1/8)' --save-as sad-t-T13.yaml --translation-scale 0.05 --save-time 0.5 --sad-version-sad-over-t

rq run -R -J sad-t2-T13 target/release/energy-mc --N 50 --filling-fraction 0.3 --sad-min-T '1/3' --well-width 1.3 --movie-time '10^(1/8)' --save-as sad-t2-T13.yaml --translation-scale 0.05 --save-time 0.5 --sad-version-sad-over-t2

rq run -R -J sad-original-T13 target/release/energy-mc --N 50 --filling-fraction 0.3 --sad-min-T '1/3' --well-width 1.3 --movie-time '10^(1/8)' --save-as sad-original-T13.yaml --translation-scale 0.05 --save-time 0.5 --sad-version-sad

# rq run -R -J sad-slow-T13 target/release/energy-mc --N 50 --filling-fraction 0.3 --sad-min-T '1/3' --well-width 1.3 --movie-time '10^(1/8)' --save-as sad-slow-original-T13.yaml --translation-scale 0.005 --save-time 0.5 --sad-version-sad

sleep 3
tail -f sad*T13.out
