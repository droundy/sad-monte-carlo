#!/bin/sh

set -ev

cargo build --release

rq run -R -J sad_1000 target/release/energy-mc --N 1000 --filling-fraction 0.22 --sad-min-T 0.1 --well-width 1.3 --movie-time '10^(1/4)' --save-as sad_1000.yaml --translation-scale 0.05 --save-time 0.5

rq run -R -J samc-1e9_1000 target/release/energy-mc --N 1000 --filling-fraction 0.22 --samc-t0 1e9 --well-width 1.3 --movie-time '10^(1/4)' --save-as samc-1e9_1000.yaml --translation-scale 0.05 --save-time 0.5

rq run -R -J samc-1e8_1000 target/release/energy-mc --N 1000 --filling-fraction 0.22 --samc-t0 1e8 --well-width 1.3 --movie-time '10^(1/4)' --save-as samc-1e8_1000.yaml --translation-scale 0.05 --save-time 0.5

rq run -R -J samc-1e7_1000 target/release/energy-mc --N 1000 --filling-fraction 0.22 --samc-t0 1e8 --well-width 1.3 --movie-time '10^(1/4)' --save-as samc-1e7_1000.yaml --translation-scale 0.05 --save-time 0.5

rq run -R -J wl_1000 target/release/energy-mc --N 1000 --filling-fraction 0.22 --wl --well-width 1.3 --movie-time '10^(1/4)' --save-as wl_1000.yaml --translation-scale 0.05 --save-time 0.5
