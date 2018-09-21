#!/bin/sh

set -ev

cargo build --release

rq run -R -J grand-samc-1e9 target/release/energy-number-mc --samc-t0 1e9 --well-width 1.3 --movie-time '10^(1/4)' --save-as grand-samc-1e9.yaml --save-time 0.5 --cell-volume 1000 --acceptance-rate 0.5
