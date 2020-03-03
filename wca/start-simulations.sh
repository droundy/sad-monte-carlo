#!/bin/sh

set -ev

rq run -J wca-1-1-0.5 -- cargo run --release --bin energy-wca -- --reduced-density 1 --N 32 --sad-min-T 0.5 --energy-bin 1.0 --max-allowed-energy 640 --acceptance-rate 0.5 --movie-time '10^(1/4)' --save-time 0.5 --save-as sad-n1.0-de1.0-minT-0.5.yaml
