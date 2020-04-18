#!/bin/sh

set -ev

cargo build --release

rq run -J wca-1-1-0.5 -- ../target/release/energy-wca --reduced-density 1 --N 32 --sad-min-T 0.5 --energy-bin 1.0 --max-allowed-energy 640 --acceptance-rate 0.5 --movie-time '10^(1/4)' --save-time 0.5 --save-as sad-n1.0-de1.0-minT-0.5.yaml

rq run -J wca-1-0.1-0.5 -- ../target/release/energy-wca --reduced-density 1 --N 32 --sad-min-T 0.5 --energy-bin 0.1 --max-allowed-energy 640 --acceptance-rate 0.5 --movie-time '10^(1/4)' --save-time 0.5 --save-as sad-n1.0-de0.1-minT-0.5.yaml

rq run -J wca-1-0.01-0.5 -- ../target/release/energy-wca --reduced-density 1 --N 32 --sad-min-T 0.5 --energy-bin 0.01 --max-allowed-energy 640 --acceptance-rate 0.5 --movie-time '10^(1/4)' --save-time 0.5 --save-as sad-n1.0-de0.01-minT-0.5.cbor

rq run -J new-wca-n1.0-de1-minT-0.5 -- ../target/release/wca-binning-energy --reduced-density 1 --N 32 --sad-min-T 0.5 --histogram-bin 1 --max-allowed-energy 640 --acceptance-rate 0.5 --movie-time '10^(1/4)' --save-time 0.5 --movie-dir new-sad-n1.0-de1-minT-0.5 --save-as new-sad-n1.0-de1-minT-0.5.yaml

rq run -J new-wl-n1.0-de1 -- ../target/release/wca-binning-energy --reduced-density 1 --N 32 --wl --wl-min-gamma 1e-4 --histogram-bin 1 --min-allowed-energy 18.5 --max-allowed-energy 640 --acceptance-rate 0.5 --movie-time '10^(1/4)' --save-time 0.5 --movie-dir new-wl-n1.0-de1 --save-as new-wl-n1.0-de1.yaml
