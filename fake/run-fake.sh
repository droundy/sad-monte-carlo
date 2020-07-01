#!/bin/bash

set -ev

cargo build --release

RUN=(../target/release/fake-transposed --translation-scale 0.05 --movie-time "10^(1/4)")

rq run -R -J linear --  ${RUN[@]//} --linear --save-as linear.yaml --min-T 0.01 --f '0.5'

rq run -R -J quadratic -- ${RUN[@]//}  --quadratic-dimensions 3 --save-as quadratic.yaml --min-T 0.01 --f '0.5'

sleep 5

tail -f *.out
