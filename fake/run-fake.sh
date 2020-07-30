#!/bin/bash

set -ev

cargo build --release --bin fake-transposed --bin fake-binning --bin erfinv-transposed

ERFINV=(../target/release/erfinv-transposed --translation-scale 0.05 --movie-time "10^(1/4)" --mean-energy 0 --N 10)

RUN=(../target/release/fake-transposed --translation-scale 0.05 --movie-time "10^(1/4)")

SAD=(../target/release/fake-binning --translation-scale 0.05 --movie-time "10^(1/4)" --histogram-bin 0.01)

rq run -R -J t-erfinv -- ${ERFINV[@]//}  --save-as t-erfinv.yaml --min-T 0.1 --f '0.5'

rq run -R -J sad-gaussian-0.01 -- ${SAD[@]//} --gaussian-sigma 0.1 --save-as sad-gaussian-0.01.yaml --sad-min-T 0.01

rq run -R -J sad-linear-0.01 --  ${SAD[@]//} --linear --save-as sad-linear-0.01.yaml --sad-min-T 0.01

rq run -R -J sad-quadratic-0.01 -- ${SAD[@]//}  --quadratic-dimensions 3 --save-as sad-quadratic-0.01.yaml --sad-min-T 0.01

rq run -R -J t-gaussian -- ${RUN[@]//}  --gaussian-sigma 0.1 --save-as t-gaussian.yaml --min-T 0.01 --f '0.5'

rq run -R -J t-linear --  ${RUN[@]//} --linear --save-as t-linear.yaml --min-T 0.01 --f '0.5'

rq run -R -J t-quadratic -- ${RUN[@]//}  --quadratic-dimensions 3 --save-as t-quadratic.yaml --min-T 0.01 --f '0.5'

sleep 5

tail -f *.out
