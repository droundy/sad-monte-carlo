#!/bin/bash

set -ev

cargo build --release --bin transposed --bin binning

ERFINV=(../target/release/transposed --translation-scale 0.05 --movie-time "10^(1/4)" --fake-erfinv-mean-energy 0 --fake-erfinv-N 10 --max-iter 1e10)

RUN=(../target/release/transposed --translation-scale 0.05 --movie-time "10^(1/4)" --max-iter 1e10)

SAD=(../target/release/binning --translation-scale 0.05 --movie-time "10^(1/4)" --max-iter 1e10)

rq run -R -J t-erfinv -- ${ERFINV[@]//}  --save-as t-erfinv.yaml --min-T 0.1 --f '0.5'

for de in 0.01 0.1 1.0; do
  rq run -R -J sad-gaussian-$de -- ${SAD[@]//} --histogram-bin $de --fake-gaussian-sigma 0.1 --save-as sad-gaussian-$de.yaml --sad-min-T 0.01

  rq run -R -J sad-linear-$de --  ${SAD[@]//} --histogram-bin $de --fake-linear --save-as sad-linear-$de.yaml --sad-min-T 0.01

  rq run -R -J sad-quadratic-$de -- ${SAD[@]//} --histogram-bin $de --fake-quadratic-dimensions 3 --save-as sad-quadratic-$de.yaml --sad-min-T 0.01
done

rq run -R -J t-gaussian -- ${RUN[@]//}  --fake-gaussian-sigma 0.1 --save-as t-gaussian.yaml --min-T 0.01 --f '0.5'

rq run -R -J t-linear --  ${RUN[@]//} --fake-linear --save-as t-linear.yaml --min-T 0.01 --f '0.5'

rq run -R -J t-quadratic -- ${RUN[@]//}  --fake-quadratic-dimensions 3 --save-as t-quadratic.yaml --min-T 0.01 --f '0.5'

sleep 5

tail -f *.out
