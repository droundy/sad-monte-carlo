#!/bin/bash

set -ev

cargo build --release --bin replicas --bin transposed --bin binning

ERFINV=(../target/release/transposed --translation-scale 0.05 --movie-time "10^(1/4)" --fake-erfinv-mean-energy 0 --fake-erfinv-N 10 --max-iter 1e10)

TRA=(../target/release/transposed --translation-scale 0.05 --movie-time "10^(1/4)" --max-iter 1e10)

REP=(../target/release/replicas --movie-time "10^(1/4)" --max-iter 1e10)

SAD=(../target/release/binning --translation-scale 0.05 --movie-time "10^(1/4)" --max-iter 1e10)

rq run -c all --max-output=30 -R -J r-erfinv -- ${REP[@]//} --save-as r-erfinv.yaml --min-T 0.01 --fake-erfinv-mean-energy 0 --fake-erfinv-N 10

rq run -c all --max-output=30 -R -J r-gaussian -- ${REP[@]//} --save-as r-gaussian.yaml --min-T 0.001 --fake-gaussian-sigma 0.1

rq run -c all --max-output=30 -R -J r-linear -- ${REP[@]//} --save-as r-linear.yaml --min-T 0.01 --fake-linear

rq run -c all --max-output=30 -R -J r-quadratic -- ${REP[@]//} --save-as r-quadratic.yaml --min-T 0.001 --fake-quadratic-dimensions 3

rq run -c all --max-output=30 -R -J r-pieces -- ${REP[@]//} --save-as r-pieces.yaml --min-T 0.001 --a 2.0 --b 4.0 --e1 1.0 --e2 3.0

rq run -R -J t-erfinv -- ${ERFINV[@]//}  --save-as t-erfinv.yaml --min-T 0.01 --f '0.5'

for de in 1 0.01 0.1; do
  rq run -R -J sad-erfinv-$de -- ${SAD[@]//} --histogram-bin $de --fake-erfinv-mean-energy 0 --fake-erfinv-N 10 --save-as sad-erfinv-$de.yaml --sad-min-T 0.01
done

for de in 0.001 0.01 0.1; do
  rq run -R -J sad-gaussian-$de -- ${SAD[@]//} --histogram-bin $de --fake-gaussian-sigma 0.1 --save-as sad-gaussian-$de.yaml --sad-min-T 0.001

  rq run -R -J sad-linear-$de --  ${SAD[@]//} --histogram-bin $de --fake-linear --save-as sad-linear-$de.yaml --sad-min-T 0.01

  rq run -R -J sad-quadratic-$de -- ${SAD[@]//} --histogram-bin $de --fake-quadratic-dimensions 3 --save-as sad-quadratic-$de.yaml --sad-min-T 0.001
  
  rq run -R -J sad-pieces-$de -- ${SAD[@]//} --histogram-bin $de --a 2.0 --b 4.0 --e1 1.0 --e2 3.0 --save-as sad-pieces-$de.yaml --sad-min-T 0.001
done

rq run -R -J t-gaussian -- ${TRA[@]//}  --fake-gaussian-sigma 0.1 --save-as t-gaussian.yaml --min-T 0.001 --f '0.5'

rq run -R -J t-linear --  ${TRA[@]//} --fake-linear --save-as t-linear.yaml --min-T 0.01 --f '0.5'

rq run -R -J t-quadratic -- ${TRA[@]//}  --fake-quadratic-dimensions 3 --save-as t-quadratic.yaml --min-T 0.001 --f '0.5'

rq run -R -J t-pieces -- ${TRA[@]//}  --a 2.0 --b 4.0 --e1 2.0 --e2 3.0 --save-as t-pieces.yaml --min-T 0.001 --f '0.5'

sleep 5

tail -f *.out
