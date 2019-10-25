set -ev

# First we need to build
cargo build --release

# Use rq to run the file

rq run -R -J 'lj-sad-31-minT0005-de005-s1' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.05 --save-as lj-sad-31-minT001-de005-s1.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3 --seed=1

rq run -R -J 'lj-sad-31-minT0005-de005-s2' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.05 --save-as lj-sad-31-minT001-de005-s2.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3 --seed=2

rq run -R -J 'lj-sad-31-minT0005-de005-s3' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.05 --save-as lj-sad-31-minT001-de005-s3.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3 --seed=3

rq run -R -J 'lj-sad-31-minT0005-de005-s4' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.05 --save-as lj-sad-31-minT001-de005-s4.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3 --seed=4

rq run -R -J 'lj-sad-31-minT0005-de001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-31-minT0005-de001.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

rq run -R -J 'lj-sad-31-minT0005-de002' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.02 --save-as lj-sad-31-minT0005-de002.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

rq run -R -J 'lj-sad-31-minT0005-de003' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.03 --save-as lj-sad-31-minT0005-de003.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

rq run -R -J 'lj-sad-31-minT0005-de004' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.04 --save-as lj-sad-31-minT0005-de004.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

rq run -R -J 'lj-samc-1e6-de001' target/release/lj-cluster --N 31 --max-allowed-energy 0 --samc-t0 1e6 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-samc-1e6-31-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

rq run -R -J 'lj-samc-1e6-de005' target/release/lj-cluster --N 31 --max-allowed-energy 0 --samc-t0 1e6 --translation-scale 0.05 --energy-bin 0.05 --save-as lj-samc-1e6-31-bin005.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

rq run -R -J 'lj-wl-31-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.4 --wl --translation-scale 0.5 --energy-bin 0.01 --save-as lj-wl-31-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

rq run -R -J 'lj-inv-t-wl-31-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.4 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-31-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

## ------ RUN ON A CONSTRAINED ENERGY RANGE ------ ## # The groundstate is -133.586422
rq run -R -J 'lj-sad-benchmark' target/release/lj-cluster --N 31 --max-allowed-energy=-110 --min-allowed-energy=-133.53 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.05 --save-as lj-sad-31-minT0005-de005-s3-max110.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3 --seed=3

rq run -R -J 'lj-samc-1e4-benchmark' target/release/lj-cluster --N 31 --max-allowed-energy=-110 --min-allowed-energy=-133.53 --samc-t0 1e4 --translation-scale 0.03 --energy-bin 0.01 --save-as lj-samc-1e4-31-benchmark-max110.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3

rq run -R -J 'lj-samc-1e5-benchmark' target/release/lj-cluster --N 31 --max-allowed-energy=-110 --min-allowed-energy=-133.53 --samc-t0 1e5 --translation-scale 0.03 --energy-bin 0.01 --save-as lj-samc-1e5-31-benchmark-max110.yaml --movie-time '10^(1/8)' --save-time 0.5 --n-radial 151 --radius 3
