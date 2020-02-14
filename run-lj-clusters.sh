set -ev

# First we need to build
cargo build --release

# Use rq to run the file

# The Benchmark Script started by David
rq run -J 'lj-31-benchmark' target/release/lj-cluster --N 31 --max-allowed-energy=-110 --min-allowed-energy=-133.53 --translation-scale 0.05 --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5 --energy-bin 0.001 --Inv-t-WL --save-as lj-31-benchmark.yaml

# The WL run scripts
rq run -R --max-output 10 -J 'lj-wl-31-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --wl --translation-scale 0.05 --energy-bin 0.01 --save-as lj-wl-31-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-wl-31-bin0001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --wl --translation-scale 0.05 --energy-bin 0.001 --save-as lj-wl-31-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-wl-31-bin01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --wl --translation-scale 0.05 --energy-bin 0.1 --save-as lj-wl-31-bin01.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

# The 1/t-WL run scripts
rq run -R --max-output 10 -J 'lj-inv-t-wl-31-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-31-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-inv-t-wl-31-bin001-58' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.58 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-31-bin001-58.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5 --seed=3

rq run -R --max-output 10 -J 'lj-inv-t-wl-31-bin001-54' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.54 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-31-bin001-54.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5 --seed=5

rq run -R --max-output 10 -J 'lj-inv-t-wl-31-bin001-52' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.52 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-31-bin001-52.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5 --seed=7

rq run -R --max-output 10 -J 'lj-inv-t-wl-31-bin0001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.001 --save-as lj-inv-t-wl-31-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-inv-t-wl-31-bin01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.1 --save-as lj-inv-t-wl-31-bin01.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

# The SAMC run scripts t0 = 1e5
rq run -R --max-output 10 -J 'lj-samc-31-1e5-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e5 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-samc-31-1e5-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-samc-31-1e5-bin0001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e5 --translation-scale 0.05 --energy-bin 0.001 --save-as lj-samc-31-1e5-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-samc-31-1e5-bin01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e5 --translation-scale 0.05 --energy-bin 0.1 --save-as lj-samc-31-1e5-bin01.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

# The SAMC run scripts t0 = 1e6
rq run -R --max-output 10 -J 'lj-samc-31-1e6-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e6 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-samc-31-1e6-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-samc-31-1e6-bin0001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e6 --translation-scale 0.05 --energy-bin 0.001 --save-as lj-samc-31-1e6-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-samc-31-1e6-bin01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e6 --translation-scale 0.05 --energy-bin 0.1 --save-as lj-samc-31-1e6-bin01.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

# The SAMC run scripts t0 = 1e7
rq run -R --max-output 10 -J 'lj-samc-31-1e7-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e7 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-samc-31-1e7-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-samc-31-1e7-bin0001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e7 --translation-scale 0.05 --energy-bin 0.001 --save-as lj-samc-31-1e7-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

rq run -R --max-output 10 -J 'lj-samc-31-1e7-bin01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e7 --translation-scale 0.05 --energy-bin 0.1 --save-as lj-samc-31-1e7-bin01.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

# The SAD run scripts
rq run -R --max-output 10 -J 'lj-sad-31-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-31-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5 --seed=3

rq run -R --max-output 10 -J 'lj-sad-31-bin0001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.001 --save-as lj-sad-31-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5 --seed=3

rq run -R --max-output 10 -J 'lj-sad-31-bin01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.1 --save-as lj-sad-31-bin01.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5 --seed=3
