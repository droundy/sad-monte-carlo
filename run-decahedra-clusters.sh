set -ev

# First we need to build
cargo build --release

# Use rq to run the file

# The Benchmark Script started by David
# rq run -J 'lj-31-benchmark' target/release/lj-cluster --N 31 --max-allowed-energy=-110 --min-allowed-energy=-133.53 --translation-scale 0.05 --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5 --energy-bin 0.001 --Inv-t-WL --save-as lj-31-benchmark.yaml

# The WL run scripts
# rq run -R --max-output 10 -J 'lj-wl-31-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --wl --translation-scale 0.05 --energy-bin 0.01 --save-as lj-wl-31-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

# The 1/t-WL run scripts
rq run -R --max-output 10 -J 'lj-inv-t-wl-75-bin001' target/release/lj-cluster --N 75 --max-allowed-energy=0 --min-allowed-energy=-397.4 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-75-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 10 -J 'lj-inv-t-wl-98-bin001' target/release/lj-cluster --N 98 --max-allowed-energy=0 --min-allowed-energy=-543.6 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-98-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 10 -J 'lj-inv-t-wl-104-bin001' target/release/lj-cluster --N 104 --max-allowed-energy=0 --min-allowed-energy=-582.0 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-104-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5


# The SAMC run scripts t0 = 1e5
# rq run -R --max-output 10 -J 'lj-samc-31-1e5-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --samc-t0 1e5 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-samc-31-1e5-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

# The SAD run scripts
rq run -R --max-output 10 -J 'lj-sad-75-bin001' target/release/lj-cluster --N 75 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-75-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 10 -J 'lj-sad-98-bin001' target/release/lj-cluster --N 98 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-98-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 10 -J 'lj-sad-104-bin001' target/release/lj-cluster --N 104 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-104-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5
