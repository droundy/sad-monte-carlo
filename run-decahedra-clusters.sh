set -ev

# First we need to build
cargo build --release

# Use rq to run the file

# MINIMUM PEAK 104 Cs -582.038429 Northby   GROUND STATE 104 C2v -582.086642 Doye2

# The Benchmark Scripts
rq run -R --max-output 10 -J 'lj-75-benchmark' target/release/lj-cluster --N 75 --max-allowed-energy=-150 --min-allowed-energy=-393 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-75-benchmark.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 20 -J 'lj-75-benchmark-tiny' target/release/lj-cluster --N 75 --max-allowed-energy=-150 --min-allowed-energy=-393 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.001 --save-as lj-75-benchmark-tiny.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 20 -J 'lj-38-benchmark' target/release/lj-cluster --N 38 --max-allowed-energy=-100 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-38-benchmark.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3

rq run -R --max-output 200 -J 'lj-38-benchmark-tiny' target/release/lj-cluster --N 38 --max-allowed-energy=-100 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.001 --save-as lj-38-benchmark-tiny.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3

# The WL run scripts
# rq run -R --max-output 10 -J 'lj-wl-31-bin001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --wl --translation-scale 0.05 --energy-bin 0.01 --save-as lj-wl-31-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2.5

# The 1/t-WL run scripts
rq run -R --max-output 10 -J 'lj-inv-t-wl-75-bin001' target/release/lj-cluster --N 75 --max-allowed-energy=0 --min-allowed-energy=-397.4 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-75-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

# The SAD run scripts

# LJ 75, 98, 104 SCRIPTS RUNNING
rq run -R --max-output 10 -J 'lj-sad-75-bin001' target/release/lj-cluster --N 75 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-75-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 500 -J 'lj-sad-75-bin0001' target/release/lj-cluster --N 75 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.001 --save-as lj-sad-75-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 10 -J 'lj-sad-104-bin001' target/release/lj-cluster --N 104 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-104-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3.5

rq run -R --max-output 500 -J 'lj-sad-104-bin0001' target/release/lj-cluster --N 104 --max-allowed-energy=0 --sad-min-T 0.005 --translation-scale 0.05 --energy-bin 0.001 --save-as lj-sad-104-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 5

# LJ 38 BIN SIZE TESTS

# ALL OTHER DATA: SAMC AND WL
# rq run -R --max-output 10 -J 'lj-samc-38-1e5-bin001' target/release/lj-cluster --N 38 --max-allowed-energy=0 --samc-t0 1e5 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-samc-38-1e5-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
# rq run -R --max-output 10 -J 'lj-samc-38-1e6-bin001' target/release/lj-cluster --N 38 --max-allowed-energy=0 --samc-t0 1e6 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-samc-38-1e6-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
# rq run -R --max-output 10 -J 'lj-samc-38-1e7-bin001' target/release/lj-cluster --N 38 --max-allowed-energy=0 --samc-t0 1e7 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-samc-38-1e7-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3

rq run -R --max-output 20 -J 'lj-wl-38-bin001' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --wl --translation-scale 0.05 --energy-bin 0.01 --save-as lj-wl-38-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3

# BIN WIDTH TESTS with 1/t-WL and SAD {0.1, 0.05, 0.01, 0.005, 0.001}
rq run -R --max-output 20 -J 'lj-sad-38-bin01' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.1 --save-as lj-sad-38-bin01.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-38-bin005' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.05 --save-as lj-sad-38-bin005.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-38-bin001' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-38-bin0005' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.005 --save-as lj-sad-38-bin0005.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-38-bin0001' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.001 --save-as lj-sad-38-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3

rq run -R --max-output 20 -J 'lj-inv-t-wl-38-bin01' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.1 --save-as lj-inv-t-wl-38-bin01.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-inv-t-wl-38-bin005' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.05 --save-as lj-inv-t-wl-38-bin005.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-inv-t-wl-38-bin001' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-38-bin001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-inv-t-wl-38-bin0005' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.005 --save-as lj-inv-t-wl-38-bin0005.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-inv-t-wl-38-bin0001' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.001 --save-as lj-inv-t-wl-38-bin0001.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3

# Rc TESTS with 1/t-WL and SAD {2, 3, 4, 5}
# rq run -R --max-output 20 -J 'lj-sad-38-bin001-r2' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001-r2.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2
# rq run -R --max-output 20 -J 'lj-sad-38-bin001-r4' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001-r4.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 4
rq run -R --max-output 20 -J 'lj-sad-38-bin001-r5' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001-r5.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 5
rq run -R --max-output 20 -J 'lj-sad-38-bin001-r6' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001-r6.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 6
rq run -R --max-output 20 -J 'lj-sad-38-bin001-r7' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001-r7.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 7

rq run -R --max-output 20 -J 'lj-sad-38-bin001-r10' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001-r10.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 10
rq run -R --max-output 20 -J 'lj-sad-38-bin001-r20' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001-r20.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 20
rq run -R --max-output 20 -J 'lj-sad-38-bin001-r30' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.05 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38-bin001-r30.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 30
#
# rq run -R --max-output 20 -J 'lj-inv-t-wl-38-bin001-r2' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-38-bin001-r2.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 2
# rq run -R --max-output 20 -J 'lj-inv-t-wl-38-bin001-r4' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-38-bin001-r4.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 4
# rq run -R --max-output 20 -J 'lj-inv-t-wl-38-bin001-r5' target/release/lj-cluster --N 38 --max-allowed-energy=0 --min-allowed-energy=-173 --Inv-t-WL --translation-scale 0.05 --energy-bin 0.01 --save-as lj-inv-t-wl-38-bin001-r5.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 5

# RUN THE SIMULATIONS THAT DAVID AND I TALKED ABOUT
# SAD FOR THE SIZES: N13 N19 N31 N38 N55 N75 N104 N135 N147

# Tmin = 0.01 energy-bin = 0.01 radius = 3Sigma

rq run -R --max-output 20 -J 'lj-sad-13' target/release/lj-cluster --N 13 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-13.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-19' target/release/lj-cluster --N 19 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-19.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-31' target/release/lj-cluster --N 31 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-31.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-38' target/release/lj-cluster --N 38 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-38.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-55' target/release/lj-cluster --N 55 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-55.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-75' target/release/lj-cluster --N 75 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-75.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-104' target/release/lj-cluster --N 104 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-104.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-135' target/release/lj-cluster --N 135 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-135.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3
rq run -R --max-output 20 -J 'lj-sad-147' target/release/lj-cluster --N 147 --max-allowed-energy=0 --sad-min-T 0.01 --translation-scale 0.05 --energy-bin 0.01 --save-as lj-sad-147.yaml --movie-time '10^(1/8)' --save-time 0.5 --radius 3

# LJ GRAND
 rq run -R --max-output 20 -J 'lj-grand' target/release/lj-cluster  --translation-scale 0.05 --radius 10 --max-N 13 --movie-time '10^(1/4)' --save-time 1 --save-as lj-grand.yaml --energy-bin 0.01 --samc-t0 1e7 --max-allowed-energy 0 --movie-name lj-grand.movie

rq run -R --max-output 20 -J 'lj-grand' --release --bin lj-grand -- --addremove-probability 0.1 --translation-scale 0.05 --radius 10 --max-N 13 --movie-time '10^(1/4)' --save-time 1 --save-as lj-grand.yaml --energy-bin 0.01 --samc-t0 1e7 --max-allowed-energy 0 --movie-name lj-grand.movie
