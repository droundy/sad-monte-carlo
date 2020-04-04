set -ev

# First we need to build
cargo build --release

rq run -R --max-output 20 -J 'ising-wl-32-minGamma-s1' target/release/ising-mc --N 32 --seed=1 --save-as ising-wl-32-minGamma-s1.yaml --movie-time '10^(1/8)' --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-2048 --max-allowed-energy=50 --wl --wl-min-gamma=1e-4

rq run -R --max-output 20 -J 'ising-wl-32-minGamma-s2' target/release/ising-mc --N 32 --seed=2 --save-as ising-wl-32-minGamma-s2.yaml --movie-time '10^(1/8)' --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-2048 --max-allowed-energy=50 --wl --wl-min-gamma=1e-4

rq run -R --max-output 20 -J 'ising-wl-32-minGamma-s3' target/release/ising-mc --N 32 --seed=3 --save-as ising-wl-32-minGamma-s3.yaml --movie-time '10^(1/8)' --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-2048 --max-allowed-energy=50 --wl --wl-min-gamma=1e-4

rq run -R --max-output 20 -J 'ising-wl-32-minGamma-s4' target/release/ising-mc --N 32 --seed=4 --save-as ising-wl-32-minGamma-s4.yaml --movie-time '10^(1/8)' --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-2048 --max-allowed-energy=50 --wl --wl-min-gamma=1e-4

rq run -R --max-output 20 -J 'ising-wl-32-minGamma-s5' target/release/ising-mc --N 32 --seed=5 --save-as ising-wl-32-minGamma-s5.yaml --movie-time '10^(1/8)' --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-2048 --max-allowed-energy=50 --wl --wl-min-gamma=1e-4

rq run -R --max-output 20 -J 'ising-wl-32-minGamma-s6' target/release/ising-mc --N 32 --seed=6 --save-as ising-wl-32-minGamma-s6.yaml --movie-time '10^(1/8)' --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-2048 --max-allowed-energy=50 --wl --wl-min-gamma=1e-4

rq run -R --max-output 20 -J 'ising-wl-32-minGamma-s7' target/release/ising-mc --N 32 --seed=7 --save-as ising-wl-32-minGamma-s7.yaml --movie-time '10^(1/8)' --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-2048 --max-allowed-energy=50 --wl --wl-min-gamma=1e-4

rq run -R --max-output 20 -J 'ising-wl-32-minGamma-s8' target/release/ising-mc --N 32 --seed=8 --save-as ising-wl-32-minGamma-s8.yaml --movie-time '10^(1/8)' --translation-scale 0.05 --save-time 0.5 --min-allowed-energy=-2048 --max-allowed-energy=50 --wl --wl-min-gamma=1e-4