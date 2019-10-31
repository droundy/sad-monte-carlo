set -ev

# First we need to build
cargo build --release

if [ "$HOSTNAME" = darcy ]; then
    echo Jobs for my laptop

    # rq run -R -J 'big-lj-benchmark-0.001-low' target/release/lj-cluster --N 31 --max-allowed-energy=-110 --min-allowed-energy=-133.579 --translation-scale 0.03 --movie-time '10^(1/8)' --save-time 2 --radius 5 --energy-bin 0.001 --Inv-t-WL --save-as big-lj-benchmark-0.001-low.yaml

    # rq run -R -J 'big-lj-sad-0.01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --translation-scale 0.05 --movie-time '10^(1/8)' --save-time 2 --radius 5 --energy-bin 0.01 --save-as big-lj-sad-0.01.yaml --sad-min-T 0.005

else
    echo Jobs for the cluster

    ## ------ RUN ON A CONSTRAINED ENERGY RANGE ------ ## # The groundstate is -133.586422
    rq run -R -J 'big-lj-sad-0.01-minT-0.01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --translation-scale 0.05 --movie-time '10^(1/8)' --save-time 2 --radius 5 --energy-bin 0.01 --save-as big-lj-sad-0.01-minT-0.01.yaml --sad-min-T 0.01

    # rq run --max-output 30 -R -J 'big-lj-sad-0.01-minT-0.01-bigmax' target/release/lj-cluster --N 31 --max-allowed-energy=200 --translation-scale 0.05 --movie-time '10^(1/8)' --movie-start 1e8 --save-time 2 --radius 5 --energy-bin 0.01 --save-as big-lj-sad-0.01-minT-0.01-bigmax.yaml --sad-min-T 0.01

    # rq run -R -J 'big-lj-sad-0.01-minT-0.02' target/release/lj-cluster --N 31 --max-allowed-energy=0 --translation-scale 0.05 --movie-time '10^(1/8)' --save-time 2 --radius 5 --energy-bin 0.01 --save-as big-lj-sad-0.01-minT-0.02.yaml --sad-min-T 0.02

    # rq run --max-output 30 -R -J 'big-lj-sad-0.001' target/release/lj-cluster --N 31 --max-allowed-energy=0 --translation-scale 0.05 --movie-time '10^(1/8)' --save-time 2 --radius 5 --energy-bin 0.001 --save-as big-lj-sad-0.001.yaml --sad-min-T 0.005

    # rq run -R -J 'big-lj-benchmark-0.001' target/release/lj-cluster --N 31 --max-allowed-energy=-110 --min-allowed-energy=-133.53 --translation-scale 0.03 --movie-time '10^(1/8)' --save-time 2 --radius 5 --energy-bin 0.001 --Inv-t-WL --save-as big-lj-benchmark-0.001.yaml

    # rq run -R -J 'big-lj-inv-0.01' target/release/lj-cluster --N 31 --max-allowed-energy=0 --min-allowed-energy=-133.53 --translation-scale 0.05 --movie-time '10^(1/8)' --save-time 2 --radius 5 --energy-bin 0.01 --Inv-t-WL --save-as big-lj-inv-0.01.yaml

    # ## ------ RUN ON A CONSTRAINED ENERGY RANGE ------ ## # The groundstate is -133.586422
    # rq run -R -J 'big-lj-benchmark' target/release/lj-cluster --N 31 --max-allowed-energy=-110 --min-allowed-energy=-133.53 --translation-scale 0.03 --movie-time '10^(1/8)' --save-time 2 --radius 5 --energy-bin 0.01 --Inv-t-WL --save-as big-lj-benchmark.yaml

fi
