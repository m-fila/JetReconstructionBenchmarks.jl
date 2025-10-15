#! /bin/sh
#
# Loop over parameters for generating benchmarks

# Abort on error
set -e

# Generate output directory - give this as the first argument to the script
output=${1:-"results/test-pp"}
mkdir -p $output

echo "Benchmarking outputs: $output"

# Input datafile list
inputs=data/events-summary-ee-pp.csv

# Iterate over backends
for backend in JetReconstruction Fastjet; do
    for strategy in N2Plain N2Tiled; do
        for radius in 0.2 0.4 1.0 1.5 2.0 3.0; do
            for algorithm in AntiKt CA Kt; do
                cmd="julia --project src/benchmark.jl --code $backend -R $radius -A $algorithm -S $strategy -m 4 --results $output $inputs"
                echo "Benchmark: $cmd"
                $cmd
            done
        done
    done
done
