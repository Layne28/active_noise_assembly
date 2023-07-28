#!/bin/bash

module load share_modules/ANACONDA/5.3_py3

Lx=100.0
Ly=100.0
nx=200
ny=200
phi=0.40

kTs=(0.0)
vas=(1.00)
taus=(0.01 0.10 0.30 1.00 10.00)
lambdas=(0.50 1.00 2.00 3.00 5.00 10.00)
potential="wca"
rcmult=$1

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for kT in "${kTs[@]}"
        do
            for va in "${vas[@]}"
            do
                input=/scratch0/laynefrechette/active-noise-assembly-results/${potential}/kT\=${kT}/phi\=${phi}/va\=${va}/tau\=${tau}/lambda\=${lambda}/Lx\=${Lx}_Ly\=${Ly}/nx\=${nx}_ny\=${ny}/seed\=1
                sbatch /home/laynefrechette/active_noise_assembly/scripts/submit_cluster.sh $input $rcmult
            done
        done
    done
done
