#!/bin/bash

module load share_modules/ANACONDA/5.3_py3

phi=$1
rc=$2

Lx=200.0
Ly=200.0
nx=400
ny=400

kTs=(0.0)
vas=(1.00)
taus=(0.01 0.10 0.30 1.00 10.00)
lambdas=(0.50 1.00 2.00 3.00 5.00 10.00)
potential="wca"

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for kT in "${kTs[@]}"
        do
            for va in "${vas[@]}"
            do
                input=/scratch0/laynefrechette/active-noise-assembly-results/${potential}/kT\=${kT}/phi\=${phi}/va\=${va}/tau\=${tau}/lambda\=${lambda}/Lx\=${Lx}_Ly\=${Ly}/nx\=${nx}_ny\=${ny}/seed\=1/prod
                sbatch /home/laynefrechette/active-noise-assembly/scripts/submit_get_cluster_hist.sh $input $rc
            done
        done
    done
done
