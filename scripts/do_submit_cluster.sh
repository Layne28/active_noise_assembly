#!/bin/bash

module load share_modules/ANACONDA/5.3_py3

Lx=200.0
Ly=200.0
nx=400
ny=400
phi=$1

kTs=(0.0)
vas=(1.00)
taus=(0.01 0.10 0.30 1.00 10.00)
lambdas=(0.50 1.00 2.00 3.00 5.00 10.00)
potential="wca"
rcmult=$2

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for kT in "${kTs[@]}"
        do
            for va in "${vas[@]}"
            do
                input=/scratch0/laynefrechette/active-noise-assembly-results/${potential}/kT\=${kT}/phi\=${phi}/va\=${va}/tau\=${tau}/lambda\=${lambda}/Lx\=${Lx}_Ly\=${Ly}/nx\=${nx}_ny\=${ny}/seed\=1
                sbatch -o "log/cluster_${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}.o%j" -e "log/cluster_${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}.e%j" scripts/submit_cluster.sh $input $rcmult
            done
        done
    done
done
