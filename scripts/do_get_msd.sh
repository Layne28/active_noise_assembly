#!/bin/bash

module load share_modules/ANACONDA/5.3_py3

Lx=50.0
Ly=50.0
nx=200
ny=200
phi=0.00

kTs=(0.0 0.1 0.3 0.5 1.0)
vas=(0.03 0.10 0.30 1.00)
taus=(0.10 1.00 10.00)
#lambdas=(1.00 5.00 10.00)
lambdas=(0.25)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for kT in "${kTs[@]}"
        do
            for va in "${vas[@]}"
            do
                input=/work/laynefrechette/lj-assembly-results/kT\=${kT}/phi\=${phi}/va\=${va}/tau\=${tau}/lambda\=${lambda}/Lx\=${Lx}_Ly\=${Ly}/nx\=${nx}_ny\=${ny}/seed\=1/prod
                sbatch /home/laynefrechette/lj_assembly/scripts/submit_get_msd.sh $input 100
                #python3 /home/laynefrechette/lj_assembly/scripts/get_msd.py $input 100
            done
        done
    done
done
