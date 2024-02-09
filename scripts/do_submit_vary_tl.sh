#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
va=$5
kT=$6
phi=$7
potential=$8
nseed=$9

#taus=(0.01 0.10 0.30 1.00 3.00 10.00)
taus=(0.300000 1.000000 3.000000 10.000000)
#lambdas=(0.50)
lambdas=(0.500000 1.000000 2.000000 3.000000 5.000000 10.000000)
#phis=(0.20 0.40 0.60)
phis=($phi)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for phi in "${phis[@]}"
        do
            input="seeds/seeds_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}.txt"
            for (( i=1; i<=$nseed; i++ ))
            do
                echo $i
                sbatch -J "${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}" -o "log/${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}.o%j" -e "log/${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}.e%j" scripts/submit_sims.sh "/home/laynefrechette/active-noise-assembly/input_files/${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}.in" $i $input
            done
        done
    done
done
