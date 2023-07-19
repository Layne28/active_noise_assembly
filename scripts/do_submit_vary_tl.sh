#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
va=$5
kT=$6
potential=$7
nseed=$8

taus=(0.01 0.10 0.30 1.00 10.00)
#lambdas=(0.50)
lambdas=(0.50 1.00 2.00 3.00 5.00 10.00)
#phis=(0.20 0.40 0.60)
phis=(0.40)

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
                sbatch -J "${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}" -o "log/${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}.o%j" -e "log/${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}.e%j" scripts/submit_sims.sh "/home/laynefrechette/active_noise_assembly/input_files/${potential}_active_noise_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}.in" $i $input
            done
        done
    done
done
