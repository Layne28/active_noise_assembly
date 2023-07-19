#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
tau=$5
lambda=$6
nseed=$7

#taus=(0.01 0.10 0.30 1.00 10.00)
#lambdas=(0.01 0.10 0.30 1.00 2.00 5.00 10.00)

#phis=(0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80)
#phis=(0.40)
phis=(0.20 0.40 0.60)
#phis=(0.20)
kTs=(0.0 0.1 0.3 0.5 1.0)
#kTs=(0.0 0.1 0.3)
#vas=(0.00 0.10 0.30 1.00 3.00)
vas=(0.00 0.03 0.10 0.30 1.00)
#vas=(0.00)

for phi in "${phis[@]}"
do
    for kT in "${kTs[@]}"
    do
        for va in "${vas[@]}"
        do
            input="seeds/seeds_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}.txt"
            for (( i=1; i<=$nseed; i++ ))
            do
                echo $i
                sbatch -J "lj_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}" -o "log/lj_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}.o%j" -e "log/lj_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}_seed=${i}.e%j" scripts/submit_sims.sh "/home/laynefrechette/lj_assembly/input_files/lj_assembly_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}.in" $i $input
            done
        done
    done
done
