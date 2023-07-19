#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
nseed=$5

taus=(0.10 1.00 10.00)
#lambdas=(1.00 2.00 5.00 10.00)
lambdas=(0.25)

phis=(0.00)
kTs=(0.0 0.1 0.3 0.5 1.0)
vas=(0.03 0.10 0.30 1.00)

for phi in "${phis[@]}"
do
    for kT in "${kTs[@]}"
    do
        for va in "${vas[@]}"
        do
            for tau in "${taus[@]}"
            do
                for lambda in "${lambdas[@]}"
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
    done
done
