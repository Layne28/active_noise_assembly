#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
va=$5
kT=$6

taus=(0.01 0.10 0.30 1.00 10.00)
lambdas=(0.50 1.00 2.00 3.00 5.00 10.00)
phis=(0.20 0.40 0.60)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for phi in "${phis[@]}"
        do
            #filename="seeds/seeds_Nx=${Nx}_nx=${nx}_va=${va}_tau=${tau}_lambda=${lambda}.txt"
            filename="seeds/seeds_kT=${kT}_phi=${phi}_va=${va}_tau=${tau}_lambda=${lambda}_Lx=${Lx}_Ly=${Ly}_nx=${nx}_ny=${ny}.txt"
            touch $filename
            for i in {1..5}
            do
                printf "$RANDOM\n" >> $filename
            done
        done
    done
done
