#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
tau=$5 #default to 10
lambda=$6 #default to 10

#phis=(0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80)
phis=(0.00 0.20 0.40 0.60)
#phis=(0.4)
kTs=(0.0 0.1 0.3 0.5 1.0)
vas=(0.00 0.03 0.10 0.30 1.00 3.00 10.00)

for phi in "${phis[@]}"
do
    for kT in "${kTs[@]}"
    do
        for va in "${vas[@]}"
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
