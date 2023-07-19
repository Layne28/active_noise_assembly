#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4

#phis=(0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80)
phis=(0.00)
kTs=(0.0 0.1 0.3 0.5 1.0)
vas=(0.00 0.03 0.10 0.30 1.00)
taus=(0.10 1.00 10.00)
lambdas=(0.25 1.00 5.00 10.00)
#lambdas=(0.25)

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
                        python3 scripts/write_input_file.py input_files --va $va --tau $tau --Lambda $lambda --kT $kT --phi $phi --Lx $Lx --Ly $Ly --nx $nx --ny $ny --particles_freq 50 --thermo_freq 50 --dt 0.002 --production_steps 500000 --particle_protocol zeros --output_dir "/work/laynefrechette/lj-assembly-results/kT=${kT}/phi=${phi}/va=${va}/tau=${tau}/lambda=${lambda}/Lx=${Lx}_Ly=${Ly}/nx=${nx}_ny=${ny}"
                done
            done
        done
    done
done
