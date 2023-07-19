#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
va=$5 
kT=$6
potential=$7
basedir=$8 #work or scratch0

dt=0.00025
freq=800
equil_steps=40000
production_steps=4000000

taus=(0.01 0.10 0.30 1.00 10.00)
lambdas=(0.50 1.00 2.00 3.00 5.00 10.00)
phis=(0.20 0.40 0.60)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for phi in "${phis[@]}"
        do
            python3 scripts/write_input_file.py input_files --va $va --tau $tau --Lambda $lambda --kT $kT --phi $phi --Lx $Lx --Ly $Ly --nx $nx --ny $ny --particles_freq $freq --thermo_freq $freq --dt $dt --equil_steps=$equil_steps --production_steps=$production_steps --particle_protocol "lattice" --potential_type $potential --output_dir "/${basedir}/laynefrechette/active-noise-assembly-results/${potential}/kT=${kT}/phi=${phi}/va=${va}/tau=${tau}/lambda=${lambda}/Lx=${Lx}_Ly=${Ly}/nx=${nx}_ny=${ny}"
        done
    done
done
