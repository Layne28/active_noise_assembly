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
randomization_steps=40000
equil_steps=400000
production_steps=4000000

taus=(0.010000 0.100000 0.300000 1.000000 3.000000 10.000000)
lambdas=(0.500000 1.000000 2.000000 3.000000 5.000000 10.000000)
phis=(0.100000 0.200000 0.400000 0.600000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for phi in "${phis[@]}"
        do
            python3 scripts/write_input_file.py input_files --va $va --tau $tau --Lambda $lambda --kT $kT --phi $phi --Lx $Lx --Ly $Ly --nx $nx --ny $ny --particles_freq $freq --thermo_freq $freq --dt $dt --equil_steps=$equil_steps --production_steps=$production_steps --particle_protocol "triangular_lattice" --nonbonded_potential_type $potential --output_dir "/${basedir}/laynefrechette/active-noise-assembly-results/${potential}/kT=${kT}/phi=${phi}/va=${va}/tau=${tau}/lambda=${lambda}/Lx=${Lx}_Ly=${Ly}/nx=${nx}_ny=${ny}"
        done
    done
done
