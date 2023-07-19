#!/bin/bash

Lx=$1
Ly=$2
nx=$3
ny=$4
tau=$5 #default to 10
lambda=$6 #default to 10
potential=$7
basedir=$8 #work or scratch0

dt=0.00025
freq=800
equil_steps=40000
production_steps=4000000
#phis=(0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80)
phis=(0.20 0.40 0.60)
kTs=(0.0 0.1 0.3 0.5 1.0)
vas=(0.00 0.03 0.10 0.30 1.00)
#vas=(0.00 0.03)

for phi in "${phis[@]}"
do
    for kT in "${kTs[@]}"
    do
        for va in "${vas[@]}"
        do
            python3 scripts/write_input_file.py input_files --va $va --tau $tau --Lambda $lambda --kT $kT --phi $phi --Lx $Lx --Ly $Ly --nx $nx --ny $ny --particles_freq $freq --thermo_freq $freq --dt $dt --equil_steps=$equil_steps --production_steps=$production_steps --particle_protocol "lattice" --potential_type $potential --output_dir "/${basedir}/laynefrechette/active-noise-assembly-results/${potential}/kT=${kT}/phi=${phi}/va=${va}/tau=${tau}/lambda=${lambda}/Lx=${Lx}_Ly=${Ly}/nx=${nx}_ny=${ny}"
        done
    done
done
