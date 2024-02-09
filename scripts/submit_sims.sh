#!/bin/bash
#Submit assembly simulations

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=96:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

in_file=$1
seed=$2
seed_file=$3

echo "Running assembly simulation with input file: $in_file"

run_dir="${HOME}/.local/bin"
echo ${run_dir}

$run_dir/ParticleSimTools $in_file $seed $seed_file

