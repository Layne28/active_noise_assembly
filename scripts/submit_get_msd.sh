#!/bin/bash
#Submit msd calculation

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=1:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

folder=$1
tmax=$2

run_dir="/home/laynefrechette/lj_assembly/scripts"

python3 $run_dir/get_msd.py $folder $tmax

