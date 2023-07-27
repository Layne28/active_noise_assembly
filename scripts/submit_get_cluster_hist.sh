#!/bin/bash
#Submit cluster calculation

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=1:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

folder=$1

run_dir="/home/laynefrechette/active_noise_assembly/scripts"

python3 $run_dir/get_cluster_hist.py $folder

