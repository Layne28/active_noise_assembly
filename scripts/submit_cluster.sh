#!/bin/bash
#Submit cluster calculation

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=10:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --exclude=compute-3-17

folder=$1
rcmult=$2

run_dir="/home/laynefrechette/active-noise-assembly/scripts"

echo $folder
python3 $run_dir/cluster.py $folder $rcmult

