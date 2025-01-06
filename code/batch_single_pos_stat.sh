#!/bin/sh

# We assume that the job is submitted from the root directory

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem-per-cpu=1GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --time 25:00:00
#SBATCH --job-name=all_C_T
#SBATCH --requeue
#SBATCH -p main
#SBATCH -e output/slurm/all_C_T-%A_%a.err
#SBATCH -o output/slurm/all_C_T-%A_%a.out

subtype="GC_AT"
pop="ALL"

echo "Running single position models for subtype ${subtype} in population ${pop}"

python code/single_position_models_posOnly.py -p ${pop} -t ${subtype} -s ""
