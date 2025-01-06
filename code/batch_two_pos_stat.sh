#!/bin/sh

# We assume that the job is submitted from the root directory

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem-per-cpu=2GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --time 10:00:00
#SBATCH --job-name=sas_CpG_T
#SBATCH --requeue
#SBATCH -p main
#SBATCH -e output/slurm/sas_CpG_T-%A_%a.err
#SBATCH -o output/slurm/sas_CpG_T-%A_%a.out

subtype="cpg_GC_AT"
pop="SAS"

echo "Running single position models for subtype ${subtype} in population ${pop}"

python code/two_position_model.py -p ${pop} -t ${subtype} -s ""
