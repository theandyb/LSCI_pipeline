#!/bin/sh

# We assume that the job is submitted from the root directory

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem-per-cpu=1GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --time 00:30:00
#SBATCH --job-name=eas_A_C
#SBATCH --requeue
#SBATCH -p main
#SBATCH -e output/slurm/eas_A_C-%A_%a.err
#SBATCH -o output/slurm/eas_A_C-%A_%a.out

subtype="AT_CG"
pop="EAS"

echo "Running single position models for subtype ${subtype} in population ${pop}"

python code/single_position_models.py -p ${pop} -s ${subtype} -f data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o output
