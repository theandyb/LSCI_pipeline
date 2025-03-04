#!/bin/sh

# We assume that the job is submitted from the root directory

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem-per-cpu=2GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --time 12:00:00
#SBATCH --job-name=single
#SBATCH --requeue
#SBATCH -p main
#SBATCH --array=0-53
#SBATCH -e output/slurm/singlePos-%A_%a.err
#SBATCH -o output/slurm/singlePos-%A_%a.out

VAR1_VALUES=("AFR" "AMR" "EAS" "EUR" "SAS" "ALL")
VAR2_VALUES=("AT_CG" "AT_GC" "AT_TA" "GC_AT" "GC_TA" "GC_CG" "cpg_GC_AT" "cpg_GC_TA" "cpg_GC_CG")

INDEX=$SLURM_ARRAY_TASK_ID

VAR1_INDEX=$(( INDEX / 9 ))
VAR2_INDEX=$(( INDEX % 9 ))

VAR1=${VAR1_VALUES[$VAR1_INDEX]}
VAR2=${VAR2_VALUES[$VAR2_INDEX]}

echo "Running single position models for subtype ${VAR2} in population ${VAR1}"

python code/single_position_models_posOnly.py -p ${VAR1} -t ${VAR2} -o output
