#!/bin/sh

#SBATCH --mail-type=fail
#SBATCH --ntasks=1
#SBATCH --mem=10000
#SBATCH --time=05:00:00
#SBATCH --job-name=extract_sing
#SBATCH --array=1-22
#SBATCH --output="slurm/singletons-%A_%a.out"
#SBATCH --error="slurm/slurmJob-%A_%a.err"

in_dir="../data/hgdp/vcf"
out_dir="../data/hgdp"

bcftools view -v snps ${in_dir}/hgdp_wgs.20190516.full.chr${SLURM_ARRAY_TASK_ID}.vcf.gz | vcftools --vcf - --singletons --out ${out_dir}/chr${SLURM_ARRAY_TASK_ID}_sing

