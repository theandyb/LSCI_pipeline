#!/bin/sh

# We assume that the job is submitted from the src/ directory

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=you@domain.tld
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --job-name=singletons_EAS
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/singleton_eas-%J.err
#SBATCH -o slurm/singleton_eas-%J.out

OUT_DIR="../output/singletons/EAS/"

bcftools view -i "%FILTER=='PASS'" -v snps ../data/reference/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${SLURM_ARRAY_TASK_ID}.recalibrated_variants.annotated.vcf.gz |\
bcftools view -i "(AC_EAS_unrel == 1 | (AC_Hom_EAS_unrel == 2 & AC_Het_EAS_unrel ==0)) & AC_EUR_unrel == 0 & AC_AFR_unrel == 0 & AC_AMR_unrel == 0 & AC_SAS_unrel == 0" |\
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" > $OUT_DIR/chr${SLURM_ARRAY_TASK_ID}.txt