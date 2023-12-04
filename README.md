# LSCI_pipeline
Code to reproduce the analyses of (...)

## Data Prep

Simple commands for downloading and processing the reference data are available in the `step0_download_ref.md`. Note that the analysis steps assume that a particular directory structure is maintained.

### Extracting 1kGP Singletons

The `src/extract_singletons_*.sh` files are example slurm jobs to extract singletons for each of the 1kGP super-populations. We suggest using SLURM or some other workload manager, as these commands do take some time to run. The chain of commands for each superpopulation are also included below (directory references assume commands are being run from the `src` directory):

#### AFR Singletons

```
for i in `seq 1 22`; do
    bcftools view -i "%FILTER=='PASS'" -v snps ../data/reference/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.annotated.vcf.gz |\
    bcftools view -i "(AC_AFR_unrel == 1 | (AC_Hom_AFR_unrel == 2 & AC_Het_AFR_unrel ==0)) & AC_EUR_unrel == 0 & AC_AMR_unrel == 0 & AC_EAS_unrel == 0 & AC_SAS_unrel == 0" |\
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" > ../output/AFR/chr${SLURM_ARRAY_TASK_ID}.txt
done
```

#### AMR Singletons

```
for i in `seq 1 22`; do
    bcftools view -i "%FILTER=='PASS'" -v snps ../data/reference/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.annotated.vcf.gz |\
    bcftools view -i "(AC_AFR_unrel == 1 | (AC_Hom_AMR_unrel == 2 & AC_Het_AMR_unrel ==0)) & AC_EUR_unrel == 0 & AC_AFR_unrel == 0 & AC_EAS_unrel == 0 & AC_SAS_unrel == 0" |\
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" > ../output/AMR/chr${SLURM_ARRAY_TASK_ID}.txt
done
```

#### EAS Singletons

```
for i in `seq 1 22`; do
    bcftools view -i "%FILTER=='PASS'" -v snps ../data/reference/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.annotated.vcf.gz |\
    bcftools view -i "(AC_EAS_unrel == 1 | (AC_Hom_EAS_unrel == 2 & AC_Het_EAS_unrel ==0)) & AC_EUR_unrel == 0 & AC_AFR_unrel == 0 & AC_AMR_unrel == 0 & AC_SAS_unrel == 0" |\
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" > ../output/EAS/chr${SLURM_ARRAY_TASK_ID}.txt
done
```

#### EUR Singletons

```
for i in `seq 1 22`; do
    bcftools view -i "%FILTER=='PASS'" -v snps ../data/reference/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.annotated.vcf.gz |\
    bcftools view -i "(AC_EUR_unrel == 1 | (AC_Hom_EUR_unrel == 2 & AC_Het_EUR_unrel ==0)) & AC_EAS_unrel == 0 & AC_AFR_unrel == 0 & AC_AMR_unrel == 0 & AC_SAS_unrel == 0" |\
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" > ../output/EUR/chr${SLURM_ARRAY_TASK_ID}.txt
done
```

#### SAS Singletons

```
for i in `seq 1 22`; do
    bcftools view -i "%FILTER=='PASS'" -v snps ../data/reference/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.annotated.vcf.gz |\
    bcftools view -i "(AC_SAS_unrel == 1 | (AC_Hom_SAS_unrel == 2 & AC_Het_SAS_unrel ==0)) & AC_EAS_unrel == 0 & AC_AFR_unrel == 0 & AC_AMR_unrel == 0 & AC_EUR_unrel == 0" |\
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" > ../output/SAS/chr${SLURM_ARRAY_TASK_ID}.txt
done
```

### Extracting HGDP Singletons


### Annotating 1kGP Singletons

The code for annotating each singleton with its subtype, nucleotide motif, etc is in the script `src/annotate_singletons.py`. The three argmuents which need to be provided are:

1. The singleton file (`-s`)
2. The location of the reference genome (`-r`)
3. The intended output (`-o`)




