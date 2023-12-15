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

We use vcftool's `singletons` option to extract singletons from the HGDP dataset. While this command can be run from the terminal, we suggest using SLURM as this task can take some time.


```

```


### Subtype Specific Files

The singleton extraction steps above yield a single file per chromosome. We rearrange the singletons into a single file per mutation subtype with the following shell commands:

```
for i in `seq 1 22`; do
echo $i
awk '{if(($3 == "A" && $4 == "C")|| ($3 == "T" && $4 == "G"))print($1"\t"$2"\t"$3)}' chr$i.txt >> AT_CG.txt
awk '{if(($3 == "A" && $4 == "G")|| ($3 == "T" && $4 == "C"))print($1"\t"$2"\t"$3)}' chr$i.txt >> AT_GC.txt
awk '{if(($3 == "A" && $4 == "T")|| ($3 == "T" && $4 == "A"))print($1"\t"$2"\t"$3)}' chr$i.txt >> AT_TA.txt
awk '{if(($3 == "C" && $4 == "A")|| ($3 == "G" && $4 == "T"))print($1"\t"$2"\t"$3)}' chr$i.txt >> all_GC_TA.txt
awk '{if(($3 == "C" && $4 == "G")|| ($3 == "G" && $4 == "C"))print($1"\t"$2"\t"$3)}' chr$i.txt >> all_GC_CG.txt
awk '{if(($3 == "C" && $4 == "T")|| ($3 == "G" && $4 == "A"))print($1"\t"$2"\t"$3)}' chr$i.txt >> all_GC_AT.txt
done 
```

This needs to be run in `/output/hgdp_singletons` and in all subdirectories of `output/singletons`.

### Identify CpG and non CpGs

In our analyses, we consided CpG and non-CpG subtypes as distinct. To do this, we need to identify which GC_NN singletons are CpGs and which are not. The script `src/annotate_cpg.py` performs this task, and can be run as follows:

```
pop="AFR"
python src/annotate_cpg.py -s output/${pop}/all_GC_AT.txt -r data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o output/singletons/${pop}/all_GC_AT_cpg.txt

python src/annotate_cpg.py -s output/${pop}/all_GC_TA.txt -r data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o output/singletons/${pop}/all_GC_TA_cpg.txt

python src/annotate_cpg.py -s output/${pop}/all_GC_CG.txt -r data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o output/singletons/${pop}/all_GC_CG_cpg.txt
```

We the generate separate CpG and Non-CpG C>N subtype files:

```
pop="AFR"
awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3)}' all_GC_AT_cpg > cpg_GC_AT.txt
awk -F, '{if($4==0)print("chr"$1"\t"$2"\t"$3)}' all_GC_AT_cpg > GC_AT.txt
awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3)}' all_GC_TA_cpg > cpg_GC_TA.txt
awk -F, '{if($4==0)print("chr"$1"\t"$2"\t"$3)}' all_GC_TA_cpg > GC_TA.txt
awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3)}' all_GC_CG_cpg > cpg_GC_CG.txt
awk -F, '{if($4==0)print("chr"$1"\t"$2"\t"$3)}' all_GC_CG_cpg > GC_CG.txt
```


### Annotating 1kGP Singletons

This is an optional step which generates some additional information for the singletons. 

The code for annotating each singleton with its subtype, nucleotide motif, etc is in the script `src/annotate_singletons.py`. The three argmuents which need to be provided are:

1. The singleton file (`-s`)
2. The location of the reference genome (`-r`)
3. The intended output (`-o`)

From the root directory, an example command would be:

```
pop="AFR"
subtype="AT_CG"
python src/annotate_singletons.py -s output/singletons/${pop}/${subtype}.txt -o output/singletons/${pop}/annotated/${subtype}.txt -r data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```





