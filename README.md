# LSCI_pipeline
Code to reproduce the analyses of (...)

## Data Prep

### Downloading data

#### 1kGP VCF Files

In the directory `data/1kgp`:

Run the following to download the 1kGP 30x deep-sequence VCF files:

```
for i in `seq 1 22`; do
  wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.annotated.vcf.gz"
done
```

We will also want to grab subject IDs and group them by their 1kGP superpopulation:

```
curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index |\
  awk -F"\t" 'NR > 24 {print($10"\t"$11)}' > metadata/subjects_populations.tsv

curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv | head -n -3 > metadata/20131219.populations.tsv
```

The script `src/append_1kgp_superpop.R` appends the super-population code to the sample list we downloaded above. The script can be run from the root directory of this project:

```
Rscript src/append_1kgp_superpop.R
```

#### HGDP Data

In the directory `data/hgdp`:

Run the following to download the HGDP VCF files (or download the files by hand from [Sanger's server](https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/)):

```
for i in `seq 1 22`; do
  wget "ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr${i}.vcf.gz"
done
```

We will also need the metadata for the samples:

```
curl https://ngs.sanger.ac.uk//production/hgdp/hgdp_wgs.20190516/metadata/hgdp_wgs.20190516.metadata.txt |\
    awk 'NR>1 {print($1"\t"$9)}' |\
    awk '{print($1) >> ("metadata/" $2 ".txt")}'
```

#### Reference Genome

From the directory `data/reference`, run the following to download human reference genome GRCh38:

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

See [this blog post](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) for why we chose this particular fasta

#### dbSNP 155

We will use dbSNP 155 to mask sites in the genome known to be variable.

```
# from the directory data/dbSNP
wget ftp://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb
bigBedToBed dbSnp155Common.bb dbSnp155Common.bed
# we only need the first three columns; subsetting yields a much smaller file size
awk '{print($1"\t"$2"\t"$3)}' dbSnp155Common.bed > dbSnp155Common_reduced.bed
rm dbSnp155Common.bb dbSnp155Common.bed
```

We use `bedtools maskfasta` to mask the variant sites in the reference genome:

```
# from the root directory of the project
bedtools maskfasta -fi data/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bed data/dbSNP/dbSnp155Common_reduced.bed -fo data/ref_genome/hg38_masked.fasta
```

### Extracting 1kGP Singletons

The `src/extract_singletons_*.sh` files are example slurm jobs to extract singletons for each of the 1kGP super-populations. We suggest using SLURM or some other workload manager, as these commands do take some time to run. The chain of commands for each super-population are also included below (directory references assume commands are being run from the `src` directory):

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
chromosome=22
population="EUROPE"
bcftools view -S metadata/${population}.txt -v snps data/hgdp/hgdp_wgs.20190516.full.chr${chromosome}.vcf.gz |\
  vcftools --vcf - --singletons --out output/singletons/hgdp/${population}/chr${chromosome}_sing
```


### Subtype Specific Files

The singleton extraction steps above yield a single file per chromosome. We rearrange the singletons into a single file per mutation subtype with the following shell commands:

```
for i in `seq 1 22`; do
echo $i
awk '{if(($3 == "A" && $4 == "C")|| ($3 == "T" && $4 == "G")) print($1"\t"$2"\t"$3) >> "AT_CG.txt"; 
  else if(($3 == "A" && $4 == "G")|| ($3 == "T" && $4 == "C")) print($1"\t"$2"\t"$3) >> "AT_GC.txt"; 
  else if(($3 == "A" && $4 == "T")|| ($3 == "T" && $4 == "A"))print($1"\t"$2"\t"$3) >> "AT_TA.txt" ;
  else if(($3 == "C" && $4 == "A")|| ($3 == "G" && $4 == "T"))print($1"\t"$2"\t"$3) >> "all_GC_TA.txt" ;
  else if(($3 == "C" && $4 == "G")|| ($3 == "G" && $4 == "C"))print($1"\t"$2"\t"$3) >> "all_GC_CG.txt";
  else if(($3 == "C" && $4 == "T")|| ($3 == "G" && $4 == "A"))print($1"\t"$2"\t"$3) >> "all_GC_AT.txt"; }' chr$i.txt 
done 
```

This needs to be run in `/output/hgdp_singletons` and in all subdirectories of `output/singletons`.

### Identify CpG and non CpGs

In our analyses, we consided CpG and non-CpG subtypes as distinct. To do this, we need to identify which GC_NN singletons are CpGs and which are not. The script `src/annotate_cpg.py` performs this task, and can be run as follows:

```
pop="SAS"
python src/annotate_cpg.py -s output/singletons/${pop}/all_GC_AT.txt -r data/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o output/singletons/${pop}/all_GC_AT_cpg.txt

python src/annotate_cpg.py -s output/singletons/${pop}/all_GC_TA.txt -r data/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o output/singletons/${pop}/all_GC_TA_cpg.txt

python src/annotate_cpg.py -s output/singletons/${pop}/all_GC_CG.txt -r data/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o output/singletons/${pop}/all_GC_CG_cpg.txt
```

We the generate separate CpG and Non-CpG C>N subtype files:

```
pop="SAS"
cd "output/singletons/$pop"

awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3) >> "cpg_GC_AT.txt";
  else if($4==0)print("chr"$1"\t"$2"\t"$3) >> "GC_AT.txt"; }' all_GC_AT_cpg.txt 

awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3) >> "cpg_GC_TA.txt";
  else if($4==0)print("chr"$1"\t"$2"\t"$3) >> "GC_TA.txt"; }' all_GC_TA_cpg.txt 
  
awk -F, '{if($4==1)print("chr"$1"\t"$2"\t"$3) >> "cpg_GC_CG.txt";
  else if($4==0)print("chr"$1"\t"$2"\t"$3) >> "GC_CG.txt"; }' all_GC_CG_cpg.txt 

cd ../../../
```

### Sample Control Observations

Run for each super-population. In our analyses, we sampled 5 control observations for each singleton.

```
# Sample controls
n_control=1
pop="AFR"
ref_genome="data/ref_genome/hg38_masked.fasta"

python src/sample_control.py -s output/singletons/${pop}/AT_CG.txt -f ${ref_genome} -o output/controls/${pop}/AT_CG.csv -t "AT_CG" -n ${n_control}
python src/sample_control.py -s output/singletons/${pop}/AT_GC.txt -f ${ref_genome} -o output/controls/${pop}/AT_GC.csv -t "AT_GC" -n ${n_control}
python src/sample_control.py -s output/singletons/${pop}/AT_TA.txt -f ${ref_genome} -o output/controls/${pop}/AT_TA.csv -t "AT_TA" -n ${n_control}

python src/sample_control.py -s output/singletons/${pop}/GC_AT.txt -f ${ref_genome} -o output/controls/${pop}/GC_AT.csv -t "GC_AT" -n ${n_control}
python src/sample_control.py -s output/singletons/${pop}/GC_TA.txt -f ${ref_genome} -o output/controls/${pop}/GC_TA.csv -t "GC_TA" -n ${n_control}
python src/sample_control.py -s output/singletons/${pop}/GC_CG.txt -f ${ref_genome} -o output/controls/${pop}/GC_CG.csv -t "GC_CG" -n ${n_control}

python src/sample_control.py -s output/singletons/${pop}/cpg_GC_AT.txt -f ${ref_genome} -o output/controls/${pop}/cpg_GC_AT.csv -t "cpg_GC_AT" -n ${n_control}
python src/sample_control.py -s output/singletons/${pop}/cpg_GC_TA.txt -f ${ref_genome} -o output/controls/${pop}/cpg_GC_TA.csv -t "cpg_GC_TA" -n ${n_control}
python src/sample_control.py -s output/singletons/${pop}/cpg_GC_CG.txt -f ${ref_genome} -o output/controls/${pop}/cpg_GC_CG.csv -t "cpg_GC_CG" -n ${n_control}
```

### Generating Position Files

These files are the input for the model fitting code. We need to separate A > N from T > N (and similarly, C > N and G > N) since relative positions need to be reversed for the T > N (G > N) relative to the A > N (C > N) positions.

#### Singletons

*Note: commands below works in zsh -- for command would need to be rewritten for bash*

```
pop="AFR"
cd output/singletons/${pop}

typeset -a subtypes
subtypes=("AT_CG" "AT_GC" "AT_TA")
for i ("$subtypes[@]"); do
  print $i
  awk -v var=$i '{if($3 == "A")print($2) >> "pos_files/"var"_"substr($1,4)".txt"; 
    else if($3 == "T")print($2) >> "pos_files/"var"_rev_"substr($1,4)".txt"; }' $i.txt
done

subtypes=("GC_AT" "GC_TA" "GC_CG" "cpg_GC_AT" "cpg_GC_TA" "cpg_GC_CG")
for i ("$subtypes[@]"); do
  print $i
  awk -v var=$i '{if($3 == "C")print($2) >> "pos_files/"var"_"substr($1,4)".txt"; 
    else if($3 == "G")print($2) >> "pos_files/"var"_rev_"substr($1,4)".txt"; }' $i.txt
done

cd ../../../
```

#### Controls

```
pop="AFR"
cd output/controls/${pop}

typeset -a subtypes
subtypes=("AT_CG" "AT_GC" "AT_TA")
for i ("$subtypes[@]"); do
  print $i
  awk -v var=$i -F, '{if($4 == "A")print($2) >> "pos_files/"var"_"substr($1,4)".txt"; 
    else if($4 == "T")print($2) >> "pos_files/"var"_rev_"substr($1,4)".txt"; }' $i.csv
done

subtypes=("GC_AT" "GC_TA" "GC_CG" "cpg_GC_AT" "cpg_GC_TA" "cpg_GC_CG")
for i ("$subtypes[@]"); do
  print $i
  awk -v var=$i -F, '{if($4 == "C")print($2) >> "pos_files/"var"_"substr($1,4)".txt"; 
    else if($4 == "G")print($2) >> "pos_files/"var"_rev_"substr($1,4)".txt"; }' $i.csv
done

cd ../../../
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

## Fitting Models

### Single Position Models

```
python src/single_pos.py 
```


