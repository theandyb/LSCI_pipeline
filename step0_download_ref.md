## Reference Genome

## 1kGP VCF Files

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

## HGDP Data

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
    awk '{print>"metadata/"$2".txt"}'
```


## Reference Genome

From the directory `data/reference`, run the following to download human reference genome GRCh38:

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

See [this blog post](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) for why we chose this particular fasta

## dbSNP 155

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
samtools faidx data/ref_genome/hg38_masked.fasta
```