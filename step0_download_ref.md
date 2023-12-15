## Reference Genome

## 1kGP VCF Files

In the directory `data/1kgp`:

Run the following to download the 1kGP 30x deep-sequence VCF files:

```
for i in `seq 1 22`; do
  wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.vcf.gz"
done
```

We will also want to grab subject IDs and group them by their 1kGP superpopulation:


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