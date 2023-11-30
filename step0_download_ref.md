## Reference Genome

## 1kGP VCF Files

In the directory `data/1kgp`:

Run the following to download the 1kGP 30x deep-sequence VCF files:

```
for i in `seq 1 22`; do
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.vcf.gz
done
```

## HGDP Data

In the directory `data/hgdp`:

Run the following to download the HGDP VCF files (or download the files by hand from [Sanger's server](https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/)):

```
for i in `seq 1 22`; do
  wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr${i}.vcf.gz
done
```
