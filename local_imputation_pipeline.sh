#!/bin/bash

##SNP calling using ATLAS on a 10 Mb region on chr20

atlas task=callNEW method=MLE infoFields=DP formatFields=GT,AD,DP,PL bam=NE1_chr20_0.05x.bam sites=1KG_chr20_10Mb_biallelic_maf0.1_sites.txt fasta=hg19chr.fasta verbose chr=chr20 out=NE1_005x_1KG_ATLASmle

##Genotype likelihood update: Beagle GL
java -Xmx50g -jar beagle.4.1.jar gl=NE1_005x_1KG_ATLASmle_MaximumLikelihood.vcf gprobs=true map=1KG_chr20_10Mb_biallelic_maf0.1.map ref=1KG_chr20_10Mb_biallelic_maf0.1.vcf window=2000 overlap=200 out=NE1_005x_1KG_ATLASmle_beagleGL20001KG

##Pre-imputation filter: GP >= 0.99
bgzip -c NE1_005x_1KG_ATLASmle_beagleGL20001KG.vcf.gz 
tabix -h NE1_005x_1KG_ATLASmle_beagleGL20001KG.vcf.gz
bcftools filter NE1_005x_1KG_ATLASmle_beagleGL20001KG.vcf.gz -i "FORMAT/GP >= 0.99" -o NE1_005x_1KG_ATLASmle_beagleGL20001KG_IGP099.vcf


##Genotype imputation: Beagle GT
java -Xmx50g -jar beagle.5.0.jar gt=NE1_005x_1KG_ATLASmle_beagleGL20001KG_IGP099.vcf gp=true map=1KG_chr20_10Mb_biallelic_maf0.1.map ref=1KG_chr20_10Mb_biallelic_maf0.1.vcf out=NE1_005x_1KG_ATLASmle_beagleGL20001KG_IGP099_beagleGT


##Post-imputation filter: GP >= 0.99
bgzip -c NE1_005x_1KG_ATLASmle_beagleGL20001KG_IGP099_beagleGT.vcf.gz
tabix -h NE1_005x_1KG_ATLASmle_beagleGL20001KG_IGP099_beagleGT.vcf.gz
bcftools filter NE1_005x_1KG_ATLASmle_beagleGL20001KG_IGP099_beagleGT -i "FORMAT/GP >= 0.99" -o NE1_005x_1KG_ATLASmle_beagleGL20001KG_IGP099_beagleGT_IIGP099.vcf

##Comparison with high-coverage NE1
java -Xmx50g -jar SnpSift.jar concordance -v NE1_highcoverage_chr20.vcf NE1_005x_1KG_ATLASmle_beagleGL20001KG_IGP099_beagleGT_IIGP099.vcf > NE1_comparison.txt
