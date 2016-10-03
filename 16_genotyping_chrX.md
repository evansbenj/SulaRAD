# Genotyping chrX is difficult

Background: my preliminary results suggest that there are some issues with genotyping chrX in males.  I have tried doing this with Haplotype caller using the default setting for ploidy (=2) and setting ploidy=1.  The former has problems because many heterozygous calls are made.  The latter has problems because the genotypes have a bias towards reference sites.

As a solution, Janet Kelso from the MPI suggested I genotype male chrX sites using the highest frequenty SNP.  This sounds like a great idea to me!  I plan to do this for the female as well so that each chrX is treated the same. I will write a script that outputs a tab delimited file from a diploid vcf file with genotyping based on the AD annotation (AD is allele depth).  For sites with an equal frequency of the ref and alt SNP, I will select the genotype randomly (and also keep track of how frequently this happens.

# First output chrX

I have some vcf files made at the NYGenome center here (on iqaluk):
```
/work/ben/2015_SulaRADtag/vcf-constitutional
```
```
-rw-rw-r-- 1 ben ben 6238245944 Sep  2 00:32 nemestrina-PM664.g.vcf.gz
-rw-rw-r-- 1 ben ben 6125321940 Sep  2 00:32 nigra-PM664.g.vcf.gz
-rw-rw-r-- 1 ben ben 6455292312 Sep  2 00:34 tonkeana-PM592.g.vcf.gz
```
Index the vcf files:
```
tabix -p vcf nemestrina-PM664.g.vcf.gz
```
Export the chrX:

```
~/tabix-0.2.6/tabix -h nemestrina-PM664.g.vcf.gz chrX > nemHiSeqchrX.vcf 
~/tabix-0.2.6/tabix -h nigra-PM664.g.vcf.gz chrX > nigraHiSeqchrX.vcf 
~/tabix-0.2.6/tabix -h tonkeana-PM592.g.vcf.gz chrX > tonkHiSeqchrX.vcf 
```

Convert from gvcf to vcf
```
/work/ben/2015_SulaRADtag/gvcftools-0.16/bin/break_blocks --ref /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --region-file /work/ben/2015_SulaRADtag/vcf-constitutional/target_interval_list_allchrs.bed < nemHiSeqchrX.vcf > nemHiSeqchrX.vcf.noblock.vcf

/work/ben/2015_SulaRADtag/gvcftools-0.16/bin/break_blocks --ref /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --region-file /work/ben/2015_SulaRADtag/vcf-constitutional/target_interval_list_allchrs.bed < tonkHiSeqchrX.vcf > tonkHiSeqchrX.vcf.noblock.vcf

/work/ben/2015_SulaRADtag/gvcftools-0.16/bin/break_blocks --ref /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --region-file /work/ben/2015_SulaRADtag/vcf-constitutional/target_interval_list_allchrs.bed < nigraHiSeqchrX.vcf > nigraHiSeqchrX.vcf.noblock.vcf
```
merge the files

```
bgzip nemHiSeqchrX.vcf.noblock.vcf
tabix -p vcf nemHiSeqchrX.vcf.noblock.vcf.gz
bgzip nigraHiSeqchrX.vcf.noblock.vcf
tabix -p vcf nigraHiSeqchrX.vcf.noblock.vcf.gz
bgzip tonkHiSeqchrX.vcf.noblock.vcf
tabix -p vcf tonkHiSeqchrX.vcf.noblock.vcf.gz

export PERL5LIB=/work/ben/vcftools/src/perl

/work/ben/vcftools/bin/vcf-merge nemHiSeqchrX.vcf.noblock.vcf.gz tonkHiSeqchrX.vcf.noblock.vcf.gz nigraHiSeqchrX.vcf.noblock.vcf.gz | bgzip -c > nem_tonk_nigra_alldiploid_chrX.vcf.gz
```

All of this appears to be bad because positions with no data are called as reference sites for some stupid reason.  So I am going to try again using GATK.
