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
```
~/tabix-0.2.6/tabix -h nemestrina-PM664.g.vcf.gz chrX > nemHiSeqchrX.vcf 
~/tabix-0.2.6/tabix -h nigra-PM664.g.vcf.gz chrX > nigraHiSeqchrX.vcf 
~/tabix-0.2.6/tabix -h tonkeana-PM592.g.vcf.gz chrX > tonkHiSeqchrX.vcf 
