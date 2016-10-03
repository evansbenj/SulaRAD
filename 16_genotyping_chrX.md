# Genotyping chrX is difficult

Background: my preliminary results suggest that there are some issues with genotyping chrX in males.  I have tried doing this with Haplotype caller using the default setting for ploidy (=2) and setting ploidy=1.  The former has problems because many heterozygous calls are made.  The latter has problems because the genotypes have a bias towards reference sites.

As a solution, Janet Kelso from the MPI suggested I genotype male chrX sites using the highest frequenty SNP.  This sounds like a great idea to me!  I plan to do this for the female as well so that each chrX is treated the same. I will write a script that outputs a tab delimited file from a diploid vcf file with genotyping based on the AD annotation (AD is allele depth).  For sites with an equal frequency of the ref and alt SNP, I will select the genotype randomly (and also keep track of how frequently this happens.

# First output chrX
