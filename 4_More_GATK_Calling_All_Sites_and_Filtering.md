# Calling invariant and variant sites with GATK and then filtering

OK so we have finished the base calling.  Now we can export a vcf file with the full data, including homozygous sites and SNPs.

I did this with this script:

```perl
# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files that makes a vcf
# file with homozygous and heterozygous high quality calls, including INDELs.  

my $status;
my $file = "recal_stampy_round2_all.bam";



my $commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";

$commandline = $commandline." -I ".$file;

$commandline = $commandline." -out_mode EMIT_ALL_CONFIDENT_SITES --genotype_likelihoods_model BOTH -o recal_stampy_allsites_round2_all.vcf";

$status = system($commandline);

```

