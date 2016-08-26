# GATK and Base Recalibration

The bwa-assisted mapping with Stampy produced sam files that had ~10% more mapped reads for all individuals. Great!  Now the next step is to use GATK to do the base recalibration and genotype calling.  I have copied below a series of scripts I used for this purpose (on info).

# Genotyping with Haplotype caller

Some issues with screen on info15-20 were resolved like this: `export TERM=xterm-color`.

I am going to write the scripts from the new analyses above the old ones. Here is the Haplotype caller commands to generate confident variants for BSQR (3_Executes_GATK_commands_Haplotypecaller.pl).

```
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("fastq/*_trimmed_sorted.realigned.bam");

my $commandline = "java -Xmx1G -jar  /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar
 -T HaplotypeCaller -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym
.fasta ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline." -L fastq/target_interval_list_autosomes.list -out_mode EMI
T_VARIANTS_ONLY -o fastq/varonly_haplotypecaller_for_BSQR.vcf";


$status = system($commandline);
```
## Filtering the set of known sites
One concern I have is that, no matter how much I filter, repetitive sites may still have mismapped reads. To me it therefore makes sense to filter the set of known variants and exclude the repetitive regions.  I also wish to exclude the sex chromosomes and mtDNA, the latter of which almost certainly has mismapped numts. So I filtered this out using a bed file I made by combining all of the aDNA repeats from repeatmasker plus the entire X, Y, and mt genomes.  This file is called `aDNA_bed_allX_allY_allM.bed` and here is the command I used to filter the file:

```
vcftools --gzvcf varonly_haplotypecaller_for_BSQR.vcf.gz --recode --out varonly_haplotypecaller_norepeats_noX_noY_noM_for_BSQR.vcf.gz --exclude-bed ../repeat_masker_chromOut/aDNA_bed_allX_allY_allM.bed
```


## Base recalibration

Then I did the base recalibration as follows:

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("*_stampy_realigned.bam");

my $commandline = "java -Xmx3G -jar  /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline." -knownSites ./nonrecal_stampy_varonly_unifiedgenotyper_ym.vcf -o recal_data.table";


$status = system($commandline);
```

## Print Reads
Then I output a new bam file with recalibrated base quality scores as follows

``` perl

#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("*_stampy_realigned.bam");

my $commandline = "java -Xmx3G -jar  /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T PrintReads -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-BQSR recal_data.table -o recal_stampy_round1_all.bam";


$status = system($commandline);
```
## Genotype the new bam file

Then I genotyped the new and improved bam file as follows

``` perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my $file = "recal_stampy_round1_all.bam";
#my @files;
   
#@files = glob("*_sorted.realigned.bam");


my $commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";

$commandline = $commandline." -I ".$file;

#foreach(@files){
#    $commandline = $commandline." -I ".$_." ";
#}

$commandline = $commandline." -o recal_stampy_varonly_round1_all.vcf";


$status = system($commandline);
```

## Compare the new variants to the Variants from the non-recalibrated genotyping

Then I used vcftools to compare the new variants to the original genotype calls. I found that the new variants were almost entirely a subset of the original calls (>99.7%), but that the original calls had a substantial fraction that were not in the new variants (~22%). These presumably are errors.  I then repeated the base recalibration using the new and improved set of known SNPs, made a new bam file, repeated the genotype calls, and output a new "round2" set of snps.  This can be compared to "round1" snps and hopefully it will be similar.

This can be accomplished by first gzipping and tabixing the vcf files as follows:

`/usr/local/tabix/bgzip recal_stampy_varonly_round2_all.vcf`

and then

`/usr/local/tabix/tabix -p vcf recal_stampy_varonly_round2_all.vcf.gz`

and then comparing two bgzipped vcf files like this:

`/usr/local/vcftools/src/perl/vcf-compare xxx.vcf.gz yyy.vcf.gz > compare.out`

Here is an example of the output when I compare round 2 to the original nonrecalibrated variants:

```
[ben@info115 good_merged_samples]$ ~/vcftools/src/perl/vcf-compare nonrecal_stampy_varonly_unifiedgenotyper_ym.vcf.gz recal_stampy_varonly_round2_all.vcf.gz > stampy_nonrecal_to_round2_vcf_compare.out
[ben@info115 good_merged_samples]$ more stampy_nonrecal_to_round2_vcf_compare.out
# This file was generated by vcf-compare.
# The command line was: vcf-compare(r953) nonrecal_stampy_varonly_unifiedgenotyper_ym.vcf.gz recal_stampy_varonly_round2_all.vcf.gz
#
#VN 'Venn-Diagram Numbers'. Use `grep ^VN | cut -f 2-` to extract this part.
#VN The columns are: 
#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN	494	recal_stampy_varonly_round2_all.vcf.gz (0.1%)
VN	103983	nonrecal_stampy_varonly_unifiedgenotyper_ym.vcf.gz (23.3%)
VN	342551	nonrecal_stampy_varonly_unifiedgenotyper_ym.vcf.gz (76.7%)	recal_stampy_varonly_round2_all.vcf.gz (99.9%)
#SN Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN	Number of REF matches:	342551
SN	Number of ALT matches:	341040
SN	Number of REF mismatches:	0
SN	Number of ALT mismatches:	1511
SN	Number of samples in GT comparison:	0
```

and here is when I get when I compare round 1 to round 2:

```
[ben@info115 good_merged_samples]$ ~/vcftools/src/perl/vcf-compare recal_stampy_varonly_round1_all.vcf.gz recal_stampy_varonly_round2_all.vcf.gz > stampy_round1_to_round2_vcf_compare.out
[ben@info115 good_merged_samples]$ more stampy_round1_to_round2_vcf_compare.out
# This file was generated by vcf-compare.
# The command line was: vcf-compare(r953) recal_stampy_varonly_round1_all.vcf.gz recal_stampy_varonly_round2_all.vcf.gz
#
#VN 'Venn-Diagram Numbers'. Use `grep ^VN | cut -f 2-` to extract this part.
#VN The columns are: 
#VN        1  .. number of sites unique to this particular combination of files
#VN        2- .. combination of files and space-separated number, a fraction of sites in the file
VN	1223	recal_stampy_varonly_round2_all.vcf.gz (0.4%)
VN	5475	recal_stampy_varonly_round1_all.vcf.gz (1.6%)
VN	341822	recal_stampy_varonly_round1_all.vcf.gz (98.4%)	recal_stampy_varonly_round2_all.vcf.gz (99.6%)
#SN Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN	Number of REF matches:	341822
SN	Number of ALT matches:	341660
SN	Number of REF mismatches:	0
SN	Number of ALT mismatches:	162
SN	Number of samples in GT comparison:	0
```

As you can see, there were 23.3% fewer snps called after base recalibration. This seems like a lot. However, the next round of recalibration produces almost dientical results. So let's go with that! 
