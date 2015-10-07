# GATK and Base Recalibration

The bwa-assisted mapping with Stampy produced sam files that had ~10% more mapped reads for all individuals. Great!  Now the next step is to use GATK to do the base recalibration and genotype calling.  I have copied below a series of scripts I used for this purpose (on info).

## Add Read groups
The bam files recovered from Stampy lacked readgroup header information; this is needed by GATK for base recalibration so I ran this script:

``` perl
#!/usr/bin/perl
# This command will add read groups to stampy bam files
# I read on the GATK forum that each sample from multiplexed samples should have 
# a different readgroup ID, so that is what I am going to do
# https://www.broadinstitute.org/gatk/guide/article?id=1317

my $status;
my @files;
my $commandline;
@files = glob("*_stampy_sorted.bam");

for (0..$#files){
    $files[$_] =~ s/\.bam$//;
}

foreach(@files){
    print $_,"\n";
    $commandline=();
    $commandline = "java -Xmx2g -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups INPUT=".$_.".bam OUTPUT=".$_."_rg.bam RGLB=".$_." RGPL=Illumina RGPU=".$_." RGSM=".$_;
    $status = system($commandline);
}
```
## Realignment
Using the bam files from bwa, samtools, and stampy, I then realigned indels using GATK in two steps.

* Step 1.

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("*_stampy_sorted_rg.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -o stampy_forIndelRealigner_ym.intervals";


$status = system($commandline);


```

* Step 2

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("*_stampy_sorted_rg.bam");

for (0..$#files){
    $files[$_] =~ s/\_stampy_sorted_rg.bam$//;
}


my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T IndelRealigner ";

foreach(@files){
    $commandline = $commandline." -I ".$_."_stampy_sorted_rg.bam ";
}

$commandline = $commandline."-R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta --targetIntervals stampy_forIndelRealigner_ym.intervals --nWayOut .stampy_realigned.bam";

$status = system($commandline);

$status = system ("rename _stampy_sorted_rg.stampy_realigned.bam _stampy_realigned.bam *_stampy_sorted_rg.stampy_realigned.bam");
$status= system ("rename _stampy_sorted_rg.stampy_realigned.bai _stampy_realigned.bai *_stampy_sorted_rg.stampy_realigned.bai");

```

The last steps were needed to rename the silly long names GATK gave the realigned files.

## Genotyping

Then I genotyped the files as follows:

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("*_stampy_realigned.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline." -out_mode EMIT_VARIANTS_ONLY -o ./nonrecal_stampy_varonly_unifiedgenotyper_ym.vcf";


$status = system($commandline);
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

my $commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta ";

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

my $commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T PrintReads -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";

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

