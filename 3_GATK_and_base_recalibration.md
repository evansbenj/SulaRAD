# GATK and Base Recalibration

The bwa-assisted mapping with Stampy produced sam files that had ~10% more mapped reads for all individuals. Great!  Now the next step is to use GATK to do the base recalibration and genotype calling.  I have copied below a series of scripts I used for this purpose (on info).

## Add Read groups
The bam files recovered from Stampy lacked readgroup header information; this is needed by GATK for base recalibration so I ran this script:

``` perl
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
Then I output a vcf file as follows

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
