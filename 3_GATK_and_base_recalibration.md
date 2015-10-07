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
