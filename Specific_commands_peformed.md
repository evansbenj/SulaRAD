# This is a specific list of commands and scripts I have performed.

## BWA, Samtools, and GATK
On info, I demultiplexed the data with Stacks.  I then combined the redundant reads within the 95 samples using the cat command. I then added the forward read from the *M. tonkeana* data.

Then I moved these fastq files into a directory called `/home/ben/2015_SulaRADtag/good_merged_samples/fastq`

Then I executed a bash script called `others_alignmnet` which looked like this:

``` bash
#!/bin/bash                                                                                                                  

path_to_bwa="/usr/local/bin"
path_to_samtools="/usr/local/bin"
path_to_data="."
path_to_chromosome="/home/ben/2015_BIO720/rhesus_genome"
chromosome="macaque_masked_chromosomes_ym.fasta"

individuals="brunescens_PF707
hecki_PF643
hecki_PF644
hecki_PF648
hecki_PF651
hecki_PM639
hecki_PM645
maura_PF615
maura_PF713
maura_PM613
maura_PM614
maura_PM616
maura_PM618
nem_Gumgum
nem_Kedurang
nem_Malay
nem_Ngasang
nem_pagensis
nem_PM664
nem_PM665
nem_Sukai_male
nigra_PF1001
nigra_PF660
nigra_PM1000
nigra_PM1003
nigrescens_PF654
ochreata_PF625
ochreata_PM571
ochreata_PM596
togeanus_PF549
togeanus_PM545"

for each_individual in $individuals
do

echo ${each_individual}
    $path_to_bwa/bwa aln $path_to_chromosome/$chromosome $path_to_data/${each_individual}.fq > $path_to_data/${each_individual}.sai
    $path_to_bwa/bwa samse -r "@RG\tID:FLOWCELL1.LANE5\tSM:${each_individual}.fq\tPL:illumina" $path_to_chromosome/$chromosome $path_to_data/${each_individual}.sai $path_to_data/${each_individual}.fq > $
path_to_data/${each_individual}.sam
    $path_to_samtools/samtools view -bt $path_to_chromosome/$chromosome -o $path_to_data/${each_individual}.bam $path_to_data/${each_individual}.sam
    $path_to_samtools/samtools sort $path_to_data/${each_individual}.bam $path_to_data/${each_individual}_sorted
    $path_to_samtools/samtools index $path_to_data/${each_individual}_sorted.bam
    rm -f $path_to_data/${each_individual}.bam $path_to_data/${each_individual}.sam $path_to_data/${each_individual}.sai
done
```

This used `bwa` to align the single end reads, then it made sorted bam files and a bam index (bai) file.  I used GATK to realign these bam files using two perl script called `1_Executes_GATK_commands_RealignerTarget.pl` and `2_Executes_GATK_commands_IndelRealigner.pl`.

Here is the `1_Executes_GATK_commands_RealignerTarget.pl`:
```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("*_sorted.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -o forIndelRealigner_ym.intervals";


$status = system($commandline);
```

Here is `2_Executes_GATK_commands_IndelRealigner.pl`:
```perl 
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("*_sorted.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T IndelRealigner ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta --targetIntervals forIndelRealigner_ym.intervals --nWayOut .realigned.bam";


$status = system($commandline);
```


## Stampy
THen I used stampy to map some more reads using these bam files as a starting point.  This was accomplished using this a bash script called `stampy_alignment.sh`

``` bash
#!/bin/bash                                                                                                                  

path_to_stampy="/usr/local/stampy"
path_to_data="."
path_to_bam="../"
path_to_genome="/home/ben/2015_BIO720/rhesus_genome"
genome="rhesus_masked"

individuals="brunescens_PF707
hecki_PF643
hecki_PF644
hecki_PF648
hecki_PF651
hecki_PM639
hecki_PM645
maura_PF615
maura_PF713
maura_PM613
maura_PM614
maura_PM616
maura_PM618
nem_Gumgum
nem_Kedurang
nem_Malay
nem_Ngasang
nem_pagensis
nem_PM664
nem_PM665
nem_Sukai_male
nigra_PF1001
nigra_PF660
nigra_PM1000
nigra_PM1003
nigrescens_PF654
ochreata_PF625
ochreata_PM571
ochreata_PM596
togeanus_PF549
togeanus_PM545"

for each_individual in $individuals
do

echo ${each_individual}
#    $path_to_stampy/stampy.py -g $path_to_genome/$genome -h $path_to_genome/$genome --substitutionrate=0.013 -t8 -M $path_to_data/${each_individual}.fq
     $path_to_stampy/stampy.py -g $path_to_genome/$genome -h $path_to_genome/$genome --substitutionrate=0.013 -t8 --bamkeepgoodreads -M $path_to_bam/${each_individual}_sorted.realigned.bam -f sam -o $pat
h_to_bam/${each_individual}_stampy.sam
done
```
Using the samtools flagstat command (after converting the sam to bam files, I had originally mapped about 60% of the reads and after stampy I managed to map about 70% of the reads.  Of note is that although the original input bam files had been realigned with GATK to maximally align indels, the new reads could still be misaligned, so this will have to be done again. I am also not sure whether the new bam files made from the stampy sam files will retain the header information I added earlier, so I may have to do this again.

Here is a script called `stampy_sam_to_bam.sh` that makes sorted bam files out of the stampy bam files and also makes a bam index for these files. 

```bash
#!/bin/bash                                                                                                                                                                                                

path_to_bwa="/usr/local/bin"
path_to_samtools="/usr/local/bin"
path_to_data="."
path_to_chromosome="/home/ben/2015_BIO720/rhesus_genome"
chromosome="macaque_masked_chromosomes_ym.fasta"

individuals="brunescens_PF707                                                                                                                                                                              
hecki_PF643                                                                                                                                                                                                
hecki_PF644                                                                                                                                                                                                
hecki_PF648                                                                                                                                                                                                
hecki_PF651                                                                                                                                                                                                
hecki_PM639                                                                                                                                                                                                
hecki_PM645                                                                                                                                                                                                
maura_PF615                                                                                                                                                                                                
maura_PF713                                                                                                                                                                                                
maura_PM613                                                                                                                                                                                                
maura_PM614                                                                                                                                                                                                
maura_PM616                                                                                                                                                                                                
maura_PM618                                                                                                                                                                                                
nem_Gumgum                                                                                                                                                                                                 
nem_Kedurang                                                                                                                                                                                               
nem_Malay                                                                                                                                                                                                  
nem_Ngasang                                                                                                                                                                                                
nem_pagensis                                                                                                                                                                                               
nem_PM664                                                                                                                                                                                                  
nem_PM665                                                                                                                                                                                                  
nem_Sukai_male                                                                                                                                                                                             
nigra_PF1001                                                                                                                                                                                               
nigra_PF660                                                                                                                                                                                                
nigra_PM1000                                                                                                                                                                                               
nigra_PM1003                                                                                                                                                                                               
nigrescens_PF654                                                                                                                                                                                           
ochreata_PF625                                                                                                                                                                                             
ochreata_PM571                                                                                                                                                                                             
ochreata_PM596                                                                                                                                                                                             
togeanus_PF549                                                                                                                                                                                             
togeanus_PM545"

for each_individual in $individuals
do

echo ${each_individual}
    $path_to_samtools/samtools view -bt $path_to_chromosome/$chromosome -o $path_to_data/${each_individual}_stampy.bam $path_to_data/${each_individual}_stampy.sam
    $path_to_samtools/samtools sort $path_to_data/${each_individual}_stampy.bam $path_to_data/${each_individual}_stampy_sorted
    $path_to_samtools/samtools index $path_to_data/${each_individual}_stampy_sorted.bam
    rm -f $path_to_data/${each_individual}_stampy.bam $path_to_data/${each_individual}_stampy.sam
done
```
