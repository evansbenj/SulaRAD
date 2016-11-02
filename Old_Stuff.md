








# STUFF BELOW WAS NOT USED!!!

and 3_Executes_GATK_commands_Haplotypecaller_for_BSQR.pl:

```
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $status;
my @files;
   
@files = glob("fastq/*_trimmed_sorted.realigned.bam");

foreach(@files){
	my $commandline = "java -Xmx1G -jar  /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -I ".$_." -L fastq/target_interval_list_autosomes.list -out_mode EMIT_VARIANTS_ONLY --emitRefConfidence GVCF -o ".$_."varonly_haplotypecaller_for_BSQR.g.vcf";
	print $commandline,"\n";
	$status = system($commandline);
}
```



### Sex chromosomes (3B_extracts_Xs.pl)
``` perl

#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make new bam files with only chrX

my $status;
my @files;
my $commandline;
   
#@files = glob("fastq/*_trimmed_sorted.realigned.bam");

#foreach(@files){
#    $commandline = "samtools view -b ".$_." chrX > ".$_."_chrX.bam";
#    print $commandline,"\n";
#    $status = system($commandline);
#    $commandline = "samtools index ".$_."_chrX.bam";
#    print $commandline,"\n";
#    $status = system($commandline);
#   # $commandline = "java -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups I=".$_."_chrX.bam O=".$_."_chrX_rg.bam RGID=FLOWCELL1.LANE5 RGLB=".$_.".fq RGPL=illumina RGPU=".$_.".fq RGSM=".$_.".fq";
#   # $status = system($commandline);
#}

my @males=("fastq/hecki_PM639_trimmed_sorted.realigned.bam_chrX.bam","fastq/hecki_PM645_trimmed_sorted.realigned.bam_chrX.bam","fastq/maura_PM613_trimmed_sorted.realigned.bam_chrX.bam","fastq/maura_PM614_trimmed_sorted.realigned.bam_chrX.bam","fastq/maura_PM616_trimmed_sorted.realigned.bam_chrX.bam","fastq/maura_PM618_trimmed_sorted.realigned.bam_chrX.bam","fastq/nem_PM664_trimmed_sorted.realigned.bam_chrX.bam","fastq/nem_PM665_trimmed_sorted.realigned.bam_chrX.bam","fastq/nem_Sukai_male_trimmed_sorted.realigned.bam_chrX.bam","fastq/nigra_PM1000_trimmed_sorted.realigned.bam_chrX.bam","fastq/nigra_PM1003_trimmed_sorted.realigned.bam_chrX.bam","fastq/ochreata_PM571_trimmed_sorted.realigned.bam_chrX.bam","fastq/ochreata_PM596_trimmed_sorted.realigned.bam_chrX.bam","fastq/togeanus_PM545_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PM561_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PM565_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PM566_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PM567_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PM582_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PM584_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PM592_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PM602_trimmed_sorted.realigned.bam_chrX.bam");
my @females=("fastq/brunescens_PF707_trimmed_sorted.realigned.bam_chrX.bam","fastq/hecki_PF643_trimmed_sorted.realigned.bam_chrX.bam","fastq/hecki_PF644_trimmed_sorted.realigned.bam_chrX.bam","fastq/hecki_PF648_trimmed_sorted.realigned.bam_chrX.bam","fastq/hecki_PF651_trimmed_sorted.realigned.bam_chrX.bam","fastq/maura_PF615_trimmed_sorted.realigned.bam_chrX.bam","fastq/maura_PF713_trimmed_sorted.realigned.bam_chrX.bam","fastq/nem_Gumgum_trimmed_sorted.realigned.bam_chrX.bam","fastq/nem_Kedurang_trimmed_sorted.realigned.bam_chrX.bam","fastq/nem_Malay_trimmed_sorted.realigned.bam_chrX.bam","fastq/nem_Ngasang_trimmed_sorted.realigned.bam_chrX.bam","fastq/nem_pagensis_trimmed_sorted.realigned.bam_chrX.bam","fastq/nigra_PF1001_trimmed_sorted.realigned.bam_chrX.bam","fastq/nigra_PF660_trimmed_sorted.realigned.bam_chrX.bam","fastq/nigrescens_PF654_trimmed_sorted.realigned.bam_chrX.bam","fastq/ochreata_PF625_trimmed_sorted.realigned.bam_chrX.bam","fastq/togeanus_PF549_trimmed_sorted.realigned.bam_chrX.bam","fastq/tonk_PF515_trimmed_sorted.realigned.bam_chrX.bam");

foreach(@females){
    my $commandline = "java -Xmx1G -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -I ".$_." -out_mode EMIT_ALL_CONFIDENT_SITES -L fastq/target_interval_list_chrX.list --emitRefConfidence GVCF -o ".$_."_xDNA_no_BSQR_all_confident_sites.g.vcf";
    print $commandline,"\n";
    $status = system($commandline);
}

foreach(@males){
    my $commandline = "java -Xmx1G -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -I ".$_." -ploidy 1 -out_mode EMIT_ALL_CONFIDENT_SITES -L fastq/target_interval_list_chrX.list --emitRefConfidence GVCF -o ".$_."_xDNA_no_BSQR_all_confident_sites.g.vcf";
    print $commandline,"\n";    
    $status = system($commandline);
}

```

## Combine the ploidy vcf files (3C_Merge_chrX_vcfs_noBSQR.pl)

```perl
#!/usr/bin/perl                                                                                        
use warnings;
use strict;

# This script will merge chrX files 

my $status;
my @files;
my $commandline;

@files = ("fastq/brunescens_PF707_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/hecki_PF643_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/hecki_PF644_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/hecki_PF648_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/hecki_PF651_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/hecki_PM639_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/hecki_PM645_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/maura_PF615_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/maura_PF713_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/maura_PM613_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/maura_PM614_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/maura_PM616_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/maura_PM618_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nem_Gumgum_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nem_Kedurang_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nem_Malay_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nem_Ngasang_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nem_PM664_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nem_PM665_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nem_Sukai_male_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nem_pagensis_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nigra_PF1001_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nigra_PF660_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nigra_PM1000_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nigra_PM1003_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/nigrescens_PF654_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/ochreata_PF625_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/ochreata_PM571_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/ochreata_PM596_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/togeanus_PF549_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/togeanus_PM545_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PF515_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PM561_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PM565_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PM566_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PM567_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PM582_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PM584_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PM592_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf","fastq/tonk_PM602_trimmed_sorted.realigned.bam_chrX.bam_xDNA_no_BSQR_all_confident_sites.g.vcf");

$commandline = "java -Xmx5g -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta ";
foreach(@files){
    $commandline = $commandline."--variant ".$_." ";
}   
$commandline = $commandline."--includeNonVariantSites -o RADseq_chrX_haplotypecaller_combined.vcf";       

print $commandline,"\n";

$status = system($commandline);

```
#  Filter combined genotypes
Filter (3D_Filter_chrX_vcfs_noBSQR.pl)

```perl
#!/usr/bin/perl                                                                                        
use warnings;
use strict;

# This script will filter chrX files 

my $status;
my $commandline;

$commandline = "java -Xmx5g -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -o /home/ben/2015_SulaRADtag/good_merged_samples/RADseq_chrX_haplotypecaller_combined_marked.vcf --variant /home/ben/2015_SulaRADtag/good_merged_samples/RADseq_chrX_haplotypecaller_combined.vcf --filterExpression \"QUAL < 30.0\" --filterName \"VeryLowQual\" -G_filter \"DP < 4\" -G_filterName \"LowCoverage\" --setFilteredGtToNocall";
print $commandline,"\n";
$status = system($commandline);


$commandline = "java -Xmx5g -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -o /home/ben/2015_SulaRADtag/good_merged_samples/RADseq_chrX_haplotypecaller_combined_filtered.vcf --variant /home/ben/2015_SulaRADtag/good_merged_samples/RADseq_chrX_haplotypecaller_combined_marked.vcf --excludeFiltered --maxNOCALLfraction 0.99";
print $commandline,"\n";
$status = system($commandline);

```
# Combine the vcf files from each chromosome from GenotypGVCF
I am going to use GATK instead of bcftools because there may be weird incompatibility issues with the gvcf files...

```
java -Xmx1G -jar  /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T CombineGVCFs -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta --variant fastq/GenotypeVCFs_noBSQR.vcf.gz --variant fastq/GenotypeVCFs_noBSQR_6_to_10.vcf.gz --variant fastq/GenotypeVCFs_noBSQR_11_to_16.vcf.gz --variant fastq/GenotypeVCFs_noBSQR_17_to_X.vcf.gz -o fastq/GenotypeVCFs_noBSQR_concat.vcf.gz
```


# Combine vcf files from each chromosome
java -cp /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta --variant fastq/GenotypeVCFs_noBSQR.vcf.gz --variant fastq/GenotypeVCFs_noBSQR_6_to_10.vcf.gz --variant fastq/GenotypeVCFs_noBSQR_11_to_16.vcf.gz --variant fastq/GenotypeVCFs_noBSQR_17_to_X.vcf.gz --outputFile fastq/GenotypeVCFs_noBSQR_concat.vcf.gz -assumeSorted

# BELOW WAS NOT USED


## Stampy (I am no longer using Stampy, so this section was skipped)
Then I used stampy to map some more reads using these bam files as a starting point.  This was accomplished using this a bash script called `stampy_alignment.sh`.  But before this is executed, we need to load python 2.7 like this (on info):

`scl enable python27 bash`



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
togeanus_PM545
tonk_PF515
tonk_PM561
tonk_PM565
tonk_PM566
tonk_PM567
tonk_PM582
tonk_PM584
tonk_PM592
tonk_PM602"

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
togeanus_PM545
tonk_PF515
tonk_PM561
tonk_PM565
tonk_PM566
tonk_PM567
tonk_PM582
tonk_PM584
tonk_PM592
tonk_PM602"

for each_individual in $individuals
do

echo ${each_individual}
    $path_to_samtools/samtools view -bt $path_to_chromosome/$chromosome -o $path_to_data/${each_individual}_stampy.bam $path_to_data/${each_individual}_stampy.sam
    $path_to_samtools/samtools sort $path_to_data/${each_individual}_stampy.bam $path_to_data/${each_individual}_stampy_sorted
    $path_to_samtools/samtools index $path_to_data/${each_individual}_stampy_sorted.bam
    rm -f $path_to_data/${each_individual}_stampy.bam $path_to_data/${each_individual}_stampy.sam
done
```
