# Calling invariant and variant sites with GATK and then filtering

OK so we have finished the base calling.  Now we can export a vcf file with the full data, including homozygous sites, SNPs, and indels.  To get indels as well (which will be used for filtering purposes in the next step), I need to  "EMIT_ALL_CONFIDENT_SITES".  I named the file this: "recal_stampy_allsites_round2_all_confident_sites.vcf".

I did this with this script:

```perl
#!/usr/bin/perl
# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files that makes a vcf
# file with homozygous and heterozygous high quality calls, including INDELs.  

my $status;
my $file = "recal_stampy_round2_all.bam";



my $commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";

$commandline = $commandline." -I ".$file;

$commandline = $commandline." -out_mode EMIT_ALL_CONFIDENT_SITES --genotype_likelihoods_model BOTH -o recal_stampy_allsites_round2_all_confident_sites.vcf";

$status = system($commandline);

```

Now I need to make a vcf file with only indels for use in filtering.  Here is a perl script that will do that:

```perl
#!/usr/bin/perl
# This script will read in a vcf file names and 
# make and execute a GATK commandline that outputs only INDELs
# in a new vcf file.  

my $status;
my $file = "recal_stampy_allsites_round2_all_confident_sites.vcf";

my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o round2_indels_only.vcf --variant ".$file." -selectType INDEL";

$status = system($commandline);

```

OK, now mark indels, bits near indels, and other low quality stuff.  This command should mark indels plus and minus 5 bp, aDNA with genotype qualities lower than 40, sites with 1/5th or more of the reads with good map qualities to other parts of the genome, and any site with a coverage less than 5.

``` perl
#!/usr/bin/perl
# This script will read in a vcf file and                                                                                       # make and execute a GATK commandline that marks INDELs and other stuff                                                     
# in this vcf file.                                                                                                                                                                     
my $status;
my $file = "recal_stampy_allsites_round2_all_confident_sites.vcf";

my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";
$commandline = $commandline." -o round2_marked.vcf --variant ".$file.";
$commandline = $commandline." --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpressio\
n \"DP < 5\" --filterName \"TooLowCoverage\" --filterExpression \"DP > 50\" --filterName \"TooHighCoverage\" --filterExpression \"QUAL < 40.0 && CHROM != \'chrY\' && CHROM != \'chrX\'\" --filterName \"LowQual_aDNA\" --filterExpression \"QUAL < 20.0 && CHROM == \'chrX\'\" --filterName \"LowQual_xDNA\" --filterExpression \"QUAL < 20.0 && CHROM == \'chrY\'\" --filterName \"LowQual_yDNA\" --mask round2_indels_only.vcf --maskName INDEL --maskExtension 5";

$status = system($commandline);

```

Now I am going to add some cool filtering for sex chromosomes based on the sex of the individual.  Eventually I want to tweak this so it also filters by coverage depending on the sex of each individual

``` perl
#!/usr/bin/perl
# This script will read in a vcf file and                                                                                       # make and execute a GATK commandline that marks INDELs and other stuff                                                         # in this vcf file.                                                                                                                                                                                                                     
my $status;
my $file = "recal_stampy_allsites_round2_all_confident_sites.vcf";
my $outfile1="round2_marked.vcf";

my @females=("brunescens_PF707_stampy_sorted","hecki_PF643_stampy_sorted","hecki_PF644_stampy_sorted","hecki_PF648_stampy_sorted","hecki_PF651_stampy_sorted","maura_PF615_stampy_sorted","maura_PF713_stampy_sorted","nem_Gumgum_stampy_sorted","nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted","nem_pagensis_stampy_sorted","nigra_PF1001_stampy_sorted","nigra_PF660_stampy_sorted","nigrescens_PF654_stampy_sorted","ochreata_PF625_stampy_sorted","togeanus_PF549_stampy_sorted","tonk_PF515_stampy_sorted");

my @males=("hecki_PM639_stampy_sorted","hecki_PM645_stampy_sorted","maura_PM613_stampy_sorted","maura_PM614_stampy_sorted","maura_PM616_stampy_sorted","maura_PM618_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted","nigra_PM1000_stampy_sorted","nigra_PM1003_stampy_sorted","ochreata_PM571_stampy_sorted","ochreata_PM596_stampy_sorted","togeanus_PM545_stampy_sorted","tonk_PM561_stampy_sorted","tonk_PM565_stampy_sorted","tonk_PM566_stampy_sorted","tonk_PM567_stampy_sorted","tonk_PM582_stampy_sorted","tonk_PM584_stampy_sorted","tonk_PM592_stampy_sorted","tonk_PM602_stampy_sorted");

my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";
$commandline = $commandline." -o ".$outfile1." --variant ".$file;

# aDNA filtering and indels                                                                                                                                                                                                             
# This includes a multiple mapping filter based on having at least 4 reads mapping to another position and where together these reads comprised at least 10% of the total reads (HARD_TO_VALIDATE)                                      
# where the depth of reads summed over all called genotypes divided by the number of alleles in all called genotypes is less than 5 (TooLowCoverage)                                                                                    
# or greater than 50 (TooHighCoverage)                                                                                                                                                                                                  
# or if a position is within 5 bp of an indel (INDEL)                                                                                                                                                                                   
$commandline = $commandline." --filterExpression \"CHROM != \'chrY\' && CHROM != \'chrX\' && MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"CHROM != \'chrY\' && CHROM != \'chrX\' && ((DP / AN) < 5)\" --filterName \"TooLowCoverage\" --filterExpression \"CHROM != \'chrY\' && CHROM != \'chrX\' && ((DP/AN) > 100)\" --filterName \"TooHighCoverage\" --filterExpression \"CHROM != \'chrY\' && CHROM != \'chrX\' && QUAL < 30.0\" --filterName \"LowQual_aDNA\" --mask round2_indels_only.vcf --maskName INDEL --maskExtension 5";

# filter the PAR from everyone, which, based on Hughes et al. 2012 Suppl Fig. 7b is everything less than position 77221 in the chrY from that paper                                                                                     
$commandline=$commandline." --filterExpression \"CHROM == \'chrY\' && POS < 77221\" --filterName \"Y_PAR\" ";

# filter all sites in which a female genotype is calledon chrY                                                                                                                                                                          
foreach (@females){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrY\'\" --filterName \"female_Y_filter_".$_."\"";
}

# filter all sites in which a male heteroz genotype is called on chrX                                                                                                                                                                   
foreach (@males){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrX\' && vc.getGenotype(\'".$_."\').isHet()\" --filterName "male_Xhet_filter_".$_."\"";
}

# filter all sites in which a male heteroz genotype is called on chrY                                                                                                                                                                   
foreach (@males){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrY\' && vc.getGenotype(\'".$_."\').isHet()\" --filterName "male_Yhet_filter\"";
}

$status = system($commandline);

```
