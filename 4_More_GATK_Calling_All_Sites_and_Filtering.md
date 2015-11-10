# Calling invariant and variant sites with GATK and then filtering

OK so we have finished the base calling.  Now we can export a vcf file with the full data, including homozygous sites, SNPs, and indels.  To get indels as well (which will be used for filtering purposes in the next step), I need to  "EMIT_ALL_CONFIDENT_SITES".  I named the file this: "recal_stampy_allsites_round2_all_confident_sites.vcf".

I did this with this script (10_Executes_GATK_commands_UnifiedGenotyper_emitALL.pl):

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

Now I need to make a vcf file with only indels for use in filtering.  Here is a perl script that will do that (11_Executes_GATK_commands_SelectVariants_output_indels.pl):

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

Now make a vcf file with only the sites with bad sex chromosome genotypes using this script (12_Executes_GATK_commands_makes_sex_mask.pl):

``` perl
#!/usr/bin/perl

# This script will read in a vcf file and                                                                                                                                                
# make and execute a GATK commandline that outputs
# a file with only bad sex chromosome genotypes.
# This file will then be used as a mask
                                                                                                                                                                     
my $status;
my $file = "recal_stampy_allsites_round2_all_confident_sites.vcf";
my $outfile1="round2_BADSEX_marked.vcf";


my @females=("brunescens_PF707_stampy_sorted","hecki_PF643_stampy_sorted","hecki_PF644_stampy_sorted","hecki_PF648_stampy_sorted","hecki_PF651_stampy_sorted","maura_PF615_stampy_sorted","maura_PF713_stampy_sorted","nem_Gumgum_stampy_sorted","nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted","nem_pagensis_stampy_sorted","nigra_PF1001_stampy_sorted","nigra_PF660_stampy_sorted","nigrescens_PF654_stampy_sorted","ochreata_PF625_stampy_sorted","togeanus_PF549_stampy_sorted","tonk_PF515_stampy_sorted");

my @males=("hecki_PM639_stampy_sorted","hecki_PM645_stampy_sorted","maura_PM613_stampy_sorted","maura_PM614_stampy_sorted","maura_PM616_stampy_sorted","maura_PM618_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted","nigra_PM1000_stampy_sorted","nigra_PM1003_stampy_sorted","ochreata_PM571_stampy_sorted","ochreata_PM596_stampy_sorted","togeanus_PM545_stampy_sorted","tonk_PM561_stampy_sorted","tonk_PM565_stampy_sorted","tonk_PM566_stampy_sorted","tonk_PM567_stampy_sorted","tonk_PM582_stampy_sorted","tonk_PM584_stampy_sorted","tonk_PM592_stampy_sorted","tonk_PM602_stampy_sorted");

my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";
$commandline = $commandline." -o ".$outfile1." --variant ".$file;

# filter all sites in which a female genotype is calledon chrY
foreach (@females){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrY\'\" --filterName \"BAD_SEX\"";
}   

# filter all sites in which a male heteroz genotype is called on chrX
foreach (@males){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrX\' && vc.getGenotype(\'".$_."\').isHet()\" --filterName \"BAD_SEX\""; 
}   

# filter all sites in which a male heteroz genotype is called on chrY
foreach (@males){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrY\' && vc.getGenotype(\'".$_."\').isHet()\" --filterName \"BAD_SEX\"";
}

$status = system($commandline);

# Now output a vcf file that contains only the sites with a BAD_SEX flag
my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o round2_BAD_SEX_only.vcf --variant ".$file." -select \'vc.isFiltered()\' -invertSelect";

$status = system($commandline);

```

Now get rid of indels and bad_sex genotypes plus a buffer of 3 bp and 200 bp respectively (13_Executes_GATK_commands_VariantFiltration_doublemask.pl):

``` perl
#!/usr/bin/perl

# This script will read in a vcf file and 
# make and execute several GATK commandlines 
# that mark and remove INDELS and BAD_SEX genotypes
# plus a buffer (3 bp and 200 bp respectively)


my $status;
my $file = "recal_stampy_allsites_round2_all_confident_sites.vcf";
my $outfile1 = "temp1.vcf";
my $outfile2 = "temp2.vcf";
my $outfile3 = "round2_filtered.vcf";

# Mark sites that map poorly and that are in or within 3 bp of an indel
my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile1." --variant ".$file." --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\"";
$commandline = $commandline." --mask round2_indels_only.vcf --maskName INDEL --maskExtension 3";
$status = system($commandline);

# Mark sites that have inappropriate sex chromosome genotypes
$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile2." --variant ".$outfile1;
$commandline = $commandline." --mask round2_BAD_SEX_only.vcf --maskName BADSEX --maskExtension 200";
print $commandline,"\n";
$status = system($commandline);

# Print out the filtered file
$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile3." --variant ".$outfile2." -select 'vc.isNotFiltered()'";
$status = system($commandline);

# now clean up the intermediate vcf files
$commandline = "rm -f ".$outfile1." ".$outfile2;
$status = system($commandline);
```
