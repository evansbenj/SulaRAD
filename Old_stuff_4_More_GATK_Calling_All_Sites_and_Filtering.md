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

Now make a vcf file with only the sites with bad sex chromosome genotypes using this script (12_Executes_GATK_commands_makes_sex_mask_VariantFiltration.pl):

``` perl
#!/usr/bin/perl

# This script will read in a vcf file and                                                                                                                                                
# make and execute a GATK commandline that outputs
# a file with only bad sex chromosome genotypes.
# This file will then be used as a mask
                                                                                                                                                                     
my $status;
my $file = "recal_stampy_allsites_round2_all_confident_sites.vcf";
my $outfile1="round2_BADSEX_marked.vcf";
my $outfile2="round2_BADSEX_only.vcf";

my @females=("brunescens_PF707_stampy_sorted","hecki_PF643_stampy_sorted","hecki_PF644_stampy_sorted","hecki_PF648_stampy_sorted","hecki_PF651_stampy_sorted","maura_PF615_stampy_sorted","maura_PF713_stampy_sorted","nem_Gumgum_stampy_sorted","nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted","nem_pagensis_stampy_sorted","nigra_PF1001_stampy_sorted","nigra_PF660_stampy_sorted","nigrescens_PF654_stampy_sorted","ochreata_PF625_stampy_sorted","togeanus_PF549_stampy_sorted","tonk_PF515_stampy_sorted");

my @males=("hecki_PM639_stampy_sorted","hecki_PM645_stampy_sorted","maura_PM613_stampy_sorted","maura_PM614_stampy_sorted","maura_PM616_stampy_sorted","maura_PM618_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted","nigra_PM1000_stampy_sorted","nigra_PM1003_stampy_sorted","ochreata_PM571_stampy_sorted","ochreata_PM596_stampy_sorted","togeanus_PM545_stampy_sorted","tonk_PM561_stampy_sorted","tonk_PM565_stampy_sorted","tonk_PM566_stampy_sorted","tonk_PM567_stampy_sorted","tonk_PM582_stampy_sorted","tonk_PM584_stampy_sorted","tonk_PM592_stampy_sorted","tonk_PM602_stampy_sorted");

my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";
$commandline = $commandline." -o ".$outfile1." --variant ".$file;

# filter all sites in which a female genotype is called on chrY
foreach (@females){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrY\' && vc.getGenotype(\'".$_."\').isHet()\" --filterName BAD_SEXfY";
    $commandline=$commandline." --filterExpression \"CHROM == \'chrY\' && vc.getGenotype(\'".$_."\').isHom()\" --filterName BAD_SEXfY";}   

# filter all sites in which a male heteroz genotype is called on chrX
foreach (@males){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrX\' && vc.getGenotype(\'".$_."\').isHet()\" --filterName BAD_SEXmX"; 
}   

# filter all sites in which a male heteroz genotype is called on chrY
foreach (@males){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrY\' && vc.getGenotype(\'".$_."\').isHet()\" --filterName BAD_SEXmY";
}
print $commandline,"\n";
$status = system($commandline);

# Now output a vcf file that contains only the sites with a BAD_SEX flag
my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile2." --variant ".$outfile1." -select \'vc.isFiltered()\'";

$status = system($commandline);
```

At this point I need to also exclude sites that have an excess of heterozygosity.  I do not want to exclude any site that departs from HWE because, due to the Wahlund effect, we expect many sites to have homozygote excess. The excess heterozygosity is likely due to duplicated regions that are slightly diverged.  We can identify the sites with excess heterozygosity like this:

for unzipped files:
```
vcftools --vcf round2_BADSEX_marked.vcf --hardy --out round2_BADSEX_marked
```


which generates a file called `round2_BADSEX_marked.hwe`, the last column of which has the probability of heterozygote excess.  We can get the entries that are less than 0.001 using this command:

```
awk -v OFS='\t' '(NR!=1) && ($8 < 0.005 ) {print $1,$2}' round2_BADSEX_marked.hwe >  bad_het_sites.txt
```

This command outputs a tab-delimted file with two columns including only the chromosome and position.  The header of the original file is skipped and not printed to the output file

Now we can use vcftools to make a vcf file with sites that have excess heterozygotes.  

for unzipped files:
```
vcftools --vcf round2_BADSEX_marked.vcf --positions bad_het_sites.txt --out bad_het_sitez --recode
```

This can be merged with the `round2_BADSEX_only.vcf` file we generated above and then used with the next per script to exclude these sites plus a buffer.



```bash
/usr/local/tabix/bgzip round2_BADSEX_only.vcf
/usr/local/tabix/tabix -p vcf round2_BADSEX_only.vcf.gz
/usr/local/tabix/bgzip bad_het_sitez.recode.vcf
/usr/local/tabix/tabix -p vcf bad_het_sitez.recode.vcf.gz
/usr/local/vcftools/src/perl/vcf-concat round2_BADSEX_only.vcf.gz bad_het_sitez.recode.vcf.gz | gzip -c > bad_sex_bad_het.vcf.gz
gunzip bad_sex_bad_het.vcf.gz
java -jar ~/picard-tools-1.131/picard.jar SortVcf I=bad_sex_bad_het.vcf O=bad_sex_bad_het_sorted.vcf SEQUENCE_DICTIONARY=/home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta.dict
mv bad_sex_bad_het_sorted.vcf bad_sex_bad_het.vcf
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
my $indel_only_file = "round2_indels_only.vcf";
my $bad_sex_file = "bad_sex_bad_het.vcf";
my $outfile1 = "temp1.vcf";
my $outfile2 = "temp2.vcf";
my $outfile3 = "round2_filtered.vcf";

# Mark sites that map poorly and that are in or within 3 bp of an indel
my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile1." --variant ".$file." --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\"";
$commandline = $commandline." --mask ".$indel_only_file." --maskName INDEL --maskExtension 3";
$status = system($commandline);

# Mark sites that have inappropriate sex chromosome genotypes
$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile2." --variant ".$outfile1;
$commandline = $commandline." --mask ".$bad_sex_file." --maskName BADSEX_badhet --maskExtension 200";
print $commandline,"\n";
$status = system($commandline);

# Print out the filtered file
$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile3." --variant ".$outfile2." -select 'vc.isNotFiltered()'";
$status = system($commandline);

# now clean up the intermediate vcf files
#$commandline = "rm -f ".$outfile1." ".$outfile2;
#$status = system($commandline);
```
OK, I am going to divide the multi-individual vcf file into individual vcf files, and then filtered individual genotypes on the basis of the depth of coverage and sex.  This will include using `awk` to insert a missing (./.) genotype for any diploid locus with less than 4X coverage and any haploid locus with less than 2X coverage.  Then I combine the filtered files back together.  This procedure allows me to not completely delete sites that need a site filteted in only one or a few individuals. The GATK option `SetFilteredGtToNoCall` did not work, which is why I used `awk`. Here is that script (14_Executes_GATK_commands_divide_mask_filter_merge):

``` perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will divide a concatenated vcf file 40 new vcf files 
# for each individual.  It will ten filter each based on depth
# and then merge them again.  I need to do individual filtering 
# because I do not want to delete an entire site because one individual
# has low coverage.

my $status;
my @files;
my @file_names;
my $n=0;
my $concat_vcf = "round2_filtered.vcf";
my $commandline;
my @females=("brunescens_PF707_stampy_sorted","hecki_PF643_stampy_sorted","hecki_PF644_stampy_sorted","hecki_PF648_stampy_sorted","hecki_PF651_stampy_sorted","maura_PF615_stampy_sorted","maura_PF713_stampy_sorted","nem_Gumgum_stampy_sorted","nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted","nem_pagensis_stampy_sorted","nigra_PF1001_stampy_sorted","nigra_PF660_stampy_sorted","nigrescens_PF654_stampy_sorted","ochreata_PF625_stampy_sorted","togeanus_PF549_stampy_sorted","tonk_PF515_stampy_sorted");
my @males=("hecki_PM639_stampy_sorted","hecki_PM645_stampy_sorted","maura_PM613_stampy_sorted","maura_PM614_stampy_sorted","maura_PM616_stampy_sorted","maura_PM618_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted","nigra_PM1000_stampy_sorted","nigra_PM1003_stampy_sorted","ochreata_PM571_stampy_sorted","ochreata_PM596_stampy_sorted","togeanus_PM545_stampy_sorted","tonk_PM561_stampy_sorted","tonk_PM565_stampy_sorted","tonk_PM566_stampy_sorted","tonk_PM567_stampy_sorted","tonk_PM582_stampy_sorted","tonk_PM584_stampy_sorted","tonk_PM592_stampy_sorted","tonk_PM602_stampy_sorted");
my $minimum_depth_of_coverage_per_diploid_genotype=4;
my $minimum_depth_of_coverage_per_haploid_genotype=2;
   
@files = glob("*_stampy_sorted.bam");

foreach(@files){
	$file_names[$n]=substr($files[$n], 0, -4);
	$n++;
}

# divide up the concatenated vcf file (uncomment for rerun)
foreach(@file_names){
	$commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta ";
	$commandline = $commandline." -V ".$concat_vcf." ";
	$commandline = $commandline." -sn ".$_." ";
	$commandline = $commandline." -o ".$_.".vcf";
	$status = system($commandline);
}

foreach(@file_names){
	if(in_array( \@females, $_) == 1){
		# filter female files; all sites are diploid so use the minimum criteria 
		$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";
		$commandline = $commandline." -o ".$_."_marked.vcf --variant ".$_.".vcf ";
		$commandline = $commandline." --filterExpression \"DP < ".$minimum_depth_of_coverage_per_diploid_genotype."\" --filterName \"LowDCoverage\"";
		$status = system($commandline);
		$commandline = "awk -v OFS=\'\\t\' \'{ if (\$7 == \"LowDCoverage\") \$10=\"./.\"; print \$0 }\' ".$_."_marked.vcf > ".$_."_filtered.vcf";
		$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
		$commandline = $commandline." -o ".$_."_filtered.vcf --variant ".$_."_marked.vcf";
		$status = system($commandline);
	}
	elsif(in_array( \@males, $_) == 1){
		# filter male files; all sites are diploid so use the minimum criteria 
		$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";
		$commandline = $commandline." -o ".$_."_marked.vcf --variant ".$_.".vcf ";
		$commandline = $commandline." --filterExpression \"DP < ".$minimum_depth_of_coverage_per_diploid_genotype." && CHROM != \'chrX\' && CHROM != \'chrY\'\" --filterName \"LowDCoverage\"";
		$commandline = $commandline." --filterExpression \"DP < ".$minimum_depth_of_coverage_per_haploid_genotype." && CHROM == \'chrX\'\" --filterName \"LowXCoverage\"";
		$commandline = $commandline." --filterExpression \"DP < ".$minimum_depth_of_coverage_per_haploid_genotype." && CHROM == \'chrY\'\" --filterName \"LowYCoverage\"";		
		$status = system($commandline);
		$commandline = "awk -v OFS=\'\\t\' \'{ if ((\$7 == \"LowDCoverage\")||(\$7 == \"LowXCoverage\")||(\$7 == \"LowYCoverage\")) \$10=\"./.\"; print \$0 }\' ".$_."_marked.vcf > ".$_."_filtered.vcf";
		$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
		$commandline = $commandline." -o ".$_."_filtered.vcf --variant ".$_."_marked.vcf";
		$status = system($commandline);
	}
	else{
		print "Problem with vcf file ".$_."\n";
	}	
}

# cleanup temporary files
foreach(@file_names){
 	$commandline = "rm -f ".$_."_marked.vcf";
 	$status = system($commandline);
}

# merge files

		$commandline = "java -Xmx5g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T CombineVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";
		foreach(@file_names){
			$commandline = $commandline." --variant ".$_."_filtered.vcf";
		}	
		$commandline = $commandline." -o final_filtered.vcf";		
		$status = system($commandline);

# cleanup more temporary files                                                                                                                                                                                                          
foreach(@file_names){
    $commandline = "rm -f ".$_.".vcf";
    $status = system($commandline);
    $commandline = "rm -f ".$_."_filtered.vcf";
    $status = system($commandline);

}

sub in_array {
my ($arr,$search_for) = @_;
foreach my $value (@$arr) {
	return 1 if $value eq $search_for;
}
 	return 0;
}
```


