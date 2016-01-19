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
my $outfile2="round2_BADSEX_only.vcf";

my @females=("brunescens_PF707_stampy_sorted","hecki_PF643_stampy_sorted","hecki_PF644_stampy_sorted","hecki_PF648_stampy_sorted","hecki_PF651_stampy_sorted","maura_PF615_stampy_sorted","maura_PF713_stampy_sorted","nem_Gumgum_stampy_sorted","nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted","nem_pagensis_stampy_sorted","nigra_PF1001_stampy_sorted","nigra_PF660_stampy_sorted","nigrescens_PF654_stampy_sorted","ochreata_PF625_stampy_sorted","togeanus_PF549_stampy_sorted","tonk_PF515_stampy_sorted");

my @males=("hecki_PM639_stampy_sorted","hecki_PM645_stampy_sorted","maura_PM613_stampy_sorted","maura_PM614_stampy_sorted","maura_PM616_stampy_sorted","maura_PM618_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted","nigra_PM1000_stampy_sorted","nigra_PM1003_stampy_sorted","ochreata_PM571_stampy_sorted","ochreata_PM596_stampy_sorted","togeanus_PM545_stampy_sorted","tonk_PM561_stampy_sorted","tonk_PM565_stampy_sorted","tonk_PM566_stampy_sorted","tonk_PM567_stampy_sorted","tonk_PM582_stampy_sorted","tonk_PM584_stampy_sorted","tonk_PM592_stampy_sorted","tonk_PM602_stampy_sorted");

my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta";
$commandline = $commandline." -o ".$outfile1." --variant ".$file;

# filter all sites in which a female genotype is calledon chrY
foreach (@females){
    $commandline=$commandline." --filterExpression \"CHROM == \'chrY\'\" --filterName BAD_SEXfY";
}   

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

```
vcftools --gzvcf XXX.vcf.gz --hardy
```

or, for unzipped files:
```
vcftools --vcf XXX.vcf --hardy
```


which generates a file called `out.hwe`, the last column of which has the probability of heterozygote excess.  We can get the entries that are less than 0.001 using this command:

```
awk -v OFS='\t' '(NR!=1) && ($8 < 0.001 ) {print $1,$2}' out.hwe >  bad_hwe_sites.txt
```

This command outputs a tab-delimted file with two columns including only the chromosome and position.  The header of the original file is skipped and not printed to the output file

Now we can use vcftools to make a vcf file with sites that have excess heterozygotes.  

for unzipped files:
```
vcftools --vcf XXX.vcf --positions bad_hwe_sites.txt --out bad_hwe_sitez --recode
```

This can be merged with the `round2_BADSEX_only.vcf` file we generated above and then used with the next per script to exclude these sites plus a buffer.



``
bgzip round2_BADSEX_only.vcf
tabix -p vcf round2_BADSEX_only.vcf.gz
bgzip bad_hwe_sitez.vcf
tabix -p vcf bad_hwe_sitez.vcf.gz
vcf-concat round2_BADSEX_only.vcf.gz bad_sitez.vcf.gz | gzip -c > bad_sex_bad_hwe.vcf.gz
gunzip bad_sex_bad_hwe.vcf.gz
``

Now get rid of indels and bad_sex genotypes plus a buffer of 3 bp and 200 bp respectively (13_Executes_GATK_commands_VariantFiltration_doublemask.pl):
(still need to change the mask file below to "bad_sex_bad_hwe.vcf.gz"

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
In the 2014 analysis, at this stage I used bed files and GATK to output subsets of the data that spanned genes, that were near, genes, and that were progressively farther from genes.  For now, I'll skip this step and analyze the whole dataset (which is mostly far from genes anyhow).

When this is done I am ready to convert the vcf file to a tab delimited file like this:

``` 
~/tabix-0.2.6/bgzip final_filtered.vcf
~/tabix-0.2.6/tabix -p vcf final_filtered.vcf.gz
zcat final_filtered.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > final_filtered.vcf.gz.tab
```

And then use a script to add to this tab deimited file data from baboons (`15_Gets_outgroup_sequence_from_axt_files_NEW2015.pl`).  Here is that script:

``` perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program reads in a tab delimited file created by 
# vcftools and then extracts an outgroup sequence
# from all axt files that are in the directory

# It will then make a new file that has the outgroup
# sequence inserted in a column after the reference sequence.

# run it like this:
# Gets_outgroup_sequence_from_axt_files_NEW2015.pl in_tabfile out_tabfile

# the main concern here is that theaxt files have gaps inserted so that the 
# number of bases don't necessarily match the difference between the coordinates
# although hopefully the coordinates match.  So I need to count the gaps in the
# rhesus sequence and adjust the coordinates of the outgroup sequence appropriately

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";


my @files = glob("*.axt");

my @macaque_coordinates;
my @ingroup;
my @outgroup;
my $coordinate;
my %macaque_coordinates_key;
my $line;
my @line;
my @axt;
my $y;
my $switch;
my $source;
my $begin;
my $chr;


while ( my $line = <DATAINPUT>) {
	@macaque_coordinates=split("	",$line);
	if ($macaque_coordinates[0] ne "#CHROM"){
		#$macaque_coordinates_key{$macaque_coordinates[0]."_".$macaque_coordinates[1]}=$macaque_coordinates[2];
		$macaque_coordinates_key{$macaque_coordinates[0]."_".$macaque_coordinates[1]}="N";
	}
}

close DATAINPUT;
print "Done with input file 1\n";

foreach(@files){
	print $_,"\n";
	unless (open DATAINPUT1, $_) {
		print "Can not find the axt files, jackass.\n";
		exit;
	}

	LINE: while ( my $line1 = <DATAINPUT1>) {	
		@axt=split(" ",$line1);
		if(defined($axt[1])){
			if($axt[1] =~ /^chr/){
				$switch=1;
				$chr=$axt[1];
				$begin=$axt[2];
				next LINE;
			}
		}
		elsif((defined($axt[0]))&&($switch == 1)){
			# we need to find out where the gaps are in the reference sequence
			@ingroup=split("",$line1);	
			$switch = 2;
			next LINE;
		}
		elsif((defined($axt[0]))&&($switch == 2)){	
			$switch = 0;
			#print "papio ",$line1;
			@outgroup=split("",$line1);
			#print "yo ",$outgroup[0]," hey ",$#outgroup," ey ",$outgroup[$#outgroup-1],"\n";
			$coordinate=$begin-1;
			for ($y = 0 ; $y < $#outgroup ; $y++ ) {
				if($ingroup[$y] ne "-"){
					$coordinate+=1;
				}
				#print $chr."_".$coordinate,"\n";
				# check if this position is a gap in the ingroup
				if(defined($macaque_coordinates_key{$chr."_".$coordinate})){
					#print $chr,"_",$coordinate,"  ",$macaque_coordinates_key{$chr."_".$coordinate}," ";
					if($outgroup[$y] ne "-"){
						$macaque_coordinates_key{$chr."_".$coordinate} = $outgroup[$y];
					}
					else{
						$macaque_coordinates_key{$chr."_".$coordinate} = "N";
					}	
					#print $macaque_coordinates_key{$chr."_".$coordinate},"\n";
					if($y == 0){
						print "begin ",$chr,"_",$coordinate,"  ",$macaque_coordinates_key{$chr."_".$coordinate}," ";
						print $macaque_coordinates_key{$chr."_".$coordinate},"\n";
					}
					elsif($y == ($#outgroup)){
						print "end ",$chr,"_",$coordinate,"  ",$macaque_coordinates_key{$chr."_".$coordinate}," ";
						print $macaque_coordinates_key{$chr."_".$coordinate},"\n";
					}
				}
			}	
		}
	}
close DATAINPUT1;	
}

# now add the outgroup data to the output file in a new column
my @data; 
unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}
while ( my $line = <DATAINPUT>) {
	@data=split(/\s+/,$line);
	if ($data[0] ne "#CHROM"){
		for ($y = 0 ; $y <= 2 ; $y++ ) {
			print OUTFILE $data[$y],"\t";
		}
		if(defined($macaque_coordinates_key{$data[0]."_".$data[1]})){
			print OUTFILE $macaque_coordinates_key{$data[0]."_".$data[1]},"\t";	
		}
		else{
			print "Missed one; problem!\n";	
			print OUTFILE "N\t";	
		}	
		for ($y = 3 ; $y < $#data ; $y++ ) {
			print OUTFILE $data[$y],"\t";
		}
		print OUTFILE $data[$#data],"\n";
	}
	else{
		for ($y = 0 ; $y <= 2 ; $y++ ) {
			print OUTFILE $data[$y],"\t";
		}
		print OUTFILE "papAnu2\t";	
		for ($y = 3 ; $y < $#data ; $y++ ) {
			print OUTFILE $data[$y],"\t";
		}
		print OUTFILE $data[$#data],"\n";
	}
}


```
This works quite well and I have also used it to add the human outgroup sequence.  A small correction to the header is then required using sed:

`sed -i -e 's/papAnu2\tpapAnu2/hg19\tpapAnu2/g' final_round2_filtered.vcf.gz_with_baboon_and_human.tab`



