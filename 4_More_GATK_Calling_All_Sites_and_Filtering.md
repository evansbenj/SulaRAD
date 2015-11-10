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
OK, now using a perl script, I filtered individual genotypes on the basis of the depth of coverage.  This inserts a missing (./.) genotype for any diploid locus with less than 4X coverage and any haploid locus with less than 2X coverage.  This script has some additional checks for bad sex genotypes that *should* be unnecessary because they should have been filtered by the previous scripts. Here is that script (14_Vcf_filter.pl):

``` perl
#!/usr/bin/env perl
use strict;
use warnings;

# This program reads in a vcf file and filters individual genotypes by replacing their
# genotypes with "./." 
# This only affects the patricular genotype in question and will not change the other
# genotypes at that position
# This was not possible to do with GATK - once can filter by an individual genotype but 
# all of the genotypes appear to be removed at any position with one filter

# to execute type Vcf_filter.pl inputfile.vcf 1111100110000111100011100110010100000000 outputfile.vcf  
# where 1111100110000111100011100110010100000000 refers to whether or not each individual 
# in the vcf file is (1) or is not (0) female

# Vcf_filter.pl round2_filtered.vcf 1111100110000111100011100110010100000000 final_filtered.vcf

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[2];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";


my @sexes = split("",$ARGV[1]);

# Because these data were collected with two difference levels of coverage,
# I am not going to impose an upper limit on depth of coverage.  This was
# done in the 2014 paper but had not much effect.  Additionally by filtering 
# based on sites with many genotype calls and with all of them being heterozygous 
# we probably get rid of most repetitive mappings
my $minimum_depth_of_coverage_per_diploid_genotype=4;
my $minimum_depth_of_coverage_per_haploid_genotype=2;
my $number_of_ref_seq_microsats=0;
my $number_of_nonref_microsats=0;
my @male_chrX_hets=();
my @male_chrY_hets=();
my @male_chrX_lowcoverage=();
my @female_chrX_lowcoverage=();
my @male_chrY_lowcoverage=();
my @female_chrY=();
my @autosomal_lowcoverage=();
my @male_Yhet_sites;
my @male_Xhet_sites;
my @female_Y_sites;
my $y;
my @columns=();
my @DP;
my $DP;
my $GT;
my @genotypes;
my @tempp;

while ( my $line = <DATAINPUT>) {
	# print all commented lines to the outfile
	if(substr($line,0,1) eq "#"){
		print OUTFILE $line;
	}
	else{
		@columns=split(/\s/,$line);
			# filter microsatellites
			if((defined($columns[3]))&&(length($columns[3]) > 1)){
				# this is a microsat in the ref seq
				# delete the genotypes in the whole line
				# print the first 8 columns
				for ($y=0; $y<= 8; $y ++){
					print OUTFILE $columns[$y],"\t";
				}
				for ($y=0; $y<= ($#sexes); $y ++){
					print OUTFILE "./.\t";
				}
				print OUTFILE "./.\n";
				$number_of_ref_seq_microsats+=1;
			}
			elsif((defined($columns[3]))&&(length($columns[3]) == 1)){	
				# this is not a microsat in the ref seq
				# first find out where the depth of coverage is			
				@DP = split(":",$columns[8]);
				$DP=-1;
				$GT=-1;
				for ($y=0; $y<= $#DP; $y += 1){
					if($DP[$y] eq "DP"){
						$DP=$y;
					}
					if($DP[$y] eq "GT"){
						$GT=$y;
					}
				}
				# now print the first 8 columns
				for ($y=0; $y<= 8; $y ++){
					print OUTFILE $columns[$y],"\t";
				}
				# now filter based on coverage and location and print or delete each passing genotype
				for ($y=9; $y<= ($#sexes +9); $y ++){
					# check each individual genotype in the ingroup
					@genotypes=split(":",$columns[$y]);
					# if a call was made, check if either genotupe is more than one base
					if($#genotypes == 0){
						# no call made, print genotype
						if($y < ($#sexes + 9)){
							print OUTFILE "./.\t";
						}	
						else{
							print OUTFILE "./.\n";
						}	
					}
					elsif(($#genotypes > 0)&&(length($genotypes[0]) > 3)){ 
						# this is a microsat in the ingroup, so delete it
						if($y < ($#sexes + 9)){
							print OUTFILE "./.\t";
						}	
						else{
							print OUTFILE "./.\n";
						}
						$number_of_nonref_microsats+=1;	
					}
					elsif(($#genotypes > 0)&&(length($genotypes[0]) == 3)){
						# a non-microsat call was made, so now treat aDNA differently from xDNA
						if(($columns[0] eq "chrX") && ($sexes[$y-9] == 1)){
							# female genotype on the X
							if($genotypes[$DP] < $minimum_depth_of_coverage_per_diploid_genotype){
								# not enough coverage
								if($y < ($#sexes + 9)){
									print OUTFILE "./.\t";
								}	
								else{
									print OUTFILE "./.\n";
								}
								$female_chrX_lowcoverage[$y-9]+=1;
							}
							else{
								if($y < ($#sexes + 9)){
									print OUTFILE $columns[$y],"\t";
								}	
								else{
									print OUTFILE $columns[$y],"\n";
								}
							}
						}	
						elsif(($columns[0] eq "chrX") &&($sexes[$y-9] == 0)){
							# male genotype on the X
							if($genotypes[$DP] < $minimum_depth_of_coverage_per_haploid_genotype){
								# not enough coverage
								if($y < ($#sexes + 9)){
									print OUTFILE "./.\t";
								}	
								else{
									print OUTFILE "./.\n";
								}
								$male_chrX_lowcoverage[$y-9]+=1;
							}
							else{
								# check if the male is heterozygous on the X
								@tempp=split("/",$genotypes[$GT]);
								if(($tempp[0] ne $tempp[1])){
									# this is a heterozygous call on the X for a male, delete it
									if($y < ($#sexes + 9)){
										print OUTFILE "./.\t";
									}	
									else{
										print OUTFILE "./.\n";
									}
									$male_chrX_hets[$y-9]+=1;
									push(@male_Xhet_sites, $columns[1]);
								}
								else{
									if($y < ($#sexes + 9)){
										print OUTFILE $columns[$y],"\t";
									}	
									else{
										print OUTFILE $columns[$y],"\n";
									}
								}
							}
						}	
						if(($columns[0] eq "chrY") && ($sexes[$y-9] == 1)){
							# possible female genotype on the Y! Delete it.
							if($y < ($#sexes + 9)){
								print OUTFILE "./.\t";
							}	
							else{
								print OUTFILE "./.\n";
							}
							if($genotypes[$GT] ne "./."){
								$female_chrY[$y-9]+=1;
								push(@female_Y_sites, $columns[1]);
							}
						}	
						elsif(($columns[0] eq "chrY") && ($sexes[$y-9] == 0)){
							# male genotype on the Y
							if($genotypes[$DP] < $minimum_depth_of_coverage_per_haploid_genotype){
								# not enough coverage
								if($y < ($#sexes + 9)){
									print OUTFILE "./.\t";
								}	
								else{
									print OUTFILE "./.\n";
								}
								$male_chrY_lowcoverage[$y-9]+=1;
							}
							else{
								# check if the male is heterozygous on the Y
								@tempp=split("/",$genotypes[$GT]);
								if(($tempp[0] ne $tempp[1])){
									# this is a heterozygous call on the Y for a male, delete it
									if($y < ($#sexes + 9)){
										print OUTFILE "./.\t";
									}	
									else{
										print OUTFILE "./.\n";
									}
									$male_chrY_hets[$y-9]+=1;
									push(@male_Yhet_sites, $columns[1]);
								}
								else{
									if($y < ($#sexes + 9)){
										print OUTFILE $columns[$y],"\t";
									}	
									else{
										print OUTFILE $columns[$y],"\n";
									}
								}
							}
						}	
						else{
							# autosomal genotype
							if($genotypes[$DP] < $minimum_depth_of_coverage_per_diploid_genotype){
								# not enough coverage
								if($y < ($#sexes + 9)){
									print OUTFILE "./.\t";
								}	
								else{
									print OUTFILE "./.\n";
								}
								$autosomal_lowcoverage[$y-9]+=1;
							}
							else{
								if($y < ($#sexes + 9)){
									print OUTFILE $columns[$y],"\t";
								}	
								else{
									print OUTFILE $columns[$y],"\n";
								}
							}
						}	
					}
				} # end for loop over each genotype within a base position
			}# end of elsif
	} # end else	
}# end while
close DATAINPUT;
close OUTFILE;

# OK print out numbers
print "male_chrX_hets\n";

for ($y=0; $y<= $#sexes; $y ++){
	if(defined($male_chrX_hets[$y])){
		print $male_chrX_hets[$y]," ";
	}
	else{
		print "0 ";
	}
}
print "\n";
print "male_chrY_hets\n";
for ($y=0; $y<= $#sexes; $y ++){
	if(defined($male_chrY_hets[$y])){
		print $male_chrY_hets[$y]," ";
	}
	else{
		print "0 ";
	}
}
print "\n";
print "male_chrX_lowcoverage\n";
for ($y=0; $y<= $#sexes; $y ++){
	if(defined($male_chrX_lowcoverage[$y])){
		print $male_chrX_lowcoverage[$y]," ";
	}
	else{
		print "0 ";
	}
}
print "\n";
print "female_chrX_lowcoverage\n";
for ($y=0; $y<= $#sexes; $y ++){
	if(defined($female_chrX_lowcoverage[$y])){
		print $female_chrX_lowcoverage[$y]," ";
	}
	else{
		print "0 ";
	}
}
print "\n";
print "male_chrY_lowcoverage\n";
for ($y=0; $y<= $#sexes; $y ++){
	if(defined($male_chrY_lowcoverage[$y])){
		print $male_chrY_lowcoverage[$y]," ";
	}
	else{
		print "0 ";
	}
}
print "\n";
print "female_chrY\n";
for ($y=0; $y<= $#sexes; $y ++){	
	if(defined($female_chrY[$y])){
		print $female_chrY[$y]," ";
	}
	else{
		print "0 ";
	}
}
print "\n";
print "autosomal_lowcoverage\n";
for ($y=0; $y<= $#sexes; $y ++){
	if(defined($autosomal_lowcoverage[$y])){
		print $autosomal_lowcoverage[$y]," ";
	}
	else{
		print "0 ";
	}
}
print "\n";
print "male_Yhet_sites\n";
print "@male_Yhet_sites\n";
print "male_Xhet_sites\n";
print "@male_Xhet_sites\n";
print "female_Y_sites\n";
print "@female_Y_sites\n";


my $outputfile2 = "male_Xhet_sites.txt";
my $outputfile3 = "male_Yhet_sites.txt";
my $outputfile4 = "female_Y_sites.txt";
my @unique=();

unless (open(OUTFILE2, ">$outputfile2"))  {
    print "I can\'t write to $outputfile2\n";
    exit;
}
print "Creating output file: $outputfile2\n";

@unique = do { my %seen; grep { !$seen{$_}++ } @male_Xhet_sites };

foreach(@unique){
    print OUTFILE2 $_,"\n";
}

@unique=();


unless (open(OUTFILE3, ">$outputfile3"))  {
    print "I can\'t write to $outputfile3\n";
    exit;
}
print "Creating output file: $outputfile3\n";

@unique = do { my %seen; grep { !$seen{$_}++ } @male_Yhet_sites };

foreach(@unique){
    print OUTFILE3 $_,"\n";
}

unless (open(OUTFILE4, ">$outputfile4"))  {
    print "I can\'t write to $outputfile4\n";
    exit;
}
print "Creating output file: $outputfile4\n";

@unique = do { my %seen; grep { !$seen{$_}++ } @female_Y_sites };

foreach(@unique){
    print OUTFILE4 $_,"\n";
}

```

When this is done I am ready to convert the vcf file to a tab delimited file like this:
``` 
cd /work/ben/macaque_RAD_TAGs/tabix-0.2.6
./bgzip /path_to/final_filtered.vcf
./tabix -p vcf /path_to/final_filtered.vcf.gz
cd /work/ben/macaque_RAD_TAGs/vcftools_0.1.9/perl/ 
zcat /path_to/final_filtered.vcf.gz | ./vcf-to-tab > /path_to/final_filtered.vcf.gz.tab
```

And then use a script to add to this tab deimited file data from baboons.  Here is that script:

``` perl

```

