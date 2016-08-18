# We now have HiSeqX data from M. nem PM664, M. tonkeana PM592, and M. nigra PF661; Cool!

In addition to the raw data, we got bam files that were already processed by Andre Corvalo at the NY Genome Center.  They are in this directory `/home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27`.

We can assess depth like this
```
samtools depth ./Sample_nemestrina-PM664/analysis/alignment_corr/nemestrina-PM664.sorted.dedup.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```

This indicated that the coverage was >40X for each of the three samples.  Cool!

Use unified genotyper to call bases:
```
java -Xmx4G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I ./Sample_nemestrina-PM664/analysis/alignment_corr/nemestrina-PM664.sorted.dedup.bam  -I ./Sample_tonkeana-PM592/analysis/alignment_corr/tonkeana-PM592.sorted.dedup.bam  -I ./Sample_nigra-PM664/analysis/alignment_corr/nigra-PM664.sorted.dedup.bam -out_mode EMIT_ALL_CONFIDENT_SITES -o nem_tonk_nigra.vcf
```
or, on iqaluk, try this:

```
java -Xmx64G -jar /work/ben/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I ./Sample_nemestrina-PM664/analysis/alignment_corr/nemestrina-PM664.sorted.dedup.bam  -I ./Sample_tonkeana-PM592/analysis/alignment_corr/tonkeana-PM592.sorted.dedup.bam  -I ./Sample_nigra-PM664/analysis/alignment_corr/nigra-PM664.sorted.dedup.bam -nt 6 -nct 16 -out_mode EMIT_ALL_CONFIDENT_SITES -o nem_tonk_nigra.vcf
```
or for each chromosome on info try this on info:

```
java -Xmx32G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I ./Sample_nemestrina-PM664/analysis/alignment_corr/nemestrina-PM664.sorted.dedup.bam  -I ./Sample_tonkeana-PM592/analysis/alignment_corr/tonkeana-PM592.sorted.dedup.bam  -I ./Sample_nigra-PM664/analysis/alignment_corr/nigra-PM664.sorted.dedup.bam -out_mode EMIT_ALL_CONFIDENT_SITES -L target_interval_list_X.list -nct 3 -nt 8 -o nem_tonk_nigra_X.vcf
```

or for each chromosome on info try this on iqaluk:
```
java -Xmx64G -jar /work/ben/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I ./Sample_nemestrina-PM664/analysis/alignment_corr/nemestrina-PM664.sorted.dedup.bam  -I ./Sample_tonkeana-PM592/analysis/alignment_corr/tonkeana-PM592.sorted.dedup.bam  -I ./Sample_nigra-PM664/analysis/alignment_corr/nigra-PM664.sorted.dedup.bam -L target_interval_list_Y.list -nt 3 -nct 8 -out_mode EMIT_ALL_CONFIDENT_SITES -o nem_tonk_nigra_Y.vcf
```
with the `target_interval_list_X.list` file looking like this: `chrX:1-153947521`.

I also want to try setting the ploidy level of chrX to haploid for males like this on info:


```
java -Xmx32G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I ./Sample_nemestrina-PM664/analysis/alignment_corr/nemestrina-PM664.sorted.dedup.bam -ploidy 1 -out_mode EMIT_ALL_CONFIDENT_SITES -L target_interval_list_X.list -nct 3 -nt 8 -o nem_haploid_X.vcf

java -Xmx32G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I ./Sample_tonkeana-PM592/analysis/alignment_corr/tonkeana-PM592.sorted.dedup.bam -ploidy 1 -out_mode EMIT_ALL_CONFIDENT_SITES -L target_interval_list_X.list -nct 3 -nt 8 -o tonk_haploid_X.vcf

# don't change the default ploidy for nigraPF660 because it is female

java -Xmx32G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I ./Sample_nigra-PM664/analysis/alignment_corr/nigra-PM664.sorted.dedup.bam  -out_mode EMIT_ALL_CONFIDENT_SITES -L target_interval_list_X.list -nct 3 -nt 8 -o nigra_diploid_X.vcf

```

And then we need to merge the vcf files:

```
~/tabix-0.2.6/bgzip nem_haploid_X.vcf
~/tabix-0.2.6/tabix -p vcf nem_haploid_X.vcf.gz

~/tabix-0.2.6/bgzip tonk_haploid_X.vcf
~/tabix-0.2.6/tabix -p vcf tonk_haploid_X.vcf.gz

~/tabix-0.2.6/bgzip nigra_diploid_X.vcf
~/tabix-0.2.6/tabix -p vcf nigra_diploid_X.vcf.gz

export PATH=$PATH:~/tabix-0.2.6/
/usr/local/vcftools/src/perl/vcf-merge nem_haploid_X.vcf.gz tonk_haploid_X.vcf.gz nigra_diploid_X.vcf.gz | bgzip -c > nem_tonk_nigra_ploidy_chrX.vcf.gz
```
Now make a tab delimited file
```
~/tabix-0.2.6/tabix -p vcf nem_tonk_nigra_ploidy_chrX.vcf.gz
zcat nem_tonk_nigra_ploidy_chrX.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > nem_tonk_nigra_ploidy_chrX.vcf.gz.tab
```

While I am running this I will develop a pipeline to (1) add the outgroup sequences to a tab-deliminted file generated from the vsf file and (2) merge this with the existing RADseq file.  It will be interesting to quantify and compare how the genotypes from the radseq data compare to the HiSeq data for the same individuals.

OK, here we go.  First subset the vcf file in progress to creat an example file:
```
cat nem_tonk_nigra.vcf | awk 'NR >= 0 && NR <= 1000000 { print }' > nem_tonk_nigra_subset.vcf
```
and, because that only sampled chr1, I did this too:
```
cat nem_tonk_nigra.vcf | awk 'NR >= 569000000 && NR <= 570000000 { print }' > nem_tonk_nigra_subset_2.vcf
```
then concatenate them:
```
cat nem_tonk_nigra_subset.vcf nem_tonk_nigra_subset_2.vcf > nem_tonk_nigra_subsett.vcf
```
OK, I am going to split up the first ten chromosomes now.  I am working in this directory: `/net/infofile4-inside/volume1/scratch/ben`

Now make a tab delimited file:
```
~/tabix-0.2.6/bgzip nem_tonk_nigra_first10chrs.vcf
~/tabix-0.2.6/tabix -p vcf nem_tonk_nigra_first10chrs.vcf.gz
```
then for each chr automate generation of a tab file with a bash script:
```
#!/bin/bash

~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr1 > HiSeqchr1.vcf 
~/tabix-0.2.6/bgzip HiSeqchr1.vcf
~/tabix-0.2.6/tabix -p vcf HiSeqchr1.vcf.gz
zcat HiSeqchr1.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr1.vcf.gz.tab
```

Before generating the tab files, we need to filter based on coverage, HWE from RADseq data, and (for sex chromosomes) sex phenotypes. I already have this from the RADseq data in this directory:
`/home/ben/2015_SulaRADtag/good_merged_samples/round2_BADSEX_only.vcf` which includes indels, heterozygous sites, bad sex chromosome sites, plus a buffer of 3, 200, and 200 bp respectively as described [here] (https://github.com/evansbenj/SulaRAD/blob/master/4_More_GATK_Calling_All_Sites_and_Filtering.md).  For the autosomes this can be done by modifying this script `13_Executes_GATK_commands_VariantFiltration_doublemask.pl`.  But for the sex chromosomes, we will need to modify this file `12_Executes_GATK_commands_makes_sex_mask_VariantFiltration.pl` and make a new mask first, and then modifying and running `13_Executes_GATK_commands_VariantFiltration_doublemask.pl`.

I'm now thinking most or all of these acrobatics (apart from the sex chromosome checks) are not worthwhile because the number of filtered sites will be very small and the amount of data we have for sites that we cannot filter based on INDELs and HWE is vast.

Instead, let's just add the outgroup sequences to the tab delimited files. Goto this directory on info: `/home/ben/2015_SulaRADtag/baboon_rhesus_alignment` and execute this file (from that directory): (16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl) like this (for each chr). 
```
/home/ben/2015_SulaRADtag/baboon_rhesus_alignment/16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl /net/infofile4-inside/volume1/scratch/ben/HiSeqchr10.vcf.gz.tab /net/infofile4-inside/volume1/scratch/ben/HiSeqchr10.vcf.gz_with_baboon.tab
```
And then use it again to add the human outgroup from this directory: `/home/ben/2015_SulaRADtag/axt_humans_rhemac2_from_UCSC` like this:
```
/home/ben/2015_SulaRADtag/axt_humans_rhemac2_from_UCSC/16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl /net/infofile4-inside/volume1/scratch/ben/HiSeqchr10.vcf.gz_with_baboon.tab /net/infofile4-inside/volume1/scratch/ben/HiSeqchr10.vcf.gz_with_baboon_and_human.tab
```

Then modify the header like this:
```
sed -i -e 's/papAnu2\tpapAnu2/hg19\tpapAnu2/g' chr1_subset.vcf.gz_with_baboon_and_human.tab
```

OK now I should have 23 files (one for each of 20 autosomes, the X, Y, and mt), each with human, baboon and rhesus.  I'd like to write a script that combines these data with the RADseq data. The format should be identical to the RADseq data but with 3 extra columns on the left side.  And for many many sites, I will need to add lines with missing data in the RADseq columns and the HiSeq genotypes. I'd like to do this separately for each chromosome to facilatate parallel processing.

I need to make an example of the RADseq data:
```
cat final_round2_filtered.vcf.gz_with_baboon_and_human.tab | awk 'NR >= 0 && NR <= 1000000 { print }' > RADseq_subset.tab
```

I am working in a temp directory called `/net/infofile4-inside/volume1/scratch` because the vcf file is very large.

I need to do the following:
* split up the vcf file into chromosomes (I can do this for about half of them now using the vcf file `nem_tonk_nigra.vcf` which is still being made)
* for aDNA, filter each one based on coverage and remove excess heterozygous sites based on RADseq data
* for sex chromosomes do the same and also filter bad sex genotypes
* then make tab delimited files, add outgroup data, and combine them with same chromosome RADseq data for the ABBABABA analysis.  I wrote a script that will do the latter (`23_Combines_tab_delimited_files.pl`):
* 
```
#!/usr/bin/env perl
use strict;
use warnings;

# This script will read in a tab delimited file with whole genome seqs generated from HiSeqX
# and combine it with another tab delimited file with RADseq data.

# I'd like to put in checks to make sure that the reference genome bases are consistent
# and also quantify how frequently the genotypes of the RADseq and HiSeq data differ for
# the same individuals

# to run the script type this:
# 23_Combines_tab_delimited_files.pl HiSeqdata.tab RADseqdata.tab Combineddata.tab
# for example
# 23_Combines_tab_delimited_files.pl chr1_subset.vcf.gz_with_baboon_and_human.tab RADseq_subsett.tab chr1_combined.tab


my $inputfile = $ARGV[0];
unless (open DATAINPUT, $inputfile) {
    print 'Can not find the input file.\n';
    exit;
}

my $inputfile2 = $ARGV[1];
unless (open DATAINPUT2, $inputfile2) {
    print 'Can not find the input file.\n';
    exit;
}

my $outputfile = $ARGV[2];
unless (open(OUTFILE, ">$outputfile"))  {
    print "I can\'t write to $outputfile\n";
    exit;
}

my %HiSeq_genotypes;
my $number_of_hiseq_genotypes=3;
my $x;
my $y;
my @namez;
my @temp;
my $previous_chr='chr0';
my $previous_position=0;

# first read in the HiSeq data
while ( my $line = <DATAINPUT>) {
    chomp($line);
    @temp=split('\t',$line);
    if($temp[0] ne '#CHROM'){ # this assumes the HiSeq genotypes start at column 6 and that the
    	# first 3 columns are chr, pos, and bp of rhesus, the 4th is bp of human, and the 5th is bp of baboon 
    	for ($y = 2 ; $y <= $#temp ; $y++ ) {
    		push( @{$HiSeq_genotypes{$temp[0]."_".$temp[1]} }, $temp[$y]); 
    		# this will save the bp of rhesus, human, and bab, and also the other genos
    	}
    	if($previous_chr eq 'chr0'){
    		# this documents the first chr and position for 
    		# which we have HiSeq data
    		$previous_chr = $temp[0];
    		$previous_position = $temp[1];
    	}	
    }
    else{
 	   @namez=@temp;
 	   print "namez @namez\n";
    }
} 	
close DATAINPUT;

# now read in the RADseq data and simultaneously print out the output
# but do it only for the chr for which we have HiSeq data; this will be 
# done for each chr  
while ( my $line = <DATAINPUT2>) {
   chomp($line);
    @temp=split('\t',$line);
    if($temp[0] ne '#CHROM'){
    	if($HiSeq_genotypes{$temp[0]."_".$temp[1]}){
    		# this position has HiSeq data
    		# before printing out this line, print any HiSeq genotypes that were
    		# after the previous RADseq genotype but before the current one
    		if($previous_position < ($temp[1]-1) ){
    			# we first need to print hiSeq data before the current position $temp[1]
    			for ($y = $previous_position; $y < $temp[1] ; $y++ ) {
    				if($HiSeq_genotypes{$temp[0]."_".$y}){
    					# print rh_chr, rh_pos, rh_bp,hu_bp, bab_bp
    					print OUTFILE $temp[0],"\t",$y,"\t",$HiSeq_genotypes{$temp[0]."_".$y}[0],"\t",$HiSeq_genotypes{$temp[0]."_".$y}[1],"\t",$HiSeq_genotypes{$temp[0]."_".$y}[2];
    					# cycle through the RADseq data, assuming genotypes start at column 6
    					for ($x = 5 ; $x <= $#temp ; $x++ ) {
    						print OUTFILE "\t\.\/\.";
    					}
    					# now print the HiSeqgenotypes
    					for ($x = $#namez-$number_of_hiseq_genotypes+1; $x <= $#namez ; $x++ ) {
    						print OUTFILE "\t",$HiSeq_genotypes{$temp[0]."_".$y}[$x-2];
    					}
    					print OUTFILE "\n";
    				}	
    			}	
    		}
     		# now print out the RADseq data then the HiSeq genotypes
    		print OUTFILE $line,"\t";
			for ($y = ($#namez-$number_of_hiseq_genotypes+1); $y < $#namez ; $y++ ) {
				print OUTFILE $HiSeq_genotypes{$temp[0]."_".$temp[1]}[$y-2],"\t";
			}
			print OUTFILE $HiSeq_genotypes{$temp[0]."_".$temp[1]}[$#namez-2],"\n";	# because first 2 columns were used for key
    		$previous_position=$temp[1];
    	}
    	elsif($temp[0] eq $previous_chr){ # this RADseq position is missing from the HiSeq data 
    									  # but still on the HiSeq chr
    		print OUTFILE $line,"\t";
			for ($y = ($#namez-$number_of_hiseq_genotypes+1); $y < $#namez ; $y++ ) {
				print OUTFILE "\.\/\.\t";
			}
			print OUTFILE "\.\/\.\n";	
    	}
   		$previous_position=$temp[1];
   		# else, we are on a different RADseq chr, so print nothing
    }
    else{
    	print OUTFILE $line,"\t";
    	# print out the names of the HiSeq genotypes
    	for ($y = ($#namez-$number_of_hiseq_genotypes+1); $y < $#namez ; $y++ ) {
    		print OUTFILE $namez[$y],"\t";
    	}
    	print OUTFILE $namez[$#namez],"\n";	
    }	
}
close DATAINPUT2;
close OUTFILE;

```

This can be run on info like this:
```
./23_Combines_tab_delimited_files.pl HiSeqchr18.vcf.gz_with_baboon_and_human.tab /home/ben/2015_SulaRADtag/good_merged_samples/final_round2_filtered.vcf.gz_with_baboon_and_human.tab chr18_HiSeq_RADseq_combined.tab
```

OK now we seem to be on the way to having the alignments by chromosome for the combined RADseq and HiSeq data.  I want to run the abbababba script on these like this:
```
/home/ben/2015_SulaRADtag/good_merged_samples/Performs_ABBA_BABA_on_populations.pl chr9_HiSeq_RADseq_combined.tab 1111100110000111100011100110010100000000010 3_6_14-41-19-20_22-42-25_32-33-34-35-36-37-38-43-40 born_nigra_tonk_HiSeqRadchr9.jk H3born_H1nigra_H2tonk_HiSeqRadchr9.stats
```
When entered this way, H3 is nem born, H1 is nigra, and H2 is tonk. Negative f_dm indicates geneflow from nem to tonk; positive means geneflow from nem to nigra

I'd like to switch around who is H1 and H2 and also do smaller windows. On iqaluk:
```
Performs_ABBA_BABA_on_populations_1mil.pl chr1_HiSeq_RADseq_combined.tab 1111100110000111100011100110010100000000010 3_6_14-41-19-20_32-33-34-35-36-37-38-43-40_22-42-25 born_tonk_nigra_HiSeqRadchr9.jk H3born_H1tonk_H2nigra_HiSeqRadchr1_1mil.stats
```
on info:
```
/home/ben/2015_SulaRADtag/good_merged_samples/Performs_ABBA_BABA_on_populations_1mil.pl chr1_HiSeq_RADseq_combined.tab 1111100110000111100011100110010100000000010 3_6_14-41-19-20_32-33-34-35-36-37-38-43-40_22-42-25 born_tonk_nigra_HiSeqRadchr1.jk H3born_H1tonk_H2nigra_HiSeqRadchr1_1mil.stat
```

for chrX we need a separate script to deal with males and females differently:

```
/home/ben/2015_SulaRADtag/good_merged_samples/Performs_ABBA_BABA_on_populations_onlychrX_window_1mil.pl chrX_HiSeq_RADseq_combined.tab 1111100110000111100011100110010100000000010 3_6_14-41-19-20_32-33-34-35-36-37-38-43-40_22-42-25 born_tonk_nigra_HiSeqRadchrX.jk H3born_H1tonk_H2nigra_HiSeqRadchrX_1mil.stat
```
