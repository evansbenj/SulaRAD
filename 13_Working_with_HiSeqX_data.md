# We now have HiSeqX data from M. nem PM664, M. tonkeana PM592, and M. nigra PM661; Cool!

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
then this might work as a bash script:
```
#!/bin/bash

~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr1 > HiSeqchr1.vcf 
~/tabix-0.2.6/bgzip HiSeqchr1.vcf
~/tabix-0.2.6/tabix -p vcf HiSeqchr1.vcf.gz
zcat HiSeqchr1.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr1.vcf.gz.tab
```

```
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr2 > HiSeqchr2.vcf # the -h flag preserves the header
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr3 > HiSeqchr3.vcf # the -h flag preserves the header
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr4 > HiSeqchr4.vcf # the -h flag preserves the header
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr5 > HiSeqchr5.vcf # the -h flag preserves the header
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr6 > HiSeqchr6.vcf # the -h flag preserves the header
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr7 > HiSeqchr7.vcf # the -h flag preserves the header
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr8 > HiSeqchr8.vcf # the -h flag preserves the header
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr9 > HiSeqchr9.vcf # the -h flag preserves the header
~/tabix-0.2.6/tabix -h nem_tonk_nigra_first10chrs.vcf.gz chr10 > HiSeqchr10.vcf # the -h flag preserves the header


~/tabix-0.2.6/bgzip HiSeqchr2.vcf
~/tabix-0.2.6/bgzip HiSeqchr3.vcf
~/tabix-0.2.6/bgzip HiSeqchr4.vcf
~/tabix-0.2.6/bgzip HiSeqchr5.vcf
~/tabix-0.2.6/bgzip HiSeqchr6.vcf
~/tabix-0.2.6/bgzip HiSeqchr7.vcf
~/tabix-0.2.6/bgzip HiSeqchr8.vcf
~/tabix-0.2.6/bgzip HiSeqchr9.vcf
~/tabix-0.2.6/bgzip HiSeqchr10.vcf


~/tabix-0.2.6/tabix -p vcf HiSeqchr2.vcf.gz
~/tabix-0.2.6/tabix -p vcf HiSeqchr3.vcf.gz
~/tabix-0.2.6/tabix -p vcf HiSeqchr4.vcf.gz
~/tabix-0.2.6/tabix -p vcf HiSeqchr5.vcf.gz
~/tabix-0.2.6/tabix -p vcf HiSeqchr6.vcf.gz
~/tabix-0.2.6/tabix -p vcf HiSeqchr7.vcf.gz
~/tabix-0.2.6/tabix -p vcf HiSeqchr8.vcf.gz
~/tabix-0.2.6/tabix -p vcf HiSeqchr9.vcf.gz
~/tabix-0.2.6/tabix -p vcf HiSeqchr10.vcf.gz


zcat HiSeqchr2.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr2.vcf.gz.tab
zcat HiSeqchr3.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr3.vcf.gz.tab
zcat HiSeqchr4.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr4.vcf.gz.tab
zcat HiSeqchr5.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr5.vcf.gz.tab
zcat HiSeqchr6.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr6.vcf.gz.tab
zcat HiSeqchr7.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr7.vcf.gz.tab
zcat HiSeqchr8.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr8.vcf.gz.tab
zcat HiSeqchr9.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr9.vcf.gz.tab
zcat HiSeqchr10.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > HiSeqchr10.vcf.gz.tab


```

And then use this script to add to this tab deimited file data from baboons (16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl). This is in this directory on info: `/home/ben/2015_SulaRADtag/baboon_rhesus_alignment` and is executed like this (from that directory):
```
/home/ben/2015_SulaRADtag/baboon_rhesus_alignment/16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/chr1_subset.vcf.gz.tab /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/chr1_subset.vcf.gz_with_baboon.tab
```
And use it again to add the human outgroup from this directory: `/home/ben/2015_SulaRADtag/axt_humans_rhemac2_from_UCSC` like this:
```
/home/ben/2015_SulaRADtag/axt_humans_rhemac2_from_UCSC/16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/chr1_subset.vcf.gz_with_baboon.tab /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/chr1_subset.vcf.gz_with_baboon_and_human.tab
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
