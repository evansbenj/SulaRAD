# ChrX

Nov 27, 2o16.  I am working with a vcf that I have already filtered to get rid of repeats.  It is in this directory:

```
/net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07
```
and this is the file name:
```
nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.all_with_one_highestdepth_allele.tab
```

This file was called for male and female individuals using ploidy = 2.

I have generated a tab file using highest depth for males and females from this file.  It is called:
```
nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.tab
```
The output from the ABBABABA script is called:
```
H1nigra_H2tonk_H3nem_chrX_highestdepth_M_and_F_norepeat.stats
```

I am also generating a diploid vcf file from the diploid calls.  It is called:
```
nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.gz.tab
```
I got the abbabbaba stats for this file like this:
```
./Performs_ABBA_BABA_on_populations_onlychrX.pl nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.gz.tab 010 3_4_1_2_3 H1nigra_H2tonk_H3nem_chrX_norepeat_ploidy2.jk H1nigra_H2tonk_H3nem_chrX_norepeat_ploidy2.stats
```

The vcf precursor to that tab file was used to identify heterozygous SNPs in males and then exclude them plus a mask of 3bp. Here is how I selected the heterozygous positions:
```
~/jre1.8.0_111/bin/java -Xmx2g -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta --variant nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.gz -select 'vc.getGenotype("nemestrina-PM664").isHet()' -select 'vc.getGenotype("tonkeana-PM592").isHet()' -o hetsites_nemtonk_chrX.vcf.gz
```
Here is how I marked and removed them plus a buffer of 3bp:
```
~/jre1.8.0_111/bin/java -Xmx2g -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -o nonrecal_filtered_chrX_final.vcf.gz_norepeat_malehets_are_marked.vcf.gz --variant nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.gz --mask hetsites_nemtonk_chrX.vcf.gz --maskName make_hets --maskExtension 3 

~/jre1.8.0_111/bin/java -Xmx2g -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta -o nonrecal_filtered_chrX_final.vcf.gz_norepeat_nomalehets.vcf.gz --variant nonrecal_filtered_chrX_final.vcf.gz_norepeat_malehets_are_marked.vcf.gz -select 'vc.isNotFiltered()'";
```

And the tab delimited file is:
```
nonrecal_filtered_chrX_final.vcf.gz_norepeat_nomalehgz.tab
```

I got the abbabbaba stats for this file like this:
```
./Performs_ABBA_BABA_on_populations_onlychrX.pl nonrecal_filtered_chrX_final.vcf.gz_norepeat_nomalehgz.tab 010 3_4_1_2_3 H1nigra_H2tonk_H3nem_chrX_norepeat_nomalehgz_ploidy2.jk H1nigra_H2tonk_H3nem_chrX_norepeat_nomalehgz_ploidy2.stats
```

 
I also generated a combined file with males (nem and tonk) genotyped with ploidy=1 and nigra with ploidy=2 (on iqaluk):

```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --variant /work/ben/2015_SulaRADtag/bam-constitutional/nemestrina-PM664.final.bam_haploid_X.g.vcf --variant /work/ben/2015_SulaRADtag/vcf-constitutional/nigra-PM664.g.vcf.gz_chrX.g.vcf --variant /work/ben/2015_SulaRADtag/bam-constitutional/tonkeana-PM592.final.bam_haploid_X.g.vcf --includeNonVariantSites -o nem_tonk_nigra_HiSeq_combined_chrX_nemtonk_ploidy1_nigra_ploidy2.vcf
```

This one needs to have the repeats removed and needs to be converted to a tab file. 



# Final files with repeat masker

So my preliminary analyses with the HiSeq and RADseq data suggested that there are some pretty serious issues with mis-mapped repeats in my data. For example, this may have caused many heterozygous sites in males on chrX.  The solution I came up with is to do the following
* Use the ploidy option in Haplotype caller to genotype the chrX of males (make it a haploid instead of a diploid)
* Filter out all repetitive elements from the genotype calls

I downloaded the chromOut.tar.gz file for rhemac2 and put this on info in this directory:
`/home/ben/2015_SulaRADtag/good_merged_samples/repeat_masker_chromOut/`

I also used a perl module called Number::Range to read in these ranges and check each site.  I had to do a local install of this module to run the perl script.  

```
mkdir ~/perl_modules
cd ~/perl_modules/Number-Range-0.12/
perl Makefile.PL INSTALL_BASE=~/perl_modules/Number-Range-0.12/
make
make install
```
And then I needed this in the top of the script:
`use lib qw(/home/ben/perl_modules/Number-Range-0.12/lib/);` to tell perl where to look

The final version of 'Performs_ABBA_BABA_on_populations_onlychrX.pl` is [here](https://github.com/evansbenj/SulaRAD/blob/master/12_More_on_visualizing_ABBABABA.md).

To run this script we need:
* A chromosome-specific tab-delimited genotype file
* the binary representation of the sex of the individuals (0 = male, 1 = female)
* the path and name of the repeatmasker file for the chr you are doing
* an output filename we don't use (`.jk`)
* an output filename we do use (`.stats`)


Here is an example of a chrX command on info:

```
/home/ben/2015_SulaRADtag/good_merged_samples/Performs_ABBA_BABA_on_populations_onlychrX.pl HiSeqchrX_ploidy.vcf.gz_with_baboon_and_humans.tab 010 3_6_1_2_3 /home/ben/2015_SulaRADtag/good_merged_samples/repeat_masker_chromOut/X/chrX.fa.out born_nigra_tonk_HiSeq_only_chrX.jk H3born_H1nigra_H2tonk_HiSeq_only_chrX_ploidy_norepeats.stats
```

I did the same for the autosomal analysis.  The final version of the script is [here](https://github.com/evansbenj/SulaRAD/blob/master/8__ABBA_BABA_on_populations.md).

Here is an example commandline:
```
/home/ben/2015_SulaRADtag/good_merged_samples/Performs_ABBA_BABA_on_populations.pl chr2_HiSeq_RADseq_combined.tab 1111100110000111100011100110010100000000010 3_6_14-41-19-20_32-33-34-35-36-37-38-43-40_22-42-25 /home/ben/2015_SulaRADtag/good_merged_samples/repeat_masker_chromOut/2/chr2.fa.out born_tonk_nigra_HiSeq_only_chr2.jk H3born_H1tonk_H2nigra_HiSeqRADseq_chr2_norepeats.stats
```

Turns out that this approach is extremely slow because the range finding perl function takes forever to load the ranges.  Instead I will just make bed files from the Repeatmasker files and then filter each vcf file before constructing the tab delimited files.  

I need a script that will take as input 3 vcf files and then do the following:
* output a subset each by chromosome 
* concatenate the three vcfs for each chr
* filter each of these vcf files using repeatmasker bed files
* output a filtered tab delimited file

# Make bed files out of repeatmasker files (24_Makes_bed_filez_out_of_repeatmasker.pl):

```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This script will read in a Repeatmasker file and make a bed file out of it

# 24_Makes_bed_filez_out_of_repeatmasker.pl repeatfile chr outfile.bed

# the second value will be the chr that is put in the first column

# ./24_Makes_bed_filez_out_of_repeatmasker.pl /projects/SulaRADTAG/perl_scripts/sliding_window/a_HiSeq_RADseq_combined_5mil_windows/repeat_masker_chromOut/4/chr4.fa.out chr4 /projects/SulaRADTAG/perl_scripts/sliding_window/a_HiSeq_RADseq_combined_5mil_windows/repeat_masker_chromOut/4/chr4.bed

my $inputfile = $ARGV[0];
unless (open DATAINPUT, $inputfile) {
    print 'Can not find the repeatmasker file.\n';
    exit;
}

my $chr = $ARGV[1];


my $outputfile = $ARGV[2];
unless (open(OUTFILE, ">$outputfile"))  {
    print "I can\'t write to $outputfile\n";
    exit;
}
my $line_number=0;
my @temp;

# print header
print OUTFILE "chrom\tchromStart\tchromEnd\n";


while ( my $line = <DATAINPUT>) {
    # begin by removing whitespaces from the silly format of Repeat masker
    $line =~ s/^\s+//;
    @temp=split(/\s+/,$line);
    $line_number+=1;
    if(($line_number == 4)&&(defined($temp[5]))&&(defined($temp[6]))){ # this is the first line
        # print to outfile
        print OUTFILE $chr,"\t",$temp[5]-1,"\t",$temp[6]-1,"\n"; # bed files start at zero
    }
    elsif(($line_number > 4)&&(defined($temp[5]))&&(defined($temp[6]))){
        # print to outfile
        print OUTFILE $chr,"\t",$temp[5]-1,"\t",$temp[6]-1,"\n";
    }
    # print counter
}   


close DATAINPUT;
close OUTFILE;
```
# Dividing up HiSeq gvcf file 

So after reading the GATK doc for the millionth time, it is clear that I need to keep the data in gvcf format and pipe it though a tool called GenotypeGVCFs to get a joint vcf file. This file can then be filtered to remove repeats and then converted to a tab delimited file.

Output a subset each by chromosome (GATK says to only use gvcf files generated by haplotypecaller, but since the original gvcf file was generated by haplotypecaller (by the NYGenomeCenter), this hopefully will be ok):

   tabix -p vcf myvcf.vcf.gz
   tabix -h myvcf.vcf.gz chr1 > chr1.vcf

Ok, now combine gVCFs for each autosome (on iqaluk; 25A_Combines_HiSeq_gVCFs.pl):

```perl
#!/usr/bin/perl                                                                                        
use warnings;
use strict;

# This script will merge aDNA gvcf files into a vcf file

my $status;
my $commandline;

#my @chr=("chr1","chr2");
my @chr=("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20");

foreach my $chr (@chr){
	$commandline = "/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --variant nemestrina-PM664.g.vcf.gz_".$chr.".g.vcf --variant tonkeana-PM592.g.vcf.gz_".$chr.".g.vcf --variant nigra-PM664.g.vcf.gz_".$chr.".g.vcf --includeNonVariantSites -o nem_tonk_nigra_HiSeq_combined_".$chr.".vcf";
	print $commandline,"\n";
	$status = system($commandline);
}

```

# Calling chrX on nem and tonk using ploidy 1 in Haplotypecaller

```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I /work/ben/2015_SulaRADtag/bam-constitutional/nemestrina-PM664.final.bam -ploidy 1 -out_mode EMIT_ALL_CONFIDENT_SITES --emitRefConfidence GVCF -L /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/target_interval_list_X.list -o nemestrina-PM664.final.bam_haploid_X.g.vcf

/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I /work/ben/2015_SulaRADtag/bam-constitutional/tonkeana-PM592.final.bam -ploidy 1 -out_mode EMIT_ALL_CONFIDENT_SITES --emitRefConfidence GVCF -L /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/target_interval_list_X.list -o tonkeana-PM592.final.bam_haploid_X.g.vcf

```

# Merge chrX HiSeq files (25B_Combines_HiSeq_gVCFs chrX.pl)
```perl
#!/usr/bin/perl                                                                                        
use warnings;
use strict;

# This script will merge xDNA gvcf files into a vcf file

my $status;
my $commandline;

#my @chr=("chr1","chr2");
my @chr=("chrX");

foreach my $chr (@chr){
	$commandline = "/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --variant /work/ben/2015_SulaRADtag/bam-constitutional/nemestrina-PM664.final.bam_haploid_X.g.vcf --variant /work/ben/2015_SulaRADtag/bam-constitutional/tonkeana-PM592.final.bam_haploid_X.g.vcf --variant /work/ben/2015_SulaRADtag/vcf-constitutional/nigra-PM664.g.vcf.gz_chrX.g.vcf --includeNonVariantSites -o nem_tonk_nigra_HiSeq_combined_".$chr.".vcf";
	print $commandline,"\n";
	$status = system($commandline);
}

```

### Filter merged vcfs (26_HiSeq_filters_combinedVCFs_chrX.pl)

```perl
#!/usr/bin/perl                                                                                        
use warnings;
use strict;

# This script will merge xDNA gvcf files into a vcf file

my $status;
my $commandline;

#my @chr=("chr1","chr2");
my @chr=("chrX");

foreach my $chr (@chr){
	$commandline = "/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --variant nem_tonk_nigra_HiSeq_combined_".$chr.".vcf -o nem_tonk_nigra_HiSeq_combined_".$chr."_marked.vcf --filterExpression \"QUAL < 30.0\" --filterName \"VeryLowQual\" -G_filter \"DP < 4\" -G_filterName \"LowCoverage\" --setFilteredGtToNocall";
	print $commandline,"\n";
	$status = system($commandline);

	$commandline = "/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -o nem_tonk_nigra_HiSeq_combined_".$chr."_filtered.vcf --variant nem_tonk_nigra_HiSeq_combined_".$chr."_marked.vcf --excludeFiltered --maxNOCALLfraction 0.99";
	print $commandline,"\n";
	$status = system($commandline);

}

```


### Filter merged vcfs and also avoid VariantFilter error (26_HiSeq_filters_combinedVCFs_chrX.pl) related to this part of the command: ` -G_filter \"DP < 4\" -G_filterName \"LowCoverage\" --setFilteredGtToNocall`

```perl
#!/usr/bin/perl                                                                                        
use warnings;
use strict;

# This script will merge xDNA gvcf files into a vcf file

my $status;
my $commandline;

#my @chr=("chr1","chr2");
my @chr=("chrX");

foreach my $chr (@chr){
	$commandline = "/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --variant nem_tonk_nigra_HiSeq_combined_".$chr.".vcf -o nem_tonk_nigra_HiSeq_combined_".$chr."_marked.vcf --filterExpression \"QUAL < 30.0\" --filterName \"VeryLowQual\"";
	print $commandline,"\n";
	$status = system($commandline);

	$commandline = "/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -o nem_tonk_nigra_HiSeq_combined_".$chr."_filtered.vcf --variant nem_tonk_nigra_HiSeq_combined_".$chr."_marked.vcf --excludeFiltered --maxNOCALLfraction 0.99";
	print $commandline,"\n";
	$status = system($commandline);

}

```



# Make a tab file
```
bgzip -c nem_tonk_nigra_HiSeq_combined_chrX_filtered.vcf > nem_tonk_nigra_HiSeq_combined_chrX_filtered.vcf.gz
tabix -p vcf nem_tonk_nigra_HiSeq_combined_chrX_filtered.vcf.gz
export PERL5LIB=/work/ben/vcftools/src/perl
zcat nem_tonk_nigra_HiSeq_combined_chrX_filtered.vcf.gz | /work/ben/vcftools/src/perl/vcf-to-tab > nem_tonk_nigra_HiSeq_combined_chrX_filtered.vcf.gz.tab
```

I suspect that the error with the genotypefilter above might be related to asterisks in the ALT field of vcf files.  This is a new [issue](http://gatkforums.broadinstitute.org/gatk/discussion/5904/asterisk-in-the-alt-feild-of-vcf). Anyhow, i still get the weird abba-baba with the haploid genotypes.  Now I am reanalyzing the data after converting the nigra seq to a haploid like this:

`sed -i 's/\/C/\//g' nem_tonk_nigra_HiSeq_combined_chrX_filtered3.vcf.gz.tab` and so on for each of the 4 bases.  I am also going to convert asterisks in this tab file (they are some, I've seen them, to a period like this:

`sed -i 's/\*/\./g' nem_tonk_nigra_HiSeq_combined_chrX_filtered3.vcf.gz.tab`





# BELOW THIS IS WRONG


This script divides up a vcf file from the NYGenome Center and then filters it using the bed files generated above (25_Splits_vcf_by_chr_and_filter_1.pl).


```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This script will take as input 3 vcf files and then do the following:

# output a subset each by chromosome
#   tabix -p vcf myvcf.vcf.gz
#   tabix -h myvcf.vcf.gz chr1 > chr1.vcf

# concatenate the three vcfs for each chr
#   export PATH=$PATH:~/tabix-0.2.6/
#   /usr/local/vcftools/src/perl/vcf-merge file1.vcf.gz file2.vcf.gz file2.vcf.gz | bgzip -c > file1_file2_file3_chr.vcf.gz

# filter each of these vcf files using repeatmasker bed files
#   /usr/local/vcftools/src/cpp/vcftools --gzvcf FILE --out OUTPUT_PREFIX --exclude-bed

# output a filtered tab delimited file for each chr

my $commandline;
my $combined_filename;
my $status;

#my @chr=("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chrX");
my @chr=("chr17","chr18","chr19","chr20","chrX");
my @vcf=("nemestrina-PM664.g.vcf.gz","nigra-PM664.g.vcf.gz","tonkeana-PM592.g.vcf.gz");


# export the path for tabix
#$commandline="export PATH=\$PATH:~/tabix-0.2.6/";
#$commandline="export PATH=\$PATH:/work/ben/vcftools/src/perl/";
#$ENV{PATH} = "$ENV{PATH}:/work/ben/vcftools/src/perl/";
# this worked for sharcnet:
# export PERL5LIB=/work/ben/vcftools/src/perl
#print $commandline,"\n";
#$status = system($commandline);


# index each vcf
#foreach my $vcf (@vcf){
#    $commandline="tabix -p vcf ".$vcf;
#    print $commandline,"\n";
#    #$status = system($commandline);
#}

# output a vcfsubset for each chr of each vcf file
#foreach my $vcf (@vcf){
#    foreach my $chr (@chr){
#        $commandline="tabix -h ".$vcf." ".$chr." > ".$vcf."_".$chr.".vcf";
#        print $commandline,"\n";
#        $status = system($commandline);
#    }
#}

# convert gvcfs to vcfs with one line per bp
foreach my $vcf (@vcf){
    foreach my $chr (@chr){
        $commandline="/work/ben/2015_SulaRADtag/gvcftools-0.16/bin/break_blocks --ref /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --region-file /work/ben/2015_SulaRADtag/vcf-constitutional/target_interval_list_allchrs.bed < ".$vcf."_".$chr.".vcf > ".$vcf."_".$chr.".noblock.vcf";
        print $commandline,"\n";
        $status = system($commandline);
    }
} 

# compress them
foreach my $vcf (@vcf){
    foreach my $chr (@chr){
        $commandline="bgzip -c ".$vcf."_".$chr.".noblock.vcf > ".$vcf."_".$chr.".noblock.vcf.gz";
        print $commandline,"\n";
        $status = system($commandline);
    }
} 

# index them
foreach my $vcf (@vcf){
    foreach my $chr (@chr){
        $commandline="tabix -p vcf ".$vcf."_".$chr.".noblock.vcf.gz";
        print $commandline,"\n";
        $status = system($commandline);
    }
} 


# merge the vcf files for each chr
foreach my $chr (@chr){
    #$commandline = "/usr/local/vcftools/src/perl/vcf-merge"; # for info
    $combined_filename=();
    $commandline = "/work/ben/vcftools/src/perl/vcf-merge"; # for sharcnet
        foreach my $vcf (@vcf){
            $commandline = $commandline." ".$vcf."_".$chr.".noblock.vcf.gz";
            $combined_filename=$combined_filename.$vcf;
        }
        $commandline = $commandline." | bgzip -c > ".$combined_filename."_".$chr.".vcf.gz";
        print $commandline,"\n";
        $status = system($commandline);
}


# filter each of the chr files
foreach my $chr (@chr){
    $commandline = "/work/ben/vcftools/src/cpp/vcftools --gzvcf ".$combined_filename."_".$chr.".vcf.gz --recode --out nemPM664_nigraPF660_tonkPM592_".$chr." --exclude-bed ./bed_files_from_repeatmasker/".$chr.".bed";
    print $commandline,"\n";
    $status = system($commandline);
}

# compress them
foreach my $chr (@chr){
        $commandline="bgzip -c nemPM664_nigraPF660_tonkPM592_".$chr.".recode.vcf > nemPM664_nigraPF660_tonkPM592_".$chr.".recode.vcf.gz";
        print $commandline,"\n";
        $status = system($commandline);
} 

# index them
foreach my $chr (@chr){
        $commandline="tabix -p vcf nemPM664_nigraPF660_tonkPM592_".$chr.".recode.vcf.gz";
        print $commandline,"\n";
        $status = system($commandline);
} 

# make tab delimited files
foreach my $chr (@chr){
    $commandline = "zcat nemPM664_nigraPF660_tonkPM592_".$chr.".recode.vcf.gz | /work/ben/vcftools/src/perl/vcf-to-tab > nemPM664_nigraPF660_tonkPM592_".$chr.".recode.vcf.gz.tab";
    print $commandline,"\n";
    $status = system($commandline);
}


```

# Using ploidy option with Haplotypecaller to call chrX on new hiseqX bam files

Today is Aug 26, 2016 and I am working with new bam files that have been BSQRed and indel realigned by the NYGenome center.  I already have vcf files they made but I want to recall chrX for nem and tonk using ploidy = 1 because they are males. I have the bam files in this directory in iqaluk: `/work/ben/2015_SulaRADtag/bam-constitutional` and I have an interval file for chrX in here: `/work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/target_interval_list_X.list`.

The ref genome is here: `/work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa`

vcf files will be in this directory on iqaluk: `/work/ben/2015_SulaRADtag/bam-constitutional`

```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I /work/ben/2015_SulaRADtag/bam-constitutional/nemestrina-PM664.final.bam -ploidy 1 -out_mode EMIT_ALL_CONFIDENT_SITES --emitRefConfidence BP_RESOLUTION -L /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/target_interval_list_X.list -o nemestrina-PM664.final.bam_haploid_X.vcf

/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa -I /work/ben/2015_SulaRADtag/bam-constitutional/tonkeana-PM592.final.bam -ploidy 1 -out_mode EMIT_ALL_CONFIDENT_SITES --emitRefConfidence BP_RESOLUTION -L /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/target_interval_list_X.list -o tonkeana-PM592.final.bam_haploid_X.vcf

```

To get vcftools to work on iqaluk, I had to type this:
```
export PERL5LIB=/work/ben/vcftools/src/perl
```

I needed to convert the nigra chrX  g.vcf to a vcf :

```
/work/ben/2015_SulaRADtag/gvcftools-0.16/bin/break_blocks --ref /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --region-file /work/ben/2015_SulaRADtag/vcf-constitutional/target_interval_list_allchrs.bed < nigra-PM664.g.vcf.gz_chrX.vcf > nigra-PM664.g.vcf.gz_chrX.noblock.vcf
```

After bgzip and tabix, I merged the chrX files:

```
/work/ben/vcftools/src/perl/vcf-merge ../bam-constitutional/nemestrina-PM664.final.bam_haploid_X.vcf.gz ../bam-constitutional/tonkeana-PM592.final.bam_haploid_X.vcf.gz nigra-PM664.g.vcf.gz_chrX.noblock.vcf.gz | bgzip -c > nemPM664tonkPM592nigraPF660_chrX_ploidy.vcf.gz
```

And then I removed the repetitive portions like this:

```
/work/ben/vcftools/src/cpp/vcftools --gzvcf nemPM664tonkPM592nigraPF660_chrX_ploidy.vcf.gz --recode --out nemPM664tonkPM592nigraPF660_chrX_ploidy_filtered.vcf.gz --exclude-bed /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/repeat_masker_chromOut/repeatmasker_bed_files/chrX.bed
```
