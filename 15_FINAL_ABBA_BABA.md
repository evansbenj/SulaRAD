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
(1) output a subset each by chromosome 
(2) concatenate the three vcfs for each chr
(3) filter each of these vcf files using repeatmasker bed files
(4) output a filtered tab delimited file
