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

The final version of 'Performs_ABBA_BABA_on_populations_onlychrX.pl` is (https://github.com/evansbenj/SulaRAD/blob/master/12_More_on_visualizing_ABBABABA.md)[here].

```
/home/ben/2015_SulaRADtag/good_merged_samples/Performs_ABBA_BABA_on_populations_onlychrX.pl HiSeqchrX_ploidy.vcf.gz_with_baboon_and_humans.tab 010 3_6_1_2_3 /home/ben/2015_SulaRADtag/good_merged_samples/repeat_masker_chromOut/X/chrX.fa.out born_nigra_tonk_HiSeq_only_chrX.jk H3born_H1nigra_H2tonk_HiSeq_only_chrX_ploidy_norepeats.stats
```
