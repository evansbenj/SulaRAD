# HiSeqX X:A polymorphism in M. nigra

Reviewer 2 was interested in knowning what the X:A ratio was in the female HiSeqX sample - very reasonable request.  I did this as follows:

```
Boot_from_tab_diverge_poly_2015.pl /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf_males_highest_depth_females_byvcf.tab 010 3_4_2 PF660_HiSeqX_nigra_poly_and_diverge.txt
```
in this directory:

```
/home/ben/2015_SulaRADtag/good_merged_samples
```

Another question concerned the number of heterozygous genotype calls in males on the X in the HiSeqX data.  This is also a very interesting request.  

This information is in these files:
```
-rw-rw-r-- 1 ben ben     970213 Nov 29  2016 /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/hetsites_nemtonk_chrX.vcf.gz.tab
-rw-rw-r-- 1 ben ben 1657155767 Nov 27  2016 /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.gz.tab

```
