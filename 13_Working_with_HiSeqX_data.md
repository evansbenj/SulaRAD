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

Now make a tab delimited file:
```
~/tabix-0.2.6/bgzip nem_tonk_nigra_subsett.vcf
~/tabix-0.2.6/tabix -p vcf nem_tonk_nigra_subsett.vcf.gz
tabix -h nem_tonk_nigra_subsett.vcf.gz chr1 > chr1_subset.vcf # this should be done for all of the chromosomes.  the -h flag preserves the header
~/tabix-0.2.6/bgzip chr1_subset.vcf
~/tabix-0.2.6/tabix -p vcf chr1_subset.vcf.gz

zcat chr1_subset.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > chr1_subset.vcf.gz.tab
```
And then use this script to add to this tab deimited file data from baboons (16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl). This is in this directory on info: `/home/ben/2015_SulaRADtag/baboon_rhesus_alignment` and is executed like this (from that directory):
```
/home/ben/2015_SulaRADtag/baboon_rhesus_alignment/16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/nem_tonk_nigra_subset.vcf.gz.tab /home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/nem_tonk_nigra_subset.vcf.gz_with_baboon.tab
```
And use it again to add the human outgroup from this directory: ``.  Then modify the header like this:
```
sed -i -e 's/papAnu2\tpapAnu2/hg19\tpapAnu2/g' nem_tonk_nigra_subset.vcf.gz_with_baboon_and_human.tab
```
