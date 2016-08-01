# We now have HiSeqX data from M. nem PM664, M. tonkeana PM592, and M. nigra PM661; Cool!

In addition to the raw data, we got bam files that were already processed by Andre Corvalo at the NY Genome Center.  From this directory `/home/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27`, I first merged them like this:

```
samtools merge nem_tonk_nigra_merged.bam ./Sample_nemestrina-PM664/analysis/alignment_corr/nemestrina-PM664.sorted.dedup.bam ./Sample_tonkeana-PM592/analysis/alignment_corr/tonkeana-PM592.sorted.dedup.bam ./Sample_nigra-PM664/analysis/alignment_corr/nigra-PM664.sorted.dedup.bam
```

Use unified genotyper to call bases (may have to add read groups first using picard)
```
java -Xmx4G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R XXX/ref.fasta -I samtools merge nem_tonk_nigra_merged.bam -out_mode EMIT_ALL_CONFIDENT_SITES -o nem_tonk_nigra.vcf
```
Make tab delimited files
```
~/tabix-0.2.6/bgzip nem_tonk_nigra.vcf
~/tabix-0.2.6/tabix -p vcf tnem_tonk_nigra.vcf.gz
zcat nem_tonk_nigra.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > nem_tonk_nigra.vcf.gz.tab
```
