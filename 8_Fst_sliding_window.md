# Analysis of introgression across the Makassar Strait

I'd like to calculate FST in a sliging window across the genome for comparisons between Borneo nemestrina and Sulawesi macaques and for comparisons between Borneo nemestrina and individual Sulawesi macaques, such as tonkeana, hecki, and maura.  Vcftools looks like the most straightforward way to do this.  I first made a file containing names of all of the individual populations I wanted to compare.  These have a file name ending with "_names" in the `/home/ben/2015_SulaRADtag/good_merged_samples` directory on info.  Using the commands below, vcftools should generate a file with the suffix `weir.fst`.


```bash
~/tabix-0.2.6/bgzip final_round2_filtered.vcf
~/tabix-0.2.6/tabix -p vcf final_round2_filtered.vcf.gz
/usr/local/vcftools/src/cpp/vcftools --gzvcf final_round2_filtered.vcf.gz --weir-fst-pop borneo_names --weir-fst-pop tonkeana_names --fst-window-size 100000 --fst-window-step 100000
```


#ANGSD analysis

Also working on the ABBA BABA test.  Using angsd:

```bash
angsd -vcf-gl /home/ben/2015_SulaRADtag/good_merged_samples/final_round2_filtered.vcf.gz -fai /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta.fai -nind 40 -doAbbababa -anc /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta
```
or with bam files (problem is that they haven't been filtered)
```
angsd -out 665_602_545 -doAbbababa 1 -doCounts 1 -b nem665_tonk602_tog545_bamfilez -anc /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta
```
```
Rscript R/jackKnife.R file=665_602_545.abbababa indNames=nem665_tonk602_tog549_bamfilez outfile=665_602_545
```

or first by subsetting the vcf file:

```
/usr/local/vcftools/src/perl/vcf-subset -c nem_PM665_stampy_sorted,tonk_PM602_stampy_sorted,togeanus_PF549_stampy_sorted /home/ben/2015_SulaRADtag/good_merged_samples/final_round2_filtered.vcf.gz | ~/tabix-0.2.6/bgzip -c > 665_602_549.vcf.gz
```

# Making a bed file to use to filter a bam file

First I used a perl script to make a bed file from my vcf files:
```
./vcf2bed.pl bad_sex_bad_het.vcf > bad_sex_bad_het.bed
./vcf2bed.pl round2_indels_only.vcf > round2_indels_only.bed
```
then I had to remove an extra column:
```
awk '{print $1 "\t" $2 "\t" $3}' bad_sex_bad_het.bed > bad_sex_bad_het_.bed
awk '{print $1 "\t" $2 "\t" $3}' round2_indels_only.bed > round2_indels_only_.vcf
mv bad_sex_bad_het_.bed bad_sex_bad_het.bed
mv round2_indels_only_.vcf round2_indels_only.bed
```
then I made bed files for flanking regions:
```
~/bedtools2/bin/bedtools flank -b 200 -i bad_sex_bad_het.bed -g rhesus_chromosome_lengths > bad_sex_bad_het_100bp_buffer.bed
~/bedtools2/bin/bedtools flank -b 3 -i round2_indels_only.vcf -g rhesus_chromosome_lengths > round2_indels_only_3bp_buffer.bed
```

