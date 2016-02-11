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
