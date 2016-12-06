# Genotyping whole mtDNA genome

I plan to use the highest depth approach as a proof of principle for genotyping the mtDNA genomes of the three HiSeqX samples. The accuracy can be checked by looking for stop codons.

Here's how I extracted the mtDNA genotypes from the gvcf files André made:

```
~/tabix-0.2.6/tabix -h Sample_nemestrina-PM664/deliverables/nemestrina-PM664.g.vcf.gz chrM > nemestrina-PM664.g.vcf.gz_chrM.g.vcf
~/tabix-0.2.6/tabix -h Sample_nigra-PM664/deliverables/nigra-PM664.vcf.gz chrM > nigra-PF660.g.vcf.gz_chrM.g.vcf
~/tabix-0.2.6/tabix -h Sample_tonkeana-PM592/deliverables/tonkeana-PM592.g.vcf.gz chrM > tonkeana-PM592.g.vcf.gz_chrM.g.vcf
```
And then genotyping with `GenotypeGVCFs`:

```

```
