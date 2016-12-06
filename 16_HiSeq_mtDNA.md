# Genotyping whole mtDNA genome

I plan to use the highest depth approach as a proof of principle for genotyping the mtDNA genomes of the three HiSeqX samples. The accuracy can be checked by looking for stop codons.

I used a script called `3.05_Executes_GATK_commands_GenotypeVCFs_noBSQR_chrM.pl` to extract the mtDNA genotypes from the three gvcf files André made. This made a chrM file called `GenotypeVCFs_noBSQR_chrM.vcf.gz`.
