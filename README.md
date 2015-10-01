# SulaRAD
This is an electronic notebook for the Sulawesi macaque RADseq project

## Questions to be addressed
There are several major questions to be addressed;
* How much gene flow has occurred across the Makassar Strait since Sulawesi was colonized by macaques?  
* What are the relationships among the Sulawesi macaques?
* What is the ratio of effective population sizes of X:A in other macaques such as M. maura and M. tonkeana? Can we confirm the estimates we previously made for M. nemestrina and/or gain more insights from other species using data from just a few individuals (e.g. M. nemestrina from Sumatra, M. ochreata, M. nigra, M. nigrescens, M. "togeanus").

## The data
We have single end reads from 95-sample multipleded RADseq library that was run on one Illumina lane.  All samples were done with at least 2X replication and many were done with 3X or 4X replication.  We have 4-6 individuals for M. maura, M. hecki, and Borneo M. nemestrina, plus 1-2 for M. nigra, M. nigrescens, M. ochreata, M. brunescens, and M. togeanus.

I also have 1 lane of RADseq data paired end reads from 9 M. tonkeana individuals (1 female and 8 males).

## The approach
I have combined the forward read from the M. tonkeana data with the single end reads from the other newer data.  I have aligned these against the rhesus genome using `bwa` and `samtools` and done base recalibration with `GATK`. This resulted in ~60% of the reads mapping to the rhesus genome.

I am now going to try this with stampy.  First I need to format the rhesus genome with these two commands (on info):

`/usr/local/stampy/stampy.py -G /home/ben/2015_BIO720/rhesus_genome/rhesus_masked /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta`

this creates a .stidx file

the next step is to build a hash table for the reference genome using this command:

`/usr/local/stampy/stampy.py -g /home/ben/2015_BIO720/rhesus_genome/rhesus_masked -H /home/ben/2015_BIO720/rhesus_genome/rhesus_masked`

This should create a .sthash file.


