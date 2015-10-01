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

I am now going to try this with stampy. To get stampy to work on info, we need to tell info to use version 2.7 of python like this:
`scl enable python27 bash`

Then we need to format the rhesus genome with these two commands (on info):

`/usr/local/stampy/stampy.py -G /home/ben/2015_BIO720/rhesus_genome/rhesus_masked /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta`

this creates a .stidx file

the next step is to build a hash table for the reference genome using this command:

`/usr/local/stampy/stampy.py -g /home/ben/2015_BIO720/rhesus_genome/rhesus_masked -H /home/ben/2015_BIO720/rhesus_genome/rhesus_masked`

This should create a .sthash file.

Then I made a bash script to do the alignments for lots of fastq files that looks like this:

```
#!/bin/bash                                                                                                                                                                                    

path_to_stampy="/usr/local/stampy"
path_to_data="."
path_to_genome="/home/ben/2015_BIO720/rhesus_genome"
genome="rhesus_masked"

individuals="brunescens_PF707                                                                                                                                                                  
hecki_PF643                                                                                                                                                                                    
hecki_PF644                                                                                                                                                                                    
hecki_PF648                                                                                                                                                                                    
hecki_PF651                                                                                                                                                                                    
hecki_PM639                                                                                                                                                                                    
hecki_PM645                                                                                                                                                                                    
maura_PF615                                                                                                                                                                                    
maura_PF713                                                                                                                                                                                    
maura_PM613                                                                                                                                                                                    
maura_PM614                                                                                                                                                                                    
maura_PM616                                                                                                                                                                                    
maura_PM618                                                                                                                                                                                    
nem_Gumgum                                                                                                                                                                                     
nem_Kedurang                                                                                                                                                                                   
nem_Malay                                                                                                                                                                                      
nem_Ngasang                                                                                                                                                                                    
nem_pagensis                                                                                                                                                                                   
nem_PM664                                                                                                                                                                                      
nem_PM665                                                                                                                                                                                      
nem_Sukai_male                                                                                                                                                                                 
nigra_PF1001                                                                                                                                                                                   
nigra_PF660                                                                                                                                                                                    
nigra_PM1000                                                                                                                                                                                   
nigra_PM1003                                                                                                                                                                                   
nigrescens_PF654                                                                                                                                                                               
ochreata_PF625                                                                                                                                                                                 
ochreata_PM571                                                                                                                                                                                 
ochreata_PM596                                                                                                                                                                                 
togeanus_PF549                                                                                                                                                                                 
togeanus_PM545"

for each_individual in $individuals
do

echo ${each_individual}
    $path_to_stampy/stampy.py -g $path_to_genome/$genome -h $path_to_genome/$genome -M $path_to_data/${each_individual}.fq
done

```


