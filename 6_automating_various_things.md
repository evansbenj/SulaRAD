# Automating addition of outgroups and calculation of popgen stats for each species

The script that adds the outgroup species can be executed for each tab file (there are 3 - one including genes plus 1000 bp on each end, one from 1001-51000 bp from genes, and one with the other bits). This needs to be done first for baboons, and then for humans and then a sed command to fix the header.  And for each of the three recal and non recal files (6 total). This can be done by pipine moving the tab delimited files into a folder and then globing them into the script I wrote for adding the outgroup.

And then we need to calculate the popgen stats and boostrap CIs for each of the species/populations for each of these files.  There are the following species:
* tonkeana
* hecki
* nigra
* nigrescent
* maura
* togeanus
* brunnescens
* ochreata
* nemestrina (Borneo)
* nemestrina (Sumatra)
* pagensis
* nemestrina (Malaysia)

I plan to do this using only the humans as outgroups, the justification being that I can get divergence data for the Y chromosome as well. Using the baboon I would not be able to do this because there is no baboon yDNA sequence.

