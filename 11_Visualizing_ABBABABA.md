# Visualizing ABBA_BABA results

OK, I am back from Ghana now and trying to get back into this analysis after a fairly long break.  Before I left I wrote this R script to plot the ABBABABA results (`X_to_A_sliding_fdm_outlier_only.R`) below. However, I just realized that the script before this (`Wrapper_for_Performs_ABBA_BABA_populations_most_inclusive.pl`) takes as input a command line that lists H3, then H1, then H2.  But this script outputs files with names in this order: H2, H1, H3.  This is very important becasue the order of the numbers in the file names defines how the ABBABABA stat is interpreted.  If correct, this means that the "j" loop below is for H1 and the i loop below is for H2.

``` R
# This R script will become an awesome plot with multiple panels that 
# illustrates the genomic regions with gene flow across the Makassar
# Strait. For each panel, I will have f_dm from the Malinksy paper plotted
# against average pairwise divergence between H2 and H3 if f_dm is positive
# or between H1 and H3 if f_dm is negative. Cool idea eh? And points with
# atypically high or low f_dm and atypically low pairwise divergence will be
# highlighted and labeled.  And there will be separate panels for comparisons
# using as H3 the Borneo samples, the Sumatra samples, and pagensis.

# (((H1,H2)H3)Outgroup)
# if f_dm is positive, this suggests gene flow from H2 and H3
# if f_dm is negative, this suggests gene flow from H1 and H3


setwd("/projects/SulaRADTAG/perl_scripts/sliding_window")
library (ggplot2)

# pdf("sliding_plot.pdf",w=8, h=4, version="1.4", bg="transparent")

# First define variables with the species names in it
species <- c("brun","hecki","maura","nigra","nigrescens","ochreata","togeanus","tonk");
nem = c("borneo", "sumatra", "pagensis")

# Define variables with the numbers that correspond to this species
species_numberz <- c("1","2-3-4-5-6-7","8-9-10-11-12-13","22-23-25","24-26","27-28-29","30-31","32-33-34-35-36-37-38-39-40");
nem_numberz = c("14-18-19-20","15-16-17","21")

# start here
borneodata <- data.frame(matrix(NA, nrow=1, ncol=1))
sumatradata <- data.frame(matrix(NA, nrow=1, ncol=1))
pagensisdata <- data.frame(matrix(NA, nrow=1, ncol=1))
alldata <- data.frame(matrix(NA, nrow=1, ncol=1))
switch<-0

# for brun as H1
#pdf("brunnescens_fdm_outlier_only.pdf")
#j=1
#for(i in 2:(length(species))) 

# for maura as H1
#pdf("maura_fdm_outlier_only.pdf")
#j=3
#for(i in c(1,2,4,5,6,7,8)) 

# for hecki as H1
#pdf("hecki_fdm_outlier_only.pdf")
#j=2
#for(i in c(1,3,4,5,6,7,8)) 

# for nigra as H1
#pdf("nigra_fdm_outlier_only.pdf")
#j=4
#for(i in c(1,2,3,5,6,7,8)) 
  
  
# for nigrescens as H1
#pdf("nigrescens_fdm_outlier_only.pdf")
#j=5
#for(i in c(1,2,3,4,6,7,8)) 

# for ochreata as H1
#pdf("ochreata_fdm_outlier_only.pdf")
#j=6
#for(i in c(1,2,3,4,5,7,8)) 

# for togeanus as H1
#pdf("togeanus_fdm_outlier_only.pdf")
#j=7
#for(i in c(1,2,3,4,5,6,8)) 
  

# for tonkeana as H1
pdf("tonkeana_fdm_outlier_only.pdf")
j=8
for(i in 1:(length(species)-1)) 
{
#  for(j in (i+1):length(species)) 
#  {
      borneodata <- borneodata[0,]
      sumatradata <- sumatradata[0,]
      pagensisdata <- pagensisdata[0,]
      borneodata<-read.table(paste(species_numberz[i],"_",species_numberz[j],"_14-18-19-20.stats",sep=""),header=TRUE)
      sumatradata<-read.table(paste(species_numberz[i],"_",species_numberz[j],"_15-16-17.stats",sep=""),header=TRUE)
      pagensisdata<-read.table(paste(species_numberz[i],"_",species_numberz[j],"_21.stats",sep=""),header=TRUE)


  
      #print(paste(k," ",i," ",j," ",species[i]," ",species[j] ))    
      
      # Read in the data
      #borneodata<-read.table("1_2-3-4-5-6-7_14-18-19-20.stats",header=TRUE)
      #sumatradata<-read.table("1_2-3-4-5-6-7_15-16-17.stats",header=TRUE)
      #pagensisdata<-read.table("1_2-3-4-5-6-7_21.stats",header=TRUE)
      
      # Get the cutoff number
      N_borneo <- as.integer(0.005*length(borneodata$chromosome[!is.na(borneodata$chromosome)]))
      N_sumatra <- as.integer(0.005*length(sumatradata$chromosome[!is.na(sumatradata$chromosome)]))
      N_pagensis <- as.integer(0.005*length(pagensisdata$chromosome[!is.na(pagensisdata$chromosome)]))
      
      print(paste("hello ",N_pagensis)) 
      
      # Make a new distance column that will have dH2H3 if f_dm is positive and dH1H3 if f_dm is negative
      borneodata$conditional_distance <- 0
      sumatradata$conditional_distance <- 0
      pagensisdata$conditional_distance <- 0
      
      # make values of the conditional_distance column equal to dH2H3 if d_fm is positive
      borneodata$conditional_distance[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH2H3) & borneodata$f_dm > 0)] <- borneodata$dH2H3[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH2H3) & borneodata$f_dm > 0)]
      sumatradata$conditional_distance[which(!is.na(sumatradata$f_dm) & !is.na(sumatradata$dH2H3) & sumatradata$f_dm > 0)] <- sumatradata$dH2H3[which(!is.na(sumatradata$f_dm) & !is.na(sumatradata$dH2H3) & sumatradata$f_dm > 0)]
      pagensisdata$conditional_distance[which(!is.na(pagensisdata$f_dm) & !is.na(pagensisdata$dH2H3) & pagensisdata$f_dm > 0)] <- pagensisdata$dH2H3[which(!is.na(pagensisdata$f_dm) & !is.na(pagensisdata$dH2H3) & pagensisdata$f_dm > 0)]
      
      # make values of the conditional_distance column equal to dH1H3 if d_fm is negative
      borneodata$conditional_distance[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH1H3) & borneodata$f_dm < 0)] <- borneodata$dH1H3[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH1H3) & borneodata$f_dm < 0)]
      sumatradata$conditional_distance[which(!is.na(sumatradata$f_dm) & !is.na(sumatradata$dH1H3) & sumatradata$f_dm < 0)] <- sumatradata$dH1H3[which(!is.na(sumatradata$f_dm) & !is.na(sumatradata$dH1H3) & sumatradata$f_dm < 0)]
      pagensisdata$conditional_distance[which(!is.na(pagensisdata$f_dm) & !is.na(pagensisdata$dH1H3) & pagensisdata$f_dm < 0)] <- pagensisdata$dH1H3[which(!is.na(pagensisdata$f_dm) & !is.na(pagensisdata$dH1H3) & pagensisdata$f_dm < 0)]
      
      # Get the lower divergence boundary
      #nth_smallest_conditional_distance_borneo <- order(borneodata$conditional_distance)[N_borneo]
      #nth_smallest_conditional_distance_sumatra <- order(sumatradata$conditional_distance)[N_sumatra]
      #nth_smallest_conditional_distance_pagensis <- order(pagensisdata$conditional_distance)[N_pagensis]
      
      # Get the upper f_dm boundary
      #nth_largest_f_dm_borneo <- order(borneodata$f_dm, decreasing = T)[N_borneo]
      #nth_largest_f_dm_sumatra <- order(sumatradata$f_dm, decreasing = T)[N_sumatra]
      #nth_largest_f_dm_pagensis <- order(pagensisdata$f_dm, decreasing = T)[N_pagensis]

      # Get the upper f_dm boundary for the absolute values
      nth_largest_f_dm_borneo <- order(abs(borneodata$f_dm), decreasing = T)[N_borneo]
      nth_largest_f_dm_sumatra <- order(abs(sumatradata$f_dm), decreasing = T)[N_sumatra]
      nth_largest_f_dm_pagensis <- order(abs(pagensisdata$f_dm), decreasing = T)[N_pagensis]
      print(paste("hello2 ",nth_largest_f_dm_pagensis, " value ",abs(pagensisdata$f_dm[nth_largest_f_dm_pagensis]))) 
            
      # Get the lower f_dm boundary
      #nth_smallest_f_dm_borneo <- order(borneodata$f_dm)[N_borneo]
      #nth_smallest_f_dm_sumatra <- order(sumatradata$f_dm)[N_sumatra]
      #nth_smallest_f_dm_pagensis <- order(pagensisdata$f_dm)[N_pagensis]
      
      # Make a new column withonly zeros that will be used for coloring and labeling points
      borneodata$m <- 0
      sumatradata$m <- 0
      pagensisdata$m <- 0
      
      # make values of the label column equal to one for high or low d_fm (conditional_distance can have any value)
      borneodata$m[!is.na(borneodata$f_dm) & !is.na(borneodata$conditional_distance) & abs(borneodata$f_dm) >= abs(borneodata$f_dm[nth_largest_f_dm_borneo])] <- 1
      #borneodata$m[!is.na(borneodata$f_dm) & !is.na(borneodata$conditional_distance) & borneodata$f_dm < borneodata$f_dm[nth_smallest_f_dm_borneo]] <- 1
      
      sumatradata$m[!is.na(sumatradata$f_dm) & !is.na(sumatradata$conditional_distance) & abs(sumatradata$f_dm) >= abs(sumatradata$f_dm[nth_largest_f_dm_sumatra])] <- 1
      #sumatradata$m[!is.na(sumatradata$f_dm) & !is.na(sumatradata$conditional_distance) & sumatradata$f_dm < sumatradata$f_dm[nth_smallest_f_dm_sumatra]] <- 1
      
      pagensisdata$m[!is.na(pagensisdata$f_dm) & !is.na(pagensisdata$conditional_distance) & abs(pagensisdata$f_dm) >= abs(pagensisdata$f_dm[nth_largest_f_dm_pagensis])] <- 1
      #pagensisdata$m[!is.na(pagensisdata$f_dm) & !is.na(pagensisdata$conditional_distance) & pagensisdata$f_dm < pagensisdata$f_dm[nth_smallest_f_dm_pagensis]] <- 1

      # make a column with blank labels for outlier points
      borneodata$labels <- ""
      sumatradata$labels <- ""
      pagensisdata$labels <- ""
      
      # now add a label to the outliers but not to the other points
      borneodata$labels[which(borneodata$m == "1")] <- paste(borneodata$chromosome[which(borneodata$m == "1")],borneodata$begin[which(borneodata$m == "1")],sep = "_")
      sumatradata$labels[which(sumatradata$m == "1")] <- paste(sumatradata$chromosome[which(sumatradata$m == "1")],sumatradata$begin[which(sumatradata$m == "1")],sep = "_")
      pagensisdata$labels[which(pagensisdata$m == "1")] <- paste(pagensisdata$chromosome[which(pagensisdata$m == "1")],pagensisdata$begin[which(pagensisdata$m == "1")],sep = "_")
      
      # Now add a column to specify the value of H3
      borneodata$H3 <- "borneo"
      sumatradata$H3 <- "sumatra"
      pagensisdata$H3 <- "pagensis"
      
      # Add a column to specify the H1 species comparison
      borneodata$comparison <- species[i]
      sumatradata$comparison <- species[i]
      pagensisdata$comparison <- species[i]

      # Combine the data
      if(switch==0){
        switch<-1
        alldata<-rbind(borneodata,sumatradata,pagensisdata)
      }
      else{
        alldata<-rbind(alldata,borneodata,sumatradata,pagensisdata) 
      }
#  } # this is for the j loop
} # this is for the i loop 

# to plot by chromosome
#data$chromosome <- factor(data$chromosome, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20"))
# d<-ggplot(data=data, aes(x=dH2H3,y=fdH2H3,group=chromosome, colour=data$m)) + geom_point(alpha = 1) + facet_wrap(~ chromosome, scales = "free_x") + theme(text = element_text(size=8),axis.text.x = element_text(angle=90, vjust=1)) + geom_text(aes(label=data$labels),hjust=0, vjust=0, size = 2)
# to plot all together

# Plot the data
d<-ggplot(data=alldata, aes(x=conditional_distance,y=f_dm, colour=m)) + geom_point(alpha = 0.5) + theme(text = element_text(size=8),axis.text.x = element_text(angle=90, vjust=1)) + geom_text(aes(label=labels),hjust=0, vjust=0, size = 2,check_overlap = TRUE, angle = 25)+facet_grid(~ H3) 
d + facet_wrap(comparison ~ H3, ncol=3, switch ="y")

#d<-ggplot(data=alldata, aes(x=conditional_distance,y=f_dm, colour=m)) +theme_bw() + geom_point(alpha = 0.5) + theme(text = element_text(size=8),axis.text.x = element_text(angle=90, vjust=1)) 
#d + facet_wrap(comparison ~ H3, ncol = 3, switch ="y") 

dev.off()



#+
#theme(strip.text.x = element_text(size=8, angle=75),
#      strip.text.y = element_text(size=12, face="bold"),
#      strip.background = element_rect(colour="red", fill="#CCCCFF"))


#outliers$chr <- data[data$m == 1,]$chromosome 
#outliers$begin <- data[data$m == 1,]$begin 
#outliers
#outliers

#dev.off()


```
