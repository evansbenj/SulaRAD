# ABBABABA fig with inset

```R
# This R script will make a plot of f_dm from the Malinksy paper plotted
# against average pairwise divergence between H2 and H3 if f_dm is positive
# or between H1 and H3 if f_dm is negative. 
setwd('/projects/SulaRADTAG/perl_scripts/sliding_window/\ \ \ \ FINAL_ABBA/FINAL_WGS_new_depth_ABBABABA/all_data_heteroz_and_homoz_WGS')
library (ggplot2)
library(scales)
pdf("aDNA_and_xDNA_ABBABABA_plot.pdf",w=8, h=4, version="1.4", bg="transparent")
borneodata <- data.frame(matrix(NA, nrow=1, ncol=1))
switch<-0
borneodata <- borneodata[0,]
#borneodata<-read.table("concat_homoz_forplot.stats",header=TRUE)
borneodata<-read.table("concat_forplot.stats",header=TRUE)
# Make a new distance for aDNA and xDNA column that will have dH2H3 if f_dm is positive and dH1H3 if f_dm is negative
borneodata$conditional_distance <- 0
# make values of the conditional_distance column equal to dH2H3 if d_fm is positive
borneodata$conditional_distance[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH2H3) & borneodata$f_dm > 0)] <- borneodata$dH2H3[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH2H3) & borneodata$f_dm > 0)]
# make values of the conditional_distance column equal to dH1H3 if d_fm is negative
borneodata$conditional_distance[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH1H3) & borneodata$f_dm < 0)] <- borneodata$dH1H3[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH1H3) & borneodata$f_dm < 0)]
# make values of the conditional_distance column equal to the average of dH1H3 and dH2H3 if d_fm is equal to zero
borneodata$conditional_distance[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH1H3) & !is.na(borneodata$dH2H3) & borneodata$f_dm == 0)] <- (borneodata$dH1H3[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH1H3) & borneodata$f_dm == 0)]+borneodata$dH2H3[which(!is.na(borneodata$f_dm) & !is.na(borneodata$dH2H3) & borneodata$f_dm == 0)])/2
# Make a new column to distinguish aDNA and xDNA
borneodata$a_or_x <- 'black'
# make values of the this column equal to xDNA for xDNA
borneodata$a_or_x[which(borneodata$chromosome == 'chrX')]  <- 'red'
# Plot the data
d<-ggplot(data=borneodata, aes(x=conditional_distance,y=f_dm, colour=a_or_x)) +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  geom_point(size=4, alpha = 0.6) + 
  # label axis 
  labs(x=expression("Conditional distance to H3"), y=expression(paste(italic(f[DM])))) +
  # italicize the y-axis label
  #theme(axis.title.y=element_text(face="italic")) +
  # legend details
  scale_colour_manual(name="Genomic region", values = c("red"="red", "black"="gray"),breaks=c("black","red"),labels=c("Autosomes", "X chromosomes"))+
  # move the legend
  theme(legend.position = c(.8, .2)) +
  # remove boxes around legend symbols
  theme(legend.key = element_blank())+
  # draw a line on the yaxis=0
  geom_hline(yintercept=0) 
  #theme(text = element_text(size=8),axis.text.x = element_text(angle=90, vjust=1)) 
d

require(grid)
dat = read.table("H3born_H1nigra_H2tonk_chrX.stats", header=TRUE)
dat$Mbp <- dat$begin/1000000
scat<-ggplot(dat, aes(x=Mbp, y=f_dm)) +
  #d<-ggplot(dat, aes(x=begin, y=f_dm, colour=color, size = color, fill=factor(color))) +
  # make it clean
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  # label axis 
  labs(x=expression("Position (Mbp)"), y=expression(paste(italic(f[DM])))) +
  # remove the legend
  theme(legend.position="none") +
  # add points
  geom_point(colour = "red", size=2, alpha = 0.6) +
  # draw a line on the y axis
  geom_hline(yintercept=0.0) + 
  # set scale
  #scale_x_continuous(breaks=c(0,100000000,150000000),labels = comma) + 
  scale_x_continuous(labels = comma) + 
  # fill the dots differently
  #scale_fill_manual(values=c("gray46","royalblue1", "red1")) + 
  # make the outline of the mtDNA dots black
  #scale_color_manual(values=c("black","royalblue1", "red1"))+
  # make the mtdna larger
  #scale_size_manual(values=c(2,2,2)) +
  # plot in seperate facets
  #facet_grid(. ~ species) +
  # remove the strips on top
  theme(strip.background = element_blank(), strip.text.x = element_blank()) 

print(scat, vp=viewport(.5, .8, .45, .3))

dev.off()

```
