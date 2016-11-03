# PCA

I am using the R package SNPrelate to generate PCA plots.  After lots of effort with ggplots, I've managed to generate exactly the plot I wanted to.  Here is the script for Sulawesi only samples:

```R
#https://github.com/zhengxwen/SNPRelate/issues/13
#http://corearray.sourceforge.net/tutorials/SNPRelate/
#https://github.com/zhengxwen/SNPRelate/wiki/Preparing-Data

#GenotypeVCFs_noBSQR_filtered_aDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_xDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_aDNA_only_no_lowcoverage_individuals.vcf

library("devtools")
library(gdsfmt)
library(SNPRelate)


#vcf.fn <- "/projects/SulaRADTAG/perl_scripts/pipeline/nonrecal_varonly_unifiedgenotyper_ym.vcf"
#vcf.fn <- "GenotypeVCFs_noBSQR_filtered_aDNA_only_no_lowcoverage_individuals.vcf.gz"
vcf.fn <- "GenotypeVCFs_noBSQR_filtered_aDNA_only_Sulawesi_only.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "test3.gds", method="biallelic.only")

genofile = openfn.gds("test3.gds", readonly=FALSE)


add.gdsn(genofile, "sample.annot", samp.annot)

snpgdsSummary("test3.gds")

pca <- snpgdsPCA(genofile, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#text(tab$EV2, tab$EV1,labels=tab$sample.id, cex= 0.4)

library(ggplot2)
#ggplot(...)+...+ theme(axis.text.x = element_text(angle=60, hjust=1))
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)

pdf("PCA_plot_all_aDNA_Sulawesi_only.pdf",w=8, h=8, version="1.4", bg="transparent")
tab$Species <- c("M. brunnescens", "M. hecki", "M. hecki", "M. hecki", "M. hecki", "M. hecki", "M. hecki", "M. maura", "M. maura", "M. maura", "M. maura", "M. maura", "M. maura", "M. nigra", "M. nigra", "M. nigrescens", "M. nigra", "M. nigrescens", "M. ochreata", "M. ochreata", "M. ochreata", "M. togeanus", "M. togeanus", "M. tonkeana", "M. tonkeana", "M. tonkeana", "M. tonkeana", "M. tonkeana", "M. tonkeana", "M. tonkeana", "M. tonkeana", "M. tonkeana")
tab$samp.color <- c("yellow", "green", "green", "green", "green", "green", "green", "orange", "orange", "orange", "orange", "orange", "orange", "purple", "purple", "pink", "purple", "pink", "brown", "brown", "brown", "gray", "gray", "black", "black", "black", "black", "black", "black", "black", "black", "black")
tab$samp.fieldid <- c("PF707","PF643","PF644","PF648","PF651","PM639","PM645","PF615","PF713","PM613","PM614","PM616","PM618","PF1001","PF660","PM1000","PM1003","PF654","PF625","PM571","PM596","PF549","PM545","PF515","PM561","PM565","PM566","PM567","PM582","PM584","PM592","PM602")

d<-ggplot(data=tab, aes(x=EV1,y=EV2, label = samp.fieldid, color = samp.color)) +
	# label axis 
 	labs(x=expression("Eigenvector 1"), y=expression("Eigenvector 2")) +
 	scale_colour_manual(name="Species", values = c("yellow"="yellow", "green"="green", "orange"="orange","purple"="purple", "pink"="pink","brown" = "brown", "gray" = "gray", "black" = "black"),breaks=c("yellow", "green", "orange","purple","pink","brown","gray","black"),labels=c("M.brunnescens", "M. hecki", "M. maura","M. nigra","M. nigrescens","M. ochreata","M. togeanus","M. tonkeana"))+
 	# add points and fieldID labels
	geom_text_repel(aes(EV1,EV2, label=(samp.fieldid))) + geom_point(size=3) + 
	# change to cleaner theme
	theme_classic(base_size = 16) +
	# make it clean
	theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
	# italicize species names
	theme(legend.text = element_text(face="italic"))+ 
	# move the legend
	theme(legend.position = c(.8, .7)) +
	# add a title
	ggtitle("PCA of Sulawesi Samples") 

dev.off()

```
