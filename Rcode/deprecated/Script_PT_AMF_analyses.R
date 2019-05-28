##April 3rd 2017
##Author: Jacynthe Masse
##R script for Mohammed Hijri's project: Premier Tech
##Last update: April 10th 2017

# Libraries and functions -------------------------------------------------

library(vegan)
library("pgirmess")  #For post hoc test with KW
library("ggplot2")
library("ape")
library(reshape2)
source("anova.1way.R")
source("summarySE.R")
source("anova.2way.R")
source("norm.her.R")
source("homeo.table.R")
source("relativeAbundance_otu.R")
source("anova.1way.table.R")

# Databases ---------------------------------------------------------------

rarefaction<-read.delim("PT_AMF_clustered_rarefaction.txt", header=T, row.names=1)
str(rarefaction)
dim(rarefaction)

alphaD<-read.delim("PT_AMF_clustered_alphaD.txt", dec=".", row.names=1, header=T)
str(alphaD)
dim(alphaD)

PT_OTU<-read.delim("PT_AMF_clustered_OTU_R.txt", row.names=1, header=T)
str(PT_OTU)
str(PT_OTU$Crop)
str(PT_OTU$Inoculation)
dim(PT_OTU)

PT_OTU_totransf<-read.delim("PT_AMF_clustered_OTU_fortrans_R.txt", row.names=1, header=T)


PT_vtx<-read.delim("PT_AMF_VTX_MA_RA_R.txt", row.names=1, header=T, dec=",")
str(PT_vtx)

PT_orders<-read.delim("PT_AMF_orders_MA_RA_R.txt", row.names=1, header=T, dec=",")
str(PT_orders)

PT_families<-read.delim("PT_AMF_families_MA_RA.txt", row.names=1, header=T, dec=",")
str(PT_families)

PT_genera<-read.delim("PT_AMF_genera_MA_RA_R.txt", row.names=1, header=T, dec=",")
str(PT_genera)

PT_colonization<-read.delim("PT_colonisation_R.txt", row.names=1, header=T, na.strings= "NA",  dec=",")
str(PT_colonization)

# Rarefaction curves ------------------------------------------------------
rarefaction<-read.delim("PT_AMF_clustered_rarefaction.txt", header=T, na.strings= "NA")

##With the basic package
(maxseq<-max(rarefaction$numsampled))
(maxotu<-max(rarefaction[,2:117], na.rm=TRUE))
cl<-rainbow(116)

plot(1, type='n', xlim=c(0,maxseq), ylim=c(0,maxotu), xlab="Number of sequences", ylab="Number of OTUs")

for (i in 2:ncol(rarefaction))
{
  lines(rarefaction$numsampled,rarefaction[,i], xlab="Number of sequences", ylab="Number of OTUs", col=cl[i], lwd=2)
}


# Bar-graph presenting the communities  ------------------------------------------------------------------
##VTX
library(reshape2)
PT_vtx<-read.delim("PT_AMF_VTX_MA_RA_R.txt", header=T,row.names=1, dec=",")
str(PT_vtx)
head(PT_vtx)
colnames(PT_vtx)

##Creating the right database by merging split databases:
##a) inoculated corn
vtx_corn<-PT_vtx[PT_vtx$Crop=="Corn",]
vtx_corn_inoculated<-vtx_corn[vtx_corn$Inoculation=="yes",]

(vtx_corn_inoculated_graph<-colMeans(vtx_corn_inoculated[,1:46]))
length(vtx_corn_inoculated_graph)
write.table(vtx_corn_inoculated_graph, file="vtx_corn_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##b) non-inoculated corn
vtx_corn_noninoculated<-vtx_corn[vtx_corn$Inoculation=="no",]

(vtx_corn_noninoculated_graph<-colMeans(vtx_corn_noninoculated[,1:46]))
length(vtx_corn_noninoculated_graph)
write.table(vtx_corn_noninoculated_graph, file="vtx_corn_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##c) inoculated wheat
vtx_wheat<-PT_vtx[PT_vtx$Crop=="Wheat",]
vtx_wheat_inoculated<-vtx_wheat[vtx_wheat$Inoculation=="yes",]

(vtx_wheat_inoculated_graph<-colMeans(vtx_wheat_inoculated[,1:46]))
length(vtx_wheat_inoculated_graph)
write.table(vtx_wheat_inoculated_graph, file="vtx_wheat_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##d) non-inoculated wheat
vtx_wheat<-PT_vtx[PT_vtx$Crop=="Wheat",]
vtx_wheat_noninoculated<-vtx_wheat[vtx_wheat$Inoculation=="no",]

(vtx_wheat_noninoculated_graph<-colMeans(vtx_wheat_noninoculated[,1:46]))
length(vtx_wheat_noninoculated_graph)
write.table(vtx_wheat_noninoculated_graph, file="vtx_wheat_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##e) inoculated soy
vtx_soy<-PT_vtx[PT_vtx$Crop=="Soy",]
vtx_soy_inoculated<-vtx_soy[vtx_soy$Inoculation=="yes",]

(vtx_soy_inoculated_graph<-colMeans(vtx_soy_inoculated[,1:46]))
length(vtx_soy_inoculated_graph)
write.table(vtx_soy_inoculated_graph, file="vtx_soy_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##f) non-inoculated soy
vtx_soy<-PT_vtx[PT_vtx$Crop=="Soy",]
vtx_soy_noninoculated<-vtx_soy[vtx_soy$Inoculation=="no",]

(vtx_soy_noninoculated_graph<-colMeans(vtx_soy_noninoculated[,1:46]))
length(vtx_soy_noninoculated_graph)
write.table(vtx_soy_noninoculated_graph, file="vtx_soy_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")


###So I merged those databases in Excel and created a new one: 
VTX_barplot<-read.delim("VTX_barplot.txt", header=T, dec=",")
str(VTX_barplot)

##Graph 1 with only the VTX numbers:
##Adjusting palette for more than 12 colors
library("RColorBrewer")
colourCount<-length(unique(VTX_barplot$VTX_all))
getPalette<-colorRampPalette(brewer.pal(9, "Set1"))
getPalette

pvtx_fill_1<-ggplot(VTX_barplot, aes(x=catergory, y=relativeAbundance, fill=just_VTX)) +
  geom_bar(stat="identity", width=0.5, position="fill") +
  scale_fill_manual(values = getPalette(colourCount))+ 
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white', colour="gray"), axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +
  coord_cartesian(ylim=c(0,1)) +
  guides(fill=guide_legend(title="Virtual taxa", ncol=12)) +
  ylab("Average relative abundance (%)") 
pvtx_fill_1

##Graph 2 with the full names of VTX:
##Adjusting palette for more than 12 colors
library("RColorBrewer")
colourCount<-length(unique(VTX_barplot$VTX_all))
getPalette<-colorRampPalette(brewer.pal(9, "Set1"))
getPalette

pvtx_fill_2<-ggplot(VTX_barplot, aes(x=catergory, y=relativeAbundance, fill=VTX_all)) +
  geom_bar(stat="identity", width=0.5, position="fill") +
  scale_fill_manual(values = getPalette(colourCount))+ 
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white', colour="gray"), axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +
  coord_cartesian(ylim=c(0,1)) +
  guides(fill=guide_legend(title=NULL, ncol=4)) +
  ylab("Average relative abundance (%)") 
pvtx_fill_2




###Per orders
##Creating the right database by merging split databases:
PT_orders<-read.delim("PT_AMF_orders_MA_RA_R.txt", row.names=1, header=T, dec=",")
##a) inoculated corn
orders_corn<-PT_orders[PT_orders$Crop=="Corn",]
orders_corn_inoculated<-orders_corn[orders_corn$Inoculation=="yes",]

(orders_corn_inoculated_graph<-colMeans(orders_corn_inoculated[,1:3]))
length(orders_corn_inoculated_graph)
write.table(orders_corn_inoculated_graph, file="orders_corn_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##b) non-inoculated corn
orders_corn_noninoculated<-orders_corn[orders_corn$Inoculation=="no",]

(orders_corn_noninoculated_graph<-colMeans(orders_corn_noninoculated[,1:3]))
length(orders_corn_noninoculated_graph)
write.table(orders_corn_noninoculated_graph, file="orders_corn_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##c) inoculated wheat
orders_wheat<-PT_orders[PT_orders$Crop=="Wheat",]
orders_wheat_inoculated<-orders_wheat[orders_wheat$Inoculation=="yes",]

(orders_wheat_inoculated_graph<-colMeans(orders_wheat_inoculated[,1:3]))
length(orders_wheat_inoculated_graph)
write.table(orders_wheat_inoculated_graph, file="orders_wheat_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##d) non-inoculated wheat
orders_wheat<-PT_orders[PT_orders$Crop=="Wheat",]
orders_wheat_noninoculated<-orders_wheat[orders_wheat$Inoculation=="no",]

(orders_wheat_noninoculated_graph<-colMeans(orders_wheat_noninoculated[,1:3]))
length(orders_wheat_noninoculated_graph)
write.table(orders_wheat_noninoculated_graph, file="orders_wheat_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##e) inoculated soy
orders_soy<-PT_orders[PT_orders$Crop=="Soy",]
orders_soy_inoculated<-orders_soy[orders_soy$Inoculation=="yes",]

(orders_soy_inoculated_graph<-colMeans(orders_soy_inoculated[,1:3]))
length(orders_soy_inoculated_graph)
write.table(orders_soy_inoculated_graph, file="orders_soy_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##f) non-inoculated soy
orders_soy<-PT_orders[PT_orders$Crop=="Soy",]
orders_soy_noninoculated<-orders_soy[orders_soy$Inoculation=="no",]

(orders_soy_noninoculated_graph<-colMeans(orders_soy_noninoculated[,1:3]))
length(orders_soy_noninoculated_graph)
write.table(orders_soy_noninoculated_graph, file="orders_soy_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

###So I merged those databases in Excel and created a new one: 
orders_barplot<-read.delim("Orders_barpot.txt", header=T, dec=",")
str(orders_barplot)


##Graph:
##Adjusting palette for more than 12 colors
library("RColorBrewer")
colourCount<-length(unique(orders_barplot$orders))
getPalette<-colorRampPalette(brewer.pal(9, "Set1"))
getPalette

porders_fill_1<-ggplot(orders_barplot, aes(x=category, y=relative_abundance, fill=orders)) +
  geom_bar(stat="identity", width=0.5, position="fill") +
  scale_fill_manual(values = getPalette(colourCount))+ 
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white', colour="gray"), axis.title.x=element_blank()) + 
  theme(legend.position="right") +
  coord_cartesian(ylim=c(0,1)) +
  guides(fill=guide_legend(title="Orders", ncol=1)) +
  ylab("Average relative abundance (%)") 
porders_fill_1


###Per Families
PT_families<-read.delim("PT_AMF_families_MA_RA.txt", row.names=1, header=T, dec=",")
str(PT_families)

##Creating the right database by merging split databases:
##a) inoculated corn
families_corn<-PT_families[PT_families$Crop=="Corn",]
families_corn_inoculated<-families_corn[families_corn$Inoculation=="yes",]

(families_corn_inoculated_graph<-colMeans(families_corn_inoculated[,1:7]))
length(families_corn_inoculated_graph)
write.table(families_corn_inoculated_graph, file="families_corn_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##b) non-inoculated corn
families_corn_noninoculated<-families_corn[families_corn$Inoculation=="no",]

(families_corn_noninoculated_graph<-colMeans(families_corn_noninoculated[,1:7]))
length(families_corn_noninoculated_graph)
write.table(families_corn_noninoculated_graph, file="families_corn_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##c) inoculated wheat
families_wheat<-PT_families[PT_families$Crop=="Wheat",]
families_wheat_inoculated<-families_wheat[families_wheat$Inoculation=="yes",]

(families_wheat_inoculated_graph<-colMeans(families_wheat_inoculated[,1:7]))
length(families_wheat_inoculated_graph)
write.table(families_wheat_inoculated_graph, file="families_wheat_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##d) non-inoculated wheat
families_wheat<-PT_families[PT_families$Crop=="Wheat",]
families_wheat_noninoculated<-families_wheat[families_wheat$Inoculation=="no",]

(families_wheat_noninoculated_graph<-colMeans(families_wheat_noninoculated[,1:7]))
length(families_wheat_noninoculated_graph)
write.table(families_wheat_noninoculated_graph, file="families_wheat_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##e) inoculated soy
families_soy<-PT_families[PT_families$Crop=="Soy",]
families_soy_inoculated<-families_soy[families_soy$Inoculation=="yes",]

(families_soy_inoculated_graph<-colMeans(families_soy_inoculated[,1:7]))
length(families_soy_inoculated_graph)
write.table(families_soy_inoculated_graph, file="families_soy_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##f) non-inoculated soy
families_soy<-PT_families[PT_families$Crop=="Soy",]
families_soy_noninoculated<-families_soy[families_soy$Inoculation=="no",]

(families_soy_noninoculated_graph<-colMeans(families_soy_noninoculated[,1:7]))
length(families_soy_noninoculated_graph)
write.table(families_soy_noninoculated_graph, file="families_soy_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

###So I merged those databases in Excel and created a new one: 
families_barplot<-read.delim("Families_barplot.txt", header=T, dec=",")
str(families_barplot)

##Adjusting palette for more than 12 colors
library("RColorBrewer")
colourCount<-length(unique(families_barplot$families))
getPalette<-colorRampPalette(brewer.pal(9, "Set1"))
getPalette

pfamilies_fill_1<-ggplot(families_barplot, aes(x=category, y=relative_abundance, fill=families)) +
  geom_bar(stat="identity", width=0.5, position="fill") +
  scale_fill_manual(values = getPalette(colourCount)) + 
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white', colour="gray"), axis.title.x=element_blank()) + 
  theme(legend.position="right") +
  coord_cartesian(ylim=c(0,1)) +
  guides(fill=guide_legend(title="Orders", ncol=1)) +
  ylab("Average relative abundance (%)") 
pfamilies_fill_1

###Per Genera
PT_genera<-read.delim("PT_AMF_genera_MA_RA_R.txt", row.names=1, header=T, dec=",")
str(PT_genera)

##Creating the right database by merging split databases:
##a) inoculated corn
genera_corn<-PT_genera[PT_genera$Crop=="Corn",]
genera_corn_inoculated<-genera_corn[genera_corn$Inoculation=="yes",]

(genera_corn_inoculated_graph<-colMeans(genera_corn_inoculated[,1:7]))
length(genera_corn_inoculated_graph)
write.table(genera_corn_inoculated_graph, file="genera_corn_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##b) non-inoculated corn
genera_corn_noninoculated<-genera_corn[genera_corn$Inoculation=="no",]

(genera_corn_noninoculated_graph<-colMeans(genera_corn_noninoculated[,1:7]))
length(genera_corn_noninoculated_graph)
write.table(genera_corn_noninoculated_graph, file="genera_corn_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##c) inoculated wheat
genera_wheat<-PT_genera[PT_genera$Crop=="Wheat",]
genera_wheat_inoculated<-genera_wheat[genera_wheat$Inoculation=="yes",]

(genera_wheat_inoculated_graph<-colMeans(genera_wheat_inoculated[,1:7]))
length(genera_wheat_inoculated_graph)
write.table(genera_wheat_inoculated_graph, file="genera_wheat_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##d) non-inoculated wheat
genera_wheat<-PT_genera[PT_genera$Crop=="Wheat",]
genera_wheat_noninoculated<-genera_wheat[genera_wheat$Inoculation=="no",]

(genera_wheat_noninoculated_graph<-colMeans(genera_wheat_noninoculated[,1:7]))
length(genera_wheat_noninoculated_graph)
write.table(genera_wheat_noninoculated_graph, file="genera_wheat_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##e) inoculated soy
genera_soy<-PT_genera[PT_genera$Crop=="Soy",]
genera_soy_inoculated<-genera_soy[genera_soy$Inoculation=="yes",]

(genera_soy_inoculated_graph<-colMeans(genera_soy_inoculated[,1:7]))
length(genera_soy_inoculated_graph)
write.table(genera_soy_inoculated_graph, file="genera_soy_inoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

##f) non-inoculated soy
genera_soy<-PT_genera[PT_genera$Crop=="Soy",]
genera_soy_noninoculated<-genera_soy[genera_soy$Inoculation=="no",]

(genera_soy_noninoculated_graph<-colMeans(genera_soy_noninoculated[,1:7]))
length(genera_soy_noninoculated_graph)
write.table(genera_soy_noninoculated_graph, file="genera_soy_noninoculated_graph", row.names=TRUE, col.names = TRUE, sep = "\t")

###So I merged those databases in Excel and created a new one: 
genera_barplot<-read.delim("genera_barplot.txt", header=T, dec=",")
str(families_barplot)

##Adjusting palette for more than 12 colors
library("RColorBrewer")
colourCount<-length(unique(genera_barplot$genera))
getPalette<-colorRampPalette(brewer.pal(6, "Set1"))
getPalette

pgenera_fill_1<-ggplot(genera_barplot, aes(x=category, y=relative_abundance, fill=genera)) +
  geom_bar(stat="identity", width=0.5, position="fill") +
  scale_fill_manual(values = getPalette(colourCount)) + 
  theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white', colour="gray"), axis.title.x=element_blank()) + 
  theme(legend.position="right") +
  coord_cartesian(ylim=c(0,1)) +
  guides(fill=guide_legend(title="Genera", ncol=1)) +
  ylab("Average relative abundance (%)") 
pgenera_fill_1


# Pie-chart - communities -------------------------------------------------

pvtx_pie<-ggplot(VTX_barplot, aes(factor(1), relativeAbundance, fill=just_VTX)) +
  geom_bar(stat="identity", position="fill", width=1) +
  scale_fill_manual(values = getPalette(colourCount)) +
  coord_polar(theta="y") + 
  facet_grid(facets=. ~ catergory) +
  #geom_errorbar(aes(ymin=value-sd1, ymax=value+sd1), size=0.7, width=.15, colour="black") + 
  theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Virtual taxa", ncol=12)) +
  theme(axis.title.y=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', colour="white"),
        panel.grid  = element_blank()) + 
  ylab("Average relative abundance (%)") 
pvtx_pie

pvtx_pie2<-ggplot(VTX_barplot, aes(factor(1), relativeAbundance, fill=VTX_all)) +
  geom_bar(stat="identity", position="fill", width=1) +
  scale_fill_manual(values = getPalette(colourCount)) +
  coord_polar(theta="y") + 
  facet_grid(facets=. ~ catergory) +
  #geom_errorbar(aes(ymin=value-sd1, ymax=value+sd1), size=0.7, width=.15, colour="black") + 
  theme(legend.position="bottom") +
  guides(fill=guide_legend(title=NULL, ncol=4)) +
  theme(axis.title.y=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', colour="white"),
        panel.grid  = element_blank()) + 
  ylab("Average relative abundance (%)") 
pvtx_pie2


porders_pie<-ggplot(orders_barplot, aes(factor(1), relative_abundance, fill=orders)) +
  geom_bar(stat="identity", position="fill", width=1) +
  scale_fill_manual(values = getPalette(colourCount)) +
  coord_polar(theta="y") + 
  facet_grid(facets=. ~ category) +
  #geom_errorbar(aes(ymin=value-sd1, ymax=value+sd1), size=0.7, width=.15, colour="black") + 
  theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Orders", ncol=7)) +
  theme(axis.title.y=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', colour="white"),
        panel.grid  = element_blank()) + 
  ylab("Average relative abundance (%)") 
porders_pie

pfamilies_pie<-ggplot(families_barplot, aes(factor(1), relative_abundance, fill=families)) +
  geom_bar(stat="identity", position="fill", width=1) +
  scale_fill_manual(values = getPalette(colourCount)) +
  coord_polar(theta="y") + 
  facet_grid(facets=. ~ category) +
  #geom_errorbar(aes(ymin=value-sd1, ymax=value+sd1), size=0.7, width=.15, colour="black") + 
  theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Families", ncol=4)) +
  theme(axis.title.y=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', colour="white"),
        panel.grid  = element_blank()) + 
  ylab("Average relative abundance (%)") 
pfamilies_pie


pgenera_pie<-ggplot(genera_barplot, aes(factor(1), relative_abundance, fill=genera)) +
  geom_bar(stat="identity", position="fill", width=1) +
  scale_fill_manual(values = getPalette(colourCount)) +
  coord_polar(theta="y") + 
  facet_grid(facets=. ~ category) +
  #geom_errorbar(aes(ymin=value-sd1, ymax=value+sd1), size=0.7, width=.15, colour="black") + 
  theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Genera", ncol=4)) +
  theme(axis.title.y=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', colour="white"),
        panel.grid  = element_blank()) + 
  ylab("Average relative abundance (%)") 
pgenera_pie




# Alpha-diversity - ANOVA - general model---------------------------------------------------
##Preparation of the data
alphaD<-read.delim("PT_AMF_clustered_alphaD.txt", dec=".", row.names=1, header=T)
str(alphaD)
dim(alphaD)

##All crops together testing inoculation
boxplot(invsimpson ~ Inoculation, data=alphaD, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"))

plot.design(invsimpson ~ Inoculation, data=alphaD)
##Probably not different, but the plot.design graph showed that plots that were not inoculated seems to have a higher alpha_diversity

boxplot(invsimpson ~ crop, data=alphaD, col=c("darkseagreen3","gold2", "royalblue"), xlab="Crop", ylab=("Inverse Simpson Index"))

plot.design(invsimpson ~ crop, data=alphaD)
#but the plot.design seems to show that soy have a higher alphaD

boxplot(invsimpson ~ Growth_stage, data=alphaD, col=c("orangered", "olivedrab1"), outline=F, xlab="Growth_stage", ylab=("Inverse Simpson Index"))

###ANOVA - general model with interactions
anov_all_interactions<-lm(invsimpson ~ crop*Inoculation*Growth_stage, data=alphaD)
anova(anov_all_interactions)
##Only crop is significantly different

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_all_interactions)

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_all_interactions))
##W = 0.94091, p-value = 6.599e-05
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(invsimpson ~ crop, data=alphaD)
#Bartlett's K-squared = 15.5, df = 2, p-value = 0.0004307 - variances not equal
bartlett.test(invsimpson ~ Inoculation, data=alphaD)
#Bartlett's K-squared = 0.7673, df = 1, p-value = 0.3811 - variances equal
bartlett.test(invsimpson ~ Growth_stage, data=alphaD)
#Bartlett's K-squared = 17.62, df = 1, p-value = 2.698e-05
##the variance are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##Transforming data because assumption 1 and 2 are not respected
alphaD$sqrt_invsimpson<-sqrt(alphaD$invsimpson)
names(alphaD)

alphaD$log_invsimpson<-log(alphaD$invsimpson)
names(alphaD)

##Do a second anova
anov_all_interactions2<-lm(sqrt_invsimpson ~ crop*Inoculation*Growth_stage, data=alphaD)
anova(anov_all_interactions2)
##Not significantly different

par(mfrow=c(2,2))
plot(anov_all_interactions2)
##a bit better

shapiro.test(resid(anov_all_interactions2))
#W = 0.98576, p-value = 0.261
##Normally distributed
bartlett.test(sqrt_invsimpson ~ crop, data=alphaD)
#Bartlett's K-squared = 5.069, df = 2, p-value = 0.0793 - variances are equal
bartlett.test(sqrt_invsimpson ~ Inoculation, data=alphaD)
#Bartlett's K-squared = 0.48669, df = 1, p-value = 0.4854 - variances are equal
bartlett.test(sqrt_invsimpson ~ Growth_stage, data=alphaD)
#BBartlett's K-squared = 11.084, df = 1, p-value = 0.0008707
##the variance are not equal

##Doing ANova on log transformed data
anov_all_interaction3<-lm(log_invsimpson ~ crop*Inoculation*Growth_stage, data=alphaD)
anova(anov_all_interaction3)  
###Still not different

par(mfrow=c(2,2))
plot(anov_all_interaction3)

shapiro.test(resid(anov_all_interaction3))  
###W = 0.99456, p-value = 0.9345   ##normally distributed
bartlett.test(log_invsimpson ~ crop, data=alphaD)
#Bartlett's K-squared = 0.78239, df = 2, p-value = 0.6762 - variances are equal
bartlett.test(log_invsimpson ~ Inoculation, data=alphaD)
#artlett's K-squared = 0.19918, df = 1, p-value = 0.6554 - variances are equal
bartlett.test(log_invsimpson ~ Growth_stage, data=alphaD)
#Bartlett's K-squared = 7.2824, df = 1, p-value = 0.006963
##the variance are not equal, but I don't think they could ever be 

##So we will take the log as the good model:

anova(anov_all_interaction3)  
#Response: log_invsimpson
#                               Df  Sum Sq Mean Sq F value   Pr(>F)    
#crop                            2  6.5848  3.2924 14.6997 2.32e-06 ***
#Inoculation                     1  0.0000  0.0000  0.0001   0.9923    
#Growth_stage                    1  0.5629  0.5629  2.5132   0.1159    
#crop:Inoculation                2  0.5088  0.2544  1.1359   0.3250    
#crop:Growth_stage               1  0.2662  0.2662  1.1884   0.2781    
#Inoculation:Growth_stage        1  0.0814  0.0814  0.3632   0.5480    
#crop:Inoculation:Growth_stage   1  0.5440  0.5440  2.4289   0.1221    
#Residuals                     106 23.7415  0.2240                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#ANOVA on crop (log inv)
##Transforming data because assumption 1 and 2 are not respected
alphaD$sqrt_invsimpson<-sqrt(alphaD$invsimpson)
names(alphaD)

alphaD$log_invsimpson<-log(alphaD$invsimpson)
names(alphaD)

##Doing ANova on log transformed data
anov_all_crop2<-lm(log_invsimpson ~ crop, data=alphaD)
anova(anov_all_crop2)  
###Still not different

par(mfrow=c(2,2))
plot(anov_all_crop2)
##Normality of residus seems better, but not homeostacity

shapiro.test(resid(anov_all_crop2))  ###Great! The residus are normally distributed
bartlett.test(log_invsimpson ~ crop, data=alphaD)  ##Variance are now homogeneous
##With the log transformation all the assumptions are met: so the final result is: 

anov_all_crop2<-lm(log_invsimpson ~ crop, data=alphaD)
anova(anov_all_crop2) 
##The alpha diversity among the different crop is different at F(2,113) = 14,474 and p=2.53e-06

##Post-hoc test
TukeyHSD(aov(anov_all_crop2),ordered=T)

#               diff         lwr       upr     p adj
#Wheat-Corn 0.1582979 -0.08542903 0.4020247 0.2751011
#Soy-Corn   0.6106707  0.34028002 0.8810614 0.0000013
#Soy-Wheat  0.4523729  0.16083952 0.7439062 0.0010204

##So soy is different from the other two crops having a higher alphaD


# alpha-diversity - anova - per crop --------------------------------------
alphaD_Corn<-alphaD[alphaD$crop=="Corn",]
head(alphaD_Corn)
dim(alphaD_Corn)

alphaD_corn_early<-alphaD_Corn[alphaD_Corn$Growth_stage=="Early",]
dim(alphaD_corn_early)

alphaD_corn_late<-alphaD_Corn[alphaD_Corn$Growth_stage=="Late",]
dim(alphaD_corn_late)

alphaD_Soy<-alphaD[alphaD$crop=="Soy",]
head(alphaD_Soy)
dim(alphaD_Soy)

alphaD_Wheat<-alphaD[alphaD$crop=="Wheat",]
head(alphaD_Wheat)
dim(alphaD_Wheat)

alphaD_Wheat_early<-alphaD_Wheat[alphaD_Wheat$Growth_stage=="Early",]
dim(alphaD_Wheat_early)

alphaD_Wheat_late<-alphaD_Wheat[alphaD_Wheat$Growth_stage=="Late",]
dim(alphaD_Wheat_late)



#######ANOVA only for Corn
##only Corn
boxplot(invsimpson ~ Inoculation, data=alphaD_Corn, col=c("skyblue1","skyblue4"), xlab="Inoculation", outline=F, ylab=("Inverse Simpson Index"), main="Corn")

boxplot(invsimpson ~ Growth_stage, data=alphaD_Corn, col=c("orangered", "olivedrab1"), outline=F, xlab="Growth stage", ylab=("Inverse Simpson Index"), main="Corn")

anov_corn_interaction<-lm(invsimpson ~ Inoculation*Growth_stage, data=alphaD_Corn)
anova(anov_corn_interaction)

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_corn_interaction)

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_corn_interaction))
##W = 0.91333, p-value = 0.0008431
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(invsimpson ~ Inoculation, data=alphaD_Corn)
#Bartlett's K-squared = 0.00048911, df = 1, p-value = 0.9824  - the variance are equal
bartlett.test(invsimpson ~ Growth_stage, data=alphaD_Corn)
#Bartlett's K-squared = 1.9233, df = 1, p-value = 0.1655  - the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##Transforming data because assumption 1 and 2 are not respected
alphaD_Corn$log_invsimpson<-log(alphaD_Corn$invsimpson)
names(alphaD_Corn)

alphaD_Corn$sqrt_invsimpson<-sqrt(alphaD_Corn$invsimpson)
names(alphaD_Corn)

##Do a second anova
anov_corn_interaction2<-lm(log_invsimpson ~ Inoculation*Growth_stage, data=alphaD_Corn)
anova(anov_corn_interaction2)
##Growth stage is significantly different

par(mfrow=c(2,2))
plot(anov_corn_interaction2)
##a bit better

shapiro.test(resid(anov_corn_interaction2))
#W = 0.97809, p-value = 0.4229 - Normally distributed
bartlett.test(log_invsimpson ~ Inoculation, data=alphaD_Corn) 
#Bartlett's K-squared = 0.047554, df = 1, p-value = 0.8274 -  Equality of variances
bartlett.test(log_invsimpson ~ Growth_stage, data=alphaD_Corn)
#Bartlett's K-squared = 1.4099, df = 1, p-value = 0.2351  - Equality of variances

##With the log transformation all the assumptions are met: so the final result is: 
anova(anov_corn_interaction2)

#Response: log_invsimpson
#                         Df Sum Sq Mean Sq F value  Pr(>F)  
#Inoculation               1 0.0910 0.09100  0.4940 0.48539  
#Growth_stage              1 0.8480 0.84802  4.6042 0.03677 *
#Inoculation:Growth_stage  1 0.4832 0.48316  2.6232 0.11160  
#Residuals                50 9.2092 0.18418                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

boxplot(invsimpson ~ Growth_stage, data=alphaD_Corn, col=c("orangered", "olivedrab1"), outline=F, xlab="Growth_stage", ylab=("Inverse Simpson Index"))


###ANOVA only for Wheat
##Only Wheat
boxplot(invsimpson ~ Inoculation, data=alphaD_Wheat, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"))

boxplot(invsimpson ~ Growth_stage, data=alphaD_Wheat, col=c("orangered", "olivedrab1"), outline=F, xlab="Growth_stage", ylab=("Inverse Simpson Index"))

##ANOVA with interactions
anov_wheat_interaction<-lm(invsimpson ~ Inoculation*Growth_stage, data=alphaD_Wheat)
anova(anov_wheat_interaction)

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_wheat_interaction)

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_wheat_interaction))
##W = 0.93149, p-value = 0.02777
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(invsimpson ~ Inoculation, data=alphaD_Wheat)
#Bartlett's K-squared = 0.049772, df = 1, p-value = 0.8235  - the variance are equal
bartlett.test(invsimpson ~ Growth_stage, data=alphaD_Corn)
#Bartlett's K-squared = 1.9233, df = 1, p-value = 0.1655  - the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##Transforming data because assumption 1 and 2 are not respected
alphaD_Wheatlog_invsimpson<-log(alphaD_Wheat$invsimpson)
names(alphaD_Wheat)

alphaD_Wheat$sqrt_invsimpson<-sqrt(alphaD_Wheat$invsimpson)
names(alphaD_Wheat)

##Do a second anova
anov_wheat_interaction2<-lm(log_invsimpson ~ Inoculation*Growth_stage, data=alphaD_Wheat)
anova(anov_wheat_interaction2)
##Nothing is different

par(mfrow=c(2,2))
plot(anov_wheat_interaction2)
##a bit better

shapiro.test(resid(anov_wheat_interaction2))
#W = 0.97286, p-value = 0.5088 - Normally distributed
bartlett.test(log_invsimpson ~ Inoculation, data=alphaD_Wheat) 
#Bartlett's K-squared = 0.15215, df = 1, p-value = 0.6965 -  Equality of variances
bartlett.test(log_invsimpson ~ Growth_stage, data=alphaD_Wheat)
#Bartlett's K-squared = 0.40722, df = 1, p-value = 0.5234 - Equality of variances

##With the log transformation all the assumptions are met: so the final result is: 
anova(anov_wheat_interaction2)

#Response: log_invsimpson
#                         Df Sum Sq  Mean Sq F value Pr(>F)
#Inoculation               1 0.0279 0.027874  0.1085 0.7440
#Growth_stage              1 0.0062 0.006234  0.0243 0.8772
#Inoculation:Growth_stage  1 0.1422 0.142217  0.5537 0.4623
#Residuals                32 8.2197 0.256866               



###ANOVA only for Soy

##Only Soy
boxplot(invsimpson ~ Inoculation, data=alphaD_Soy, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"))


anov_soy_inoculation<-lm(invsimpson ~ Inoculation, data=alphaD_Soy)
anova(anov_soy_inoculation)
##Not significantly different

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_soy_inoculation)
##Doesn't look good: Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: probably because the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_soy_inoculation))
##W = 0.92675, p-value = 0.06485
##Residus are normally distributed

##Assumption 2: Homoegeneity of the variances
bartlett.test(invsimpson ~ Inoculation, data=alphaD_Soy)
#Bartlett's K-squared = 0.045629, df = 1, p-value = 0.8309
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 
##So all the assumptions are met, the final anova is:  

anov_soy_inoculation<-lm(invsimpson ~ Inoculation, data=alphaD_Soy)
anova(anov_soy_inoculation)
#Response: invsimpson
#             Df  Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1   7.525  7.5253  1.0479 0.3162
#Residuals   24 172.356  7.1815

###The alpha-diversity is not significantly different between inoculated and none-inoculated in the soy fields at F(1,24)=1.0479 and p=0.3162




# alpha-diversity - anova - per crop and growth stage --------------------------------------
##Early-corn
alphaD_corn_early<-alphaD_Corn[alphaD_Corn$Growth_stage=="Early",]
dim(alphaD_corn_early)

boxplot(invsimpson ~ Inoculation, data=alphaD_corn_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"), outline=F)


anov_corn_early<-lm(invsimpson ~ Inoculation, data=alphaD_corn_early)
anova(anov_corn_early)

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_corn_early)

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_corn_early))
##W = 0.81984, p-value = 0.0003905
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(invsimpson ~ Inoculation, data=alphaD_corn_early)
#Bartlett's K-squared = 4.8036, df = 1, p-value = 0.0284  - the variance are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##Transforming data because assumption 1 and 2 are not respected
alphaD_corn_early$log_invsimpson<-log(alphaD_corn_early$invsimpson)
names(alphaD_corn_early)

alphaD_corn_early$sqrt_invsimpson<-sqrt(alphaD_corn_early$invsimpson)
names(alphaD_corn_early)

##Do a second anova
anov_corn_early2<-lm(log_invsimpson ~ Inoculation, data=alphaD_corn_early)
anova(anov_corn_early2)
##Nothiring is different

par(mfrow=c(2,2))
plot(anov_corn_early2)
##a bit better

shapiro.test(resid(anov_corn_early2))
#W = 0.94833, p-value = 0.2119 - Normally distributed
bartlett.test(log_invsimpson ~ Inoculation, data=alphaD_corn_early) 
#Bartlett's K-squared = 1.2018, df = 1, p-value = 0.273 -  Equality of variances

##With the log transformation all the assumptions are met: so the final result is: 
anova(anov_corn_early2)
#Response: log_invsimpson
#               Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1 0.4823 0.48231  2.1613 0.1545
#Residuals   24 5.3558 0.22316 

##So the alphaD is not different between inoculated and non-inoculated corn fields at the early stage growth

####Corn - late
alphaD_corn_late<-alphaD_Corn[alphaD_Corn$Growth_stage=="Late",]
dim(alphaD_corn_late)

par(mfrow=c(1,1))
boxplot(invsimpson ~ Inoculation, data=alphaD_corn_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"), outline=F)

anov_corn_late<-lm(invsimpson ~ Inoculation, data=alphaD_corn_late)
anova(anov_corn_late)

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_corn_late)

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_corn_late))
##W = 0.95137, p-value = 0.2144
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(invsimpson ~ Inoculation, data=alphaD_corn_late)
#Bartlett's K-squared = 0.43898, df = 1, p-value = 0.5076  - the variance are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##all the conditions are met so we can analyse the first model: 
 
anova(anov_corn_late)
#Response: invsimpson
#             Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1  0.926  0.9256  0.6411 0.4306
#Residuals   26 37.536  1.4437

##So the alphaD is not different between inoculated and non-inoculated corn fields at the late  stage growth

###Early - wheat
alphaD_Wheat_early<-alphaD_Wheat[alphaD_Wheat$Growth_stage=="Early",]
dim(alphaD_Wheat_early)

boxplot(invsimpson ~ Inoculation, data=alphaD_Wheat_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"), outline=F)

anov_wheat_early<-lm(invsimpson ~ Inoculation, data=alphaD_Wheat_early)
anova(anov_wheat_early)

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_wheat_early)

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_wheat_early))
##W = 0.88348, p-value = 0.0247
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(invsimpson ~ Inoculation, data=alphaD_Wheat_early)
#Bartlett's K-squared = 0.2848, df = 1, p-value = 0.5936  - the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##Transforming data because assumption 1 and 2 are not respected
alphaD_Wheat_early$log_invsimpson<-log(alphaD_Wheat_early$invsimpson)
names(alphaD_Wheat_early)

alphaD_Wheat_early$sqrt_invsimpson<-sqrt(alphaD_Wheat_early$invsimpson)
names(alphaD_Wheat_early)

##Do a second anova
anov_wheat_early2<-lm(log_invsimpson ~ Inoculation, data=alphaD_Wheat_early)
anova(anov_wheat_early2)
##Nothiring is different

par(mfrow=c(2,2))
plot(anov_wheat_early2)
##a bit better

shapiro.test(resid(anov_wheat_early2))
#W = 0.9438, p-value = 0.3083- Normally distributed
bartlett.test(log_invsimpson ~ Inoculation, data=alphaD_Wheat_early) 
#Bartlett's K-squared = 0.10488, df = 1, p-value = 0.7461-  Equality of variances

##With the log transformation all the assumptions are met: so the final result is: 
anova(anov_wheat_early2)
#Response: log_invsimpson
#             Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1 0.0300 0.02999  0.1007 0.7549
#Residuals   17 5.0635 0.29785 

##So the alphaD is not different between inoculated and non-inoculated corn fields at the early stage growth


###Wheat - late
alphaD_Wheat_late<-alphaD_Wheat[alphaD_Wheat$Growth_stage=="Late",]
dim(alphaD_Wheat_late)

par(mfrow=c(1,1))
boxplot(invsimpson ~ Inoculation, data=alphaD_Wheat_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"), outline=F)

anov_wheat_late<-lm(invsimpson ~ Inoculation, data=alphaD_Wheat_late)
anova(anov_wheat_late)

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_wheat_late)

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_wheat_late))
##W = 0.97107, p-value = 0.8368
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(invsimpson ~ Inoculation, data=alphaD_Wheat_late)
#Bartlett's K-squared = 0.94926, df = 1, p-value = 0.3299 - the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

## all the assumptions are met, so we can interpret the result obtained by the first anova 
anova(anov_wheat_late)
#Response: invsimpson
#             Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1  2.612  2.6122  1.0428 0.3234
#Residuals   15 37.575  2.5050 

##So the alphaD is not different between inoculated and non-inoculated wheat fields at the late stage growth


# Alpha-diversity - graph boxplot -----------------------------------------

###Graph - inoculation general
boxplot(invsimpson ~ Inoculation, data=alphaD, col=c("skyblue1","skyblue4"), outline=F, xlab="Inoculation", ylab=("Inverse Simpson Index"))

##Graph - crop general
boxplot(invsimpson ~ crop, data=alphaD, col=c("darkseagreen3","gold2", "royalblue"), outline=F, xlab="Crop", ylab=("Inverse Simpson Index"))

###Graph - growing stage general
boxplot(invsimpson ~ Growth_stage, data=alphaD, col=c("orangered", "olivedrab1"), outline=F, xlab="Growth stage", ylab=("Inverse Simpson Index"))

##Graph - inoculation
##only Corn
boxplot(invsimpson ~ Inoculation, data=alphaD_Corn, col=c("skyblue1","skyblue4"), outline=F, xlab="Inoculation", ylab="Inverse Simpson Index", main="Corn")

boxplot(invsimpson ~ Growth_stage, data=alphaD_Corn, col=c("orangered", "olivedrab1"), outline=F, xlab="Growth stage", ylab=("Inverse Simpson Index"), main="Corn")

##Only wheat
boxplot(invsimpson ~ Inoculation, data=alphaD_Wheat, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"),outline=F, main="Wheat")

boxplot(invsimpson ~ Growth_stage, data=alphaD_Wheat, col=c("orangered", "olivedrab1"), outline=F, xlab="Growth_stage", ylab=("Inverse Simpson Index"),  main="Wheat")


##Only Soy
boxplot(invsimpson ~ Inoculation, data=alphaD_Soy, col=c("skyblue1","skyblue4"), outline=F, xlab="Inoculation", ylab="Inverse Simpson Index", main="Soy")

###Early-corn
boxplot(invsimpson ~ Inoculation, data=alphaD_corn_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"), outline=F, main="Early-corn")

##Late-corn
boxplot(invsimpson ~ Inoculation, data=alphaD_corn_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"), outline=F, main="Late-corn")

##Early-wheat
boxplot(invsimpson ~ Inoculation, data=alphaD_Wheat_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"), outline=F, main="Early-wheat")

##Late-wheat
boxplot(invsimpson ~ Inoculation, data=alphaD_Wheat_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Inverse Simpson Index"), outline=F, main="Late-wheat")


# Beta-diversity - OTU raw abundances----------------------------------------------------------
PT_OTU<-read.delim("PT_AMF_clustered_OTU_R.txt", row.names=1, header=T, dec=",")
dim(PT_OTU)
str(PT_OTU$Crop)
str(PT_OTU$OTU3_04026447528258e.32)

PT_OTU_corn<-PT_OTU[PT_OTU$Crop=="Corn",]
dim(PT_OTU_corn)
str(PT_OTU_corn$Crop)
str(PT_OTU_corn$OTU3_04026447528258e.32)

PT_OTU_soy<-PT_OTU[PT_OTU$Crop=="Soy",]
dim(PT_OTU_soy)
str(PT_OTU_soy$Crop)
str(PT_OTU_soy$OTU3_04026447528258e.32)

PT_OTU_wheat<-PT_OTU[PT_OTU$Crop=="Wheat",]
dim(PT_OTU_wheat)

###General with interactions: 
PT_OTU_justotu<-PT_OTU[,1:408]
PT_OTU_justotu_hel<-decostand(PT_OTU_justotu, "hel")

(adonis(formula=PT_OTU_justotu_hel ~ Crop*Inoculation*Growth_stage, data=PT_OTU, permutations=999, method="bray"))

#                               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Crop                            2    9.5867  4.7933  41.703 0.41715  0.001 ***
#Inoculation                     1    0.0733  0.0733   0.637 0.00319  0.766    
#Growth_stage                    1    0.3354  0.3354   2.918 0.01459  0.009 ** 
#Crop:Inoculation                2    0.3622  0.1811   1.576 0.01576  0.081 .  
#Crop:Growth_stage               1    0.2335  0.2335   2.031 0.01016  0.032 *  
#Inoculation:Growth_stage        1    0.1256  0.1256   1.093 0.00546  0.333    
#Crop:Inoculation:Growth_stage   1    0.0811  0.0811   0.705 0.00353  0.689    
#Residuals                     106   12.1836  0.1149         0.53015           
#Total                         115   22.9813                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Crop generate different microbial communities; growth stage generate different microbial community and those two factors interact together

###Crop: 
(adonis(formula=PT_OTU_justotu ~ Crop, data=PT_OTU, permutations=999, method="bray"))
#           Df SumsOfSqs MeanSqs F.Model       R2   Pr(>F)    
#Crop        2    10.762  5.3808  26.948  0.32293  0.001 ***
#Residuals 113    22.563  0.1997          0.67707           
#Total     115    33.324                  1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###General - Inoculation: 
(adonis(formula=PT_OTU_justotu ~ Inoculation, data=PT_OTU, permutations=9999, method="bray"))

###Corn

##Let's test the assumptions:
#1) Normal distribution: 
norm.her(PT_OTU_corn[,1:408], 0.01)
##Cannot do it cause some otu are absent in the corn fields

##Let's see if a log transformation can help
PT_OTU_corn_log<-decostand(PT_OTU_corn[,1:408], "log")
norm.her(PT_OTU_corn_log, 0.01)
##Cannot do it cause some otu are absent in the corn fields. Let's try to see the graph

#2) Linearity: yes

#3) Homogeneity of Variances: 
homeo.table(PT_OTU_corn[,1:408], PT_OTU_corn$Inoculation, 0.05)
##Doesn't work

##Let's do a permanova then!
##First let's devide the database
#WGRF_16S_otu_factors<-WGRF_otu_16S[,1:5]
PT_OTU_corn_justOTU<-PT_OTU_corn[,1:408]
PT_OTU_corn_justOTU_hel<-decostand(PT_OTU_corn_justOTU, "hel")

##Permanova on untransformed data 
(Permanova_OTU_corn_1<-adonis(formula=PT_OTU_corn_justOTU ~ Inoculation, data=PT_OTU_corn, permutations=999, method="bray"))

#             Df  SumsOfSqs MeanSqs   F.Model      R2   Pr(>F)
#Inoculation  1    0.1197   0.11970   0.80368   0.01522  0.541
#Residuals   52    7.7452   0.14894             0.98478       
#Total       53    7.8649                       1.00000    

##Permanova on hellinger-transformed data 
(Permanova_OTU_corn_2<-adonis(formula=PT_OTU_corn_justOTU_hel ~ Inoculation, data=PT_OTU_corn, permutations=999, method="bray"))

##Not really different, so let's keep the untransformed data:

#             Df  SumsOfSqs   MeanSqs   F.Model      R2       Pr(>F)
#Inoculation  1    0.1288   0.12879     1.2191      0.02291  0.266
#Residuals   52    5.4936   0.10565                0.97709       
#Total       53    5.6224                         1.00000  

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated corn fields

###Soy

##Let's test the assumptions:
#1) Normal distribution: Doesn't work

#2) Linearity: yes

#3) Homogeneity of Variances: 
##Doesn't work

##Let's do a permanova then!
##First let's devide the database
PT_OTU_soy_justOTU<-PT_OTU_soy[,1:408]
PT_OTU_soy_justOTU_hel<-decostand(PT_OTU_soy_justOTU, "hel")

##Permanova on untransformed data 
(Permanova_OTU_soy_1<-adonis(formula=PT_OTU_soy_justOTU ~ Inoculation, data=PT_OTU_soy, permutations=99999, method="bray"))

#             Df  SumsOfSqs MeanSqs   F.Model      R2   Pr(>F)  
#Inoculation  1    0.3456   0.34561  1.6051     0.06269 0.09854 .
#Residuals   24    5.1677   0.21532             0.93731          
#Total       25    5.5133                       1.00000          
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 


##Permanova on hellinger-transformed data 
(Permanova_OTU_soy_2<-adonis(formula=PT_OTU_soy_justOTU_hel ~ Inoculation, data=PT_OTU_soy, permutations=99999, method="bray"))

#             Df    SumsOfSqs   MeanSqs   F.Model      R2     Pr(>F)
#Inoculation  1      0.1860     0.18605   1.3775      0.05428 0.1528
#Residuals   24      3.2415     0.13506                0.94572       
#Total       25      3.4275                           1.00000 

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated soy fields


###Wheat

##Let's test the assumptions:
#1) Normal distribution: Doesn't work

#2) Linearity: yes

#3) Homogeneity of Variances: 
##Doesn't work

##Let's do a permanova then!
##First let's devide the database
PT_OTU_wheat_justOTU<-PT_OTU_wheat[,1:408]
PT_OTU_wheat_justOTU_hel<-decostand(PT_OTU_wheat_justOTU, "hel")

##Permanova on untransformed data 
(Permanova_OTU_wheat_1<-adonis(formula=PT_OTU_wheat_justOTU ~ Inoculation, data=PT_OTU_wheat, permutations=99999, method="bray"))

#             Df  SumsOfSqs MeanSqs   F.Model      R2     Pr(>F)
#Inoculation  1    0.2194   0.21942   0.83216     0.02389 0.5377
#Residuals   34    8.9652   0.26368               0.97611       
#Total       35    9.1846                         1.00000 


##Permanova on hellinger-transformed data 
(Permanova_OTU_wheat_2<-adonis(formula=PT_OTU_wheat_justOTU_hel ~ Inoculation, data=PT_OTU_wheat, permutations=99999, method="bray"))

#            Df  SumsOfSqs  MeanSqs   F.Model    R2     Pr(>F)
#Inoculation  1    0.1260   0.12601   1.0155     0.029   0.3943
#Residuals   34    4.2186    0.12408             0.971       
#Total       35    4.3446                       1.000 
  
##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated wheat fields


# Beta-diversity | permanova | relative abundance -------------------------

PT_OTU_totransf<-read.delim("PT_AMF_clustered_OTU_fortrans_R.txt", row.names=1, header=T)

##Relative abundance
PT_OTU_rel<-relativeAbundance_otu(PT_OTU_totransf)

##Import this database in Excel to have a database with the OTUs and 
write.table(PT_OTU_rel, file="OTU_table_relabundance.txt", row.names=TRUE, col.names = TRUE, sep = "\t")

##In excel, I flipped the database to have the samples as rows and the OTU and I added the treatments

OTU_relabundance<-read.delim("OTU_table_relAbundance.txt", row.names=1, header=T)

##Permanova with interactions

PT_OTUrel_justotu<-OTU_relabundance[,1:408]
PT_OTUrel_justotu_hel<-decostand(PT_OTUrel_justotu, "hel")

(adonis(formula=PT_OTUrel_justotu_hel ~ Crop*Inoculation*Growth_stage, data=OTU_relabundance, permutations=999, method="bray"))

#                               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Crop                            2    9.5867  4.7933  41.703 0.41715  0.001 ***
#Inoculation                     1    0.0733  0.0733   0.637 0.00319  0.759    
#Growth_stage                    1    0.3354  0.3354   2.918 0.01459  0.010 ** 
#Crop:Inoculation                2    0.3622  0.1811   1.576 0.01576  0.085 .  
#Crop:Growth_stage               1    0.2335  0.2335   2.031 0.01016  0.052 .  
#Inoculation:Growth_stage        1    0.1256  0.1256   1.093 0.00546  0.361    
#Crop:Inoculation:Growth_stage   1    0.0811  0.0811   0.705 0.00353  0.687    
#Residuals                     106   12.1836  0.1149         0.53015           
#Total                         115   22.9813                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##So there is a difference between crop and growth stage, but no interactions. 

##However, the pcoa reveals that the impact of growth stage could also only be due to the fact that soybeans have only been sampled at early stage. So we will re-run a permanova with only the wheat and the corn samples

OTU_relabundance<-read.delim("OTU_table_relAbundance.txt", row.names=1, header=T)

##Permanova on only wheat and corn samples

OTU_relabundance_WC<-OTU_relabundance[!OTU_relabundance$Crop=="Soy",]

OTUrel_WC_justotu<-OTU_relabundance_WC[,1:408]
OTUrel_WC_justotu_hel<-decostand(OTUrel_WC_justotu, "hel")

(adonis(formula=OTUrel_WC_justotu_hel ~ Crop*Inoculation*Growth_stage, data=OTU_relabundance_WC, permutations=999, method="bray"))

#                               Df SumsOfSqs MeanSqs F.Model   R2    Pr(>F)    
#Crop                           1    4.0009  4.0009  36.689 0.28643  0.001 ***
#Inoculation                    1    0.0655  0.0655   0.601 0.00469  0.820    
#Growth_stage                   1    0.3362  0.3362   3.083 0.02407  0.007 ** 
#Crop:Inoculation               1    0.1830  0.1830   1.678 0.01310  0.093 .  
#Crop:Growth_stage              1    0.2335  0.2335   2.141 0.01671  0.035 *  
#Inoculation:Growth_stage       1    0.1256  0.1256   1.152 0.00899  0.269    
#Crop:Inoculation:Growth_stage  1    0.0811  0.0811   0.743 0.00580  0.640    
#Residuals                     82    8.9421  0.1091         0.64019           
#Total                         89   13.9680                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




OTU_RA_corn<-OTU_relabundance[OTU_relabundance$Crop=="Corn",]
dim(OTU_RA_corn)
str(OTU_RA_corn$Crop)
str(OTU_RA_corn$OTU3_04026447528258e.32)

OTU_RA_soy<-OTU_relabundance[OTU_relabundance$Crop=="Soy",]
dim(OTU_RA_soy)
str(OTU_RA_soy$Crop)
str(OTU_RA_soy$OTU3_04026447528258e.32)

OTU_RA_wheat<-OTU_relabundance[OTU_relabundance$Crop=="Wheat",]
dim(OTU_RA_wheat)
str(OTU_RA_wheat$Crop)
str(OTU_RA_wheat$OTU3_04026447528258e.32)

###Corn

##Let's test the assumptions:
#1) Normal distribution: 
norm.her(OTU_RA_corn[,1:408], 0.01)
##Cannot do it cause some otu are absent in the corn fields

##Let's see if a log transformation can help
PT_OTU_corn_log<-decostand(PT_OTU_corn[,1:408], "log")
norm.her(PT_OTU_corn_log, 0.01)
##Cannot do it cause some otu are absent in the corn fields. Let's try to see the graph

#2) Linearity: yes

#3) Homogeneity of Variances: 
homeo.table(PT_OTU_corn[,1:408], PT_OTU_corn$Inoculation, 0.05)
##Doesn't work

##Let's do a permanova then!
##First let's devide the database
#WGRF_16S_otu_factors<-WGRF_otu_16S[,1:5]
OTU_RA_corn_justOTU<-OTU_RA_corn[,1:408]
OTU_RA_corn_justOTU_hel<-decostand(OTU_RA_corn_justOTU, "hel")

##Permanova on untransformed data 
(Permanova_OTU_rel_corn_1<-adonis(formula=OTU_RA_corn_justOTU ~ Inoculation, data=OTU_RA_corn, permutations=999, method="bray"))

#             Df SumsOfSqs  MeanSqs   F.Model       R2    Pr(>F)
#Inoculation  1    0.1112   0.11119    0.8472   0.01603   0.465
#Residuals   52    6.8246   0.13124             0.98397       
#Total       53    6.9358                       1.00000                       

##Permanova on hellinger-transformed data 
(Permanova_OTU_rel_corn_2<-adonis(formula=OTU_RA_corn_justOTU_hel ~ Inoculation, data=OTU_RA_corn, permutations=9999, method="bray"))

##Not really different, so let's keep the untransformed data

#             Df SumsOfSqs  MeanSqs   F.Model      R2   Pr(>F)
#Inoculation  1    0.1288   0.12879   1.2191  0.02291   0.2495
#Residuals   52    5.4936   0.10565           0.97709       
#Total       53    5.6224                     1.00000  

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated corn fields

###Soy

##Let's test the assumptions:
#1) Normal distribution: Doesn't work

#2) Linearity: yes

#3) Homogeneity of Variances: 
##Doesn't work

##Let's do a permanova then!
##First let's devide the database
OTU_RA_soy_justOTU<-OTU_RA_soy[,1:408]
OTU_RA_soy_justOTU_hel<-decostand(OTU_RA_soy_justOTU, "hel")

##Permanova on untransformed data 
(Permanova_OTU_rel_soy_1<-adonis(formula=OTU_RA_soy_justOTU ~ Inoculation, data=OTU_RA_soy, permutations=99999, method="bray"))

#             Df SumsOfSqs   MeanSqs  F.Model      R2   Pr(>F)
#Inoculation  1    0.2980   0.29799   1.4349  0.05642   0.1658
#Residuals   24    4.9841   0.20767           0.94358       
#Total       25    5.2821                     1.00000 


##Permanova on hellinger-transformed data 
(Permanova_OTUP_rel_soy_2<-adonis(formula=OTU_RA_soy_justOTU_hel ~ Inoculation, data=OTU_RA_soy, permutations=99999, method="bray"))

#             Df  SumsOfSqs   MeanSqs   F.Model       R2    Pr(>F)
#Inoculation  1      0.1860   0.18605    1.3775   0.05428   0.1524
#Residuals   24      3.2415   0.13506              0.94572       
#Total       25     3.4275                        1.00000                       

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated soy fields


###Wheat

##Let's test the assumptions:
#1) Normal distribution: Doesn't work

#2) Linearity: yes

#3) Homogeneity of Variances: 
##Doesn't work

##Let's do a permanova then!
##First let's devide the database
OTU_RA_wheat_justOTU<-OTU_RA_wheat[,1:408]
OTU_RA_wheat_justOTU_hel<-decostand(OTU_RA_wheat_justOTU, "hel")

##Permanova on untransformed data 
(Permanova_OTU_rel_wheat_1<-adonis(formula=OTU_RA_wheat_justOTU ~ Inoculation, data=OTU_RA_wheat, permutations=99999, method="bray"))

#             Df  SumsOfSqs   MeanSqs   F.Model      R2   Pr(>F)
#Inoculation  1      0.1708   0.17080   1.0224  0.02919   0.3568
#Residuals   34      5.6803   0.16707           0.97081       
#Total       35      5.8511                   1.00000


##Permanova on hellinger-transformed data 
(Permanova_OTU_rel_wheat_2<-adonis(formula=OTU_RA_wheat_justOTU_hel ~ Inoculation, data=OTU_RA_wheat, permutations=99999, method="bray"))

#             Df SumsOfSqs  MeanSqs   F.Model    R2   Pr(>F)
#Inoculation  1    0.1260   0.12601   1.0155  0.029   0.3903
#Residuals   34    4.2186   0.12408            0.971       
#Total       35    4.3446                   1.000 

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated wheat fields


# Beta-diversity | ANOVA table which OTU are different early-late ---------
OTU_relabundance<-read.delim("OTU_table_relAbundance.txt", row.names=1, header=T)
PT_OTUrel_corn<-OTU_relabundance[OTU_relabundance$Crop=="Corn",]


source("anova.1way.table.R")
anova_corn_table<-anova.1way.table(PT_OTUrel_corn[,1:408], PT_OTUrel_corn[,411],0.1)
View(anova_corn_table$resultat)

anova_corn_table_all<-anova_corn_table$resultat
class(anova_corn_table_all)
anova_corn_table_all<-as.data.frame(anova_corn_table_all)
(anova_corn_diff<-anova_corn_table_all[anova_corn_table_all$significance=="***TRUE***",])
dim(anova_corn_diff)

write.table(anova_corn_diff, file="anova_corn_diff", row.names=TRUE, col.names = TRUE, sep = "\t")

##So these 34 OTUs are different in early and late corn:
#Fvalue               pvalue significance
#OTU3_48355458501253e.31 2.82950191809067   0.0985436282206871   ***TRUE***
#OTU2_18035888669434e.32 11.4346209912191  0.00137533442639937   ***TRUE***
#OTU3_69348918298545e.32 15.4167922032178 0.000254941834878473   ***TRUE***
#OTU2_37057280116204e.31  3.6498762186276   0.0615946741111836   ***TRUE***
#OTU7_40221647409179e.32 12.7448770383698 0.000778463164500005   ***TRUE***
#OTU8_9421130128138e.32  4.36293208583371    0.041641471895664   ***TRUE***
#OTU1_41112008488097e.32 3.29100465909932   0.0754317029836086   ***TRUE***
#OTU2_63044671981184e.32 3.23225742515342   0.0780079133309209   ***TRUE***
#OTU2_57476572471429e.32  6.1988556192459   0.0160185821601681   ***TRUE***
#OTU4_10313003274393e.32 11.7382853423323  0.00120372642570359   ***TRUE***
#OTU5_79684382547281e.32  2.8920732757427   0.0949903008905494   ***TRUE***
#OTU1_56076894714513e.32 6.46141248353767   0.0140432325174662   ***TRUE***
#OTU5_04137900712937e.32 16.3367919423595 0.000175774757594038   ***TRUE***
#OTU9_86677944498954e.32 2.88095451966022    0.095611086555226   ***TRUE***
#OTU3_84871044258901e.31 7.75729205160152  0.00744604179886037   ***TRUE***
#OTU9_67096257934678e.32 15.3986933169079 0.000256829326454157   ***TRUE***
#OTU8_27852383179735e.32 6.53531929295868   0.0135351737504309   ***TRUE***
#OTU4_2824197729902e.32  6.13405047422809   0.0165503262713679   ***TRUE***
#OTU3_02855277831146e.32 5.30995297553775     0.02523036560269   ***TRUE***
#OTU9_21383855904061e.31 17.5835882850632 0.000107203039659978   ***TRUE***
#OTU2_82149776952224e.32  3.1182194609158   0.0832918763886221   ***TRUE***
#OTU2_07951248057061e.32  10.870091171814  0.00176608187396189   ***TRUE***
#OTU7_45568607458033e.32 8.57316195356144  0.00505218711185589   ***TRUE***
#OTU9_7162688448805e.32   7.1107867701681   0.0101879440422324   ***TRUE***
#OTU1_1045954652759e.32    5.424167999608   0.0237806401918363   ***TRUE***
#OTU5_90681954151361e.32 4.62047420641644   0.0362661569223906   ***TRUE***
#OTU5_31488713438041e.32 5.92711018401864    0.018378049714472   ***TRUE***
#OTU8_59407853036458e.32 5.81548920086583    0.019452250234918   ***TRUE***
#OTU8_07046339897152e.32  3.1812385441173   0.0803244475973694   ***TRUE***
#OTU2_89567166830228e.32 8.01551905826048  0.00657987442331571   ***TRUE***
#OTU2_29059676461656e.32 2.89147518870957   0.0950235783705599   ***TRUE***
#OTU4_90719687014263e.31 8.23233030870287  0.00593483988951569   ***TRUE***
#OTU1_15169030982265e.32  3.0852876809304   0.0848908748842753   ***TRUE***
#OTU3_19163376518745e.32  3.7391562894299    0.058602978439117   ***TRUE***


###Wheat: 
PT_OTUrel_wheat<-OTU_relabundance[OTU_relabundance$Crop=="Wheat",]

anova_wheat_table<-anova.1way.table(PT_OTUrel_wheat[,1:408], PT_OTUrel_wheat[,411],0.05)
View(anova_wheat_table$resultat)

anova_wheat_table_all<-anova_wheat_table$resultat
class(anova_wheat_table_all)
anova_wheat_table_all<-as.data.frame(anova_wheat_table_all)
(anova_wheat_diff<-anova_wheat_table_all[anova_wheat_table_all$significance=="***TRUE***",])
dim(anova_wheat_diff)

###So those 33 OTUs are different in early and late stage:
#Fvalue               pvalue significance
#OTU2_93133048647181e.32 7.77690621959544  0.00860526452923422   ***TRUE***
#OTU3_04026447528258e.32 8.60266184909651   0.0059708082233479   ***TRUE***
#OTU5_07156797541889e.32 7.10053601553698   0.0116995801461441   ***TRUE***
#OTU6_22033025644747e.32 4.86101478397915   0.0343362584246818   ***TRUE***
#OTU1_47033348157962e.32 8.87470954061443  0.00530482524456403   ***TRUE***
#OTU7_28828978033793e.32 6.06551831633206   0.0190062016722376   ***TRUE***
#OTU7_16498050959203e.32 5.21318046985293   0.0287881966300598   ***TRUE***
#OTU7_93234946063446e.32 6.44834750017504   0.0158477420576521   ***TRUE***
#OTU2_37057280116204e.31 4.55563523222325   0.0401040836248967   ***TRUE***
#OTU7_14192426647688e.32 4.45477948075214   0.0422367180177737   ***TRUE***
#OTU7_55726595050142e.32 12.1194508686706  0.00139054628901975   ***TRUE***
#OTU7_11676234027564e.32 5.67265231163383   0.0229713657339749   ***TRUE***
#OTU3_85420992683549e.32  5.7154596470618   0.0224985195944962   ***TRUE***
#OTU3_13149445883173e.32 7.35730120692028   0.0104027650107204   ***TRUE***
#OTU1_23951953357309e.32 4.39244260241733   0.0436171406370754   ***TRUE***
#OTU9_51677721475722e.32 6.23930880800867   0.0174950543668689   ***TRUE***
#OTU1_57912308107679e.32 4.76982829815811   0.0359569018394461   ***TRUE***
#OTU7_65631960428976e.32 7.15824055876209   0.0113936419218883   ***TRUE***
#OTU1_56076894714513e.32  4.1384714526767    0.049778979406192   ***TRUE***
#OTU4_63275406933015e.32 4.80697179214565   0.0352868449790304   ***TRUE***
#OTU3_94197835481386e.32 4.67410017170101   0.0377490733031509   ***TRUE***
#OTU3_9780919216923e.32   4.3067980583783    0.045595485329464   ***TRUE***
#OTU9_95475015671996e.32 5.46645151586859   0.0254057963805616   ***TRUE***
#OTU6_80633105527058e.31 5.20811209907112   0.0288607154878485   ***TRUE***
#OTU1_12332810494806e.32 14.6160547798383 0.000536009697362636   ***TRUE***
#OTU3_48214014301998e.32 11.0979433553337   0.0020916438202747   ***TRUE***
#OTU5_11820026326735e.32 12.2968316070789   0.0012968520172952   ***TRUE***
#OTU3_15650557947173e.32 7.18726994873571   0.0112430077906031   ***TRUE***
#OTU7_96215417158983e.31 4.79402658005132   0.0355188045303784   ***TRUE***
#OTU5_9647444296936e.32  7.36672443133888   0.0103582229028835   ***TRUE***
#OTU5_86636051000625e.32 10.9575977391072   0.0022143078194478   ***TRUE***
#OTU8_20890754066958e.32 7.57277970370861  0.00943389289153928   ***TRUE***
#OTU6_56566645999028e.32 4.75944955219754    0.036146617543641   ***TRUE***

write.table(anova_wheat_diff, file="anova_wheat_diff", row.names=TRUE, col.names = TRUE, sep = "\t")


# Beta-diversity | permanova | vtx ----------------------------------------

PT_vtx<-read.delim("PT_AMF_VTX_MA_RA_R.txt", row.names=1, header=T, dec=",")
str(PT_vtx)

vtx_corn<-PT_vtx[PT_vtx$Crop=="Corn",]
dim(vtx_corn)
str(vtx_corn$Crop)

vtx_soy<-PT_vtx[PT_vtx$Crop=="Soy",]
dim(vtx_soy)
str(vtx_soy$Crop)

vtx_wheat<-PT_vtx[PT_vtx$Crop=="Wheat",]
dim(vtx_wheat)
str(vtx_wheat$Crop)

##Let's do a permanova on corn!!
##First let's devide the database
vtx_corn_justOTU<-vtx_corn[,1:46]
vtx_corn_justOTU_hel<-decostand(vtx_corn_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=vtx_corn_justOTU ~ Inoculation, data=vtx_corn, permutations=9999, method="bray"))

###Not significantly different


##Permanova on hellinger-transformed data 
(adonis(formula=vtx_corn_justOTU_hel ~ Inoculation, data=vtx_corn, permutations=9999, method="bray"))

#not significantly different

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated corn fields

###Permanova - soy
##First let's devide the database
vtx_soy_justOTU<-vtx_soy[,1:46]
vtx_soy_justOTU_hel<-decostand(vtx_soy_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=vtx_soy_justOTU ~ Inoculation, data=vtx_soy, permutations=9999, method="bray"))

###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=vtx_soy_justOTU_hel ~ Inoculation, data=vtx_soy, permutations=9999, method="bray"))

#not significantly different

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated soy fields

###Permanova - wheat
##First let's devide the database
vtx_wheat_justOTU<-vtx_wheat[,1:46]
vtx_wheat_justOTU_hel<-decostand(vtx_wheat_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=vtx_wheat_justOTU ~ Inoculation, data=vtx_wheat, permutations=9999, method="bray"))
###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=vtx_wheat_justOTU_hel ~ Inoculation, data=vtx_wheat, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated wheat fields


# Beta-diversity | permanova | orders -------------------------------------

PT_orders<-read.delim("PT_AMF_orders_MA_RA_R.txt", row.names=1, header=T, dec=",")
str(PT_orders)

orders_corn<-PT_orders[PT_orders$Crop=="Corn",]
dim(orders_corn)
str(orders_corn$Crop)

orders_soy<-PT_orders[PT_orders$Crop=="Soy",]
dim(orders_soy)
str(orders_soy$Crop)

orders_wheat<-PT_orders[PT_orders$Crop=="Wheat",]
dim(orders_wheat)
str(orders_wheat$Crop)

##Let's do a permanova on corn!!
##First let's split the database
orders_corn_justOTU<-orders_corn[,1:3]
orders_corn_justOTU_hel<-decostand(order_corn_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=orders_corn_justOTU ~ Inoculation, data=orders_corn, permutations=9999, method="bray"))
###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=orders_corn_justOTU_hel ~ Inoculation, data=orders_corn, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the orders level is not different in inoculated or none-inoculated corn fields

###Permanova - soy
##First let's split the database
orders_soy_justOTU<-orders_soy[,1:3]
orders_soy_justOTU_hel<-decostand(orders_soy_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=orders_soy_justOTU ~ Inoculation, data=orders_soy, permutations=9999, method="bray"))
###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=orders_soy_justOTU_hel ~ Inoculation, data=orders_soy, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the OTU level is not different in inoculated or none-inoculated soy fields

###Permanova - wheat
##First let's devide the database
orders_wheat_justOTU<-orders_wheat[,1:3]
orders_wheat_justOTU_hel<-decostand(orders_wheat_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=orders_wheat_justOTU ~ Inoculation, data=orders_wheat, permutations=9999, method="bray"))
###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=orders_wheat_justOTU_hel ~ Inoculation, data=orders_wheat, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the orders level is not different in inoculated or none-inoculated wheat fields

# Beta-diversity | permanova | families -------------------------------------

PT_families<-read.delim("PT_AMF_families_MA_RA.txt", row.names=1, header=T, dec=",")
str(PT_families)

families_corn<-PT_families[PT_families$Crop=="Corn",]
dim(families_corn)
str(families_corn$Crop)

families_soy<-PT_families[PT_families$Crop=="Soy",]
dim(families_soy)
str(families_soy$Crop)

families_wheat<-PT_families[PT_families$Crop=="Wheat",]
dim(families_wheat)
str(families_wheat$Crop)

##Let's do a permanova on corn!!
##First let's split the database
families_corn_justOTU<-families_corn[,1:7]
families_corn_justOTU_hel<-decostand(families_corn_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=families_corn_justOTU ~ Inoculation, data=families_corn, permutations=9999, method="bray"))
###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=families_corn_justOTU_hel ~ Inoculation, data=families_corn, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the families level is not different in inoculated or none-inoculated corn fields

###Permanova - soy
##First let's split the database
families_soy_justOTU<-families_soy[,1:7]
families_soy_justOTU_hel<-decostand(families_soy_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=families_soy_justOTU ~ Inoculation, data=families_soy, permutations=9999, method="bray"))
###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=families_soy_justOTU_hel ~ Inoculation, data=families_soy, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the families level is not different in inoculated or none-inoculated soy fields

###Permanova - wheat
##First let's devide the database
families_wheat_justOTU<-families_wheat[,1:7]
families_wheat_justOTU_hel<-decostand(families_wheat_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=families_wheat_justOTU ~ Inoculation, data=families_wheat, permutations=99999, method="bray"))

#             Df SumsOfSqs  MeanSqs   F.Model      R2   Pr(>F)  
#Inoculation  1   0.09762 0.097624    2.671   0.07284   0.07775 .
#Residuals   34   1.24271 0.036550            0.92716          
#Total       35   1.34033                     1.00000          
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Permanova on hellinger-transformed data 
(adonis(formula=families_wheat_justOTU_hel ~ Inoculation, data=families_wheat, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the families level is not different in inoculated or none-inoculated wheat fields 

# Beta-diversity | permanova | genera -------------------------------------

PT_genera<-read.delim("PT_AMF_genera_MA_RA_R.txt", row.names=1, header=T, dec=",")
str(PT_genera)

genera_corn<-PT_genera[PT_genera$Crop=="Corn",]
dim(genera_corn)

genera_soy<-PT_genera[PT_genera$Crop=="Soy",]
dim(genera_soy)

genera_wheat<-PT_genera[PT_genera$Crop=="Wheat",]
dim(genera_wheat)

##Let's do a permanova on corn!!
##First let's split the database
genera_corn_justOTU<-genera_corn[,1:7]
genera_corn_justOTU_hel<-decostand(genera_corn_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=genera_corn_justOTU ~ Inoculation, data=genera_corn, permutations=9999, method="bray"))
###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=genera_corn_justOTU_hel ~ Inoculation, data=genera_corn, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the genera level is not different in inoculated or none-inoculated corn fields

###Permanova - soy
##First let's split the database
genera_soy_justOTU<-genera_soy[,1:7]
genera_soy_justOTU_hel<-decostand(genera_soy_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=genera_soy_justOTU ~ Inoculation, data=genera_soy, permutations=9999, method="bray"))
###Not significantly different

##Permanova on hellinger-transformed data 
(adonis(formula=genera_soy_justOTU_hel ~ Inoculation, data=genera_soy, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the families level is not different in inoculated or none-inoculated soy fields

###Permanova - wheat
##First let's devide the database
genera_wheat_justOTU<-genera_wheat[,1:7]
genera_wheat_justOTU_hel<-decostand(families_wheat_justOTU, "hel")

##Permanova on untransformed data 
(adonis(formula=genera_wheat_justOTU ~ Inoculation, data=genera_wheat, permutations=99999, method="bray"))

#Terms added sequentially (first to last)

#             Df SumsOfSqs  MeanSqs F.Model      R2  Pr(>F)  
#Inoculation  1   0.09762 0.097624  2.671   0.07284 0.07815 .
#Residuals   34   1.24271 0.036550           0.92716          
#Total       35   1.34033                   1.00000          
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Permanova on hellinger-transformed data 
(adonis(formula=genera_wheat_justOTU_hel ~ Inoculation, data=genera_wheat, permutations=9999, method="bray"))
#not significantly different

##Same answer with or without the transformation: The structure of the community at the genera level is not different in inoculated or none-inoculated wheat fields 


#ANOVAs - VTX00113 and VTX00114 - all treatments together -------------------------------------

##On vtx 00114
PT_vtx<-read.delim("PT_AMF_VTX_MA_RA_R.txt", row.names=1, header=T, dec=",")
str(PT_vtx)

##Vizualizing data
boxplot(PT_vtx$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=PT_vtx, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("abundance of Glomeraceae Glomus VTX 00114"), outline=F)

plot.design(PT_vtx$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=PT_vtx)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00114

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00114_inoculation<-lm(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=PT_vtx)
anova(anov_vtx00114_inoculation)
##Not significantly different at F(1,114)=2.0043 and p=0.1596

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation)
##Look rather good: Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation))
#W = 0.886, p-value = 6.035e-08
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=PT_vtx)
#Bartlett's K-squared = 1.7107, df = 1, p-value = 0.1909
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##Transforming data because assumption 1 and 2 are not respected
PT_vtx$sqrt_Glomeraceae_Glomus_sp._VTX00114<-sqrt(PT_vtx$Glomeraceae_Glomus_sp._VTX00114)
names(PT_vtx)

PT_vtx$log_Glomeraceae_Glomus_sp._VTX00114<-log(PT_vtx$Glomeraceae_Glomus_sp._VTX00114)
names(PT_vtx)

##Do a second anova
##replace -Inf by NA
PT_vtx$log_Glomeraceae_Glomus_sp._VTX00114[PT_vtx$log_Glomeraceae_Glomus_sp._VTX00114=="-Inf"]<-NA

anov_vtx00114_inoculation2<-lm(log_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=PT_vtx)
anova(anov_vtx00114_inoculation2)
##Not significantly different

par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation2)
##Not better

shapiro.test(resid(anov_vtx00114_inoculation2)) ##Still not normally distributed
bartlett.test(log_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=PT_vtx)  ##Equality of variances

##Doing ANova on sqrt transformed data
anov_vtx114_inoculation3<-lm(sqrt_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=PT_vtx)
anova(anov_vtx114_inoculation3)  
###Still not different

par(mfrow=c(2,2))
plot(anov_vtx114_inoculation3)

shapiro.test(resid(anov_vtx114_inoculation3))  ###The residus are bit normally distributed
bartlett.test(log_invsimpson ~ Inoculation, data=alphaD)  ##Equality of variances
##The residus are not normally distributed with both sqrt and log transformation, so let's do a non-parametric anova

anova.1way(PT_vtx$Glomeraceae_Glomus_sp._VTX00114 ~ PT_vtx$Inoculation, nperm=999)
##The abundance of vtx000114 is not different in inoculated or non-inoculated fields (F(1,114)=2.004; p=0.173)


##On vtx 00113
##Vizualizing data
boxplot(PT_vtx$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=PT_vtx, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("abundance of Glomeraceae Glomus VTX 00113"), outline=F)

plot.design(PT_vtx$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=PT_vtx)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00113_inoculation<-lm(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=PT_vtx)
anova(anov_vtx00113_inoculation)
##Not significantly different at F(1,114)=0.11531 and p=0.6963

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation))
#W = 0.94232, p-value = 8.204e-05
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=PT_vtx)
#Bartlett's K-squared = 0.013889, df = 1, p-value = 0.9062
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##Transforming data because assumption 1 and 2 are not respected
PT_vtx$sqrt_Glomeraceae_Glomus_sp._VTX00113<-sqrt(PT_vtx$Glomeraceae_Glomus_sp._VTX00113)
names(PT_vtx)

PT_vtx$log_Glomeraceae_Glomus_sp._VTX00113<-log(PT_vtx$Glomeraceae_Glomus_sp._VTX00113)
names(PT_vtx)

##Do a second anova

anov_vtx00113_inoculation2<-lm(log_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=PT_vtx)
anova(anov_vtx00113_inoculation2)
##Not significantly different

par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation2)
##Worst for normality

shapiro.test(resid(anov_vtx00113_inoculation2)) ##Still not normally distributed
bartlett.test(log_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=PT_vtx)  ##Equality of variances

##Doing ANova on sqrt transformed data
anov_vtx113_inoculation3<-lm(sqrt_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=PT_vtx)
anova(anov_vtx113_inoculation3)  
###Still not different

par(mfrow=c(2,2))
plot(anov_vtx113_inoculation3)

shapiro.test(resid(anov_vtx113_inoculation3))  ###The residus are not normally distributed
bartlett.test(sqrt_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=PT_vtx)  ##Equality of variances
##The residus are not normally distributed with both sqrt and log transformation, so let's do a non-parametric anova

anova.1way(PT_vtx$Glomeraceae_Glomus_sp._VTX00113 ~ PT_vtx$Inoculation, nperm=999)
##The abundance of vtx000113 is not different in inoculated or non-inoculated fields (F(1,114)=0.6963; p=0.686)


# ANOVAs - VTX00113 -split by crop and growing stage --------
##Let's redo-it by slitting the database per growing stage
##Corn
vtx_corn<-PT_vtx[PT_vtx$Crop=="Corn",]
vtx_corn_early<-vtx_corn[vtx_corn$Growing_stage=="Early",]
dim(vtx_corn_early)

vtx_corn_late<-vtx_corn[vtx_corn$Growing_stage=="Late",]
dim(vtx_corn_late)

##Soy
###Cannot do it, cause we only have one growing stage for soy

##Wheat:
vtx_wheat<-PT_vtx[PT_vtx$Crop=="Wheat",]
vtx_wheat_early<-vtx_wheat[vtx_wheat$Growing_stage=="Early",]
dim(vtx_wheat_early)

vtx_wheat_late<-vtx_wheat[vtx_wheat$Growing_stage=="Late",]
dim(vtx_wheat_late)

###VTX00113
###Corn-early
##On vtx 00113
##Vizualizing data
boxplot(vtx_corn_early$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_corn_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00113"))

plot.design(vtx_corn_early$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_corn_early)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00113_inoculation_CE<-lm(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_corn_early)
anova(anov_vtx00113_inoculation_CE)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_CE)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_CE))
#W = 0.94928, p-value = 0.2231
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_corn_early)
#Bartlett's K-squared = 0.26588, df = 1, p-value = 0.6061
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

#All the assumptions are respected, so the first anova can be interpreted: 
anova(anov_vtx00113_inoculation_CE)
#Response: Glomeraceae_Glomus_sp._VTX00113
#             Df  Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1    15.3   15.35  0.0203 0.8879
#Residuals   24 18157.2  756.55 

##The relative abundance of VTX00113 is not different in inoculated and non-inoculated in early-corn fields at F(1,24)=0.0203 and p=0.8879

###Corn-late
##On vtx 00113
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_corn_late$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_corn_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00113"))

plot.design(vtx_corn_late$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_corn_late)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00113_inoculation_CL<-lm(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_corn_late)
anova(anov_vtx00113_inoculation_CL)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_CL)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_CL))
#W = 0.96659, p-value = 0.4928
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_corn_late)
#Bartlett's K-squared = 0.54844, df = 1, p-value = 0.459
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

#All the assumptions are respected, so the first anova can be interpreted: 
anova(anov_vtx00113_inoculation_CL)
#Response: Glomeraceae_Glomus_sp._VTX00113
#             Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1    4.1    4.11   0.011 0.9172
#Residuals   26 9676.8  372.19

##The relative abundance of VTX00113 is not different in inoculated and non-inoculated in late-corn fields at F(1,26)=0.011 and p=0.9172

###Wheat-early
##On vtx 00113
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_wheat_early$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00113"))

plot.design(vtx_wheat_early$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_early)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00113_inoculation_WE<-lm(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_early)
anova(anov_vtx00113_inoculation_WE)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_WE)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_WE))
#W = 0.86214, p-value = 0.01064
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_early)
#Bartlett's K-squared = 2.4052, df = 1, p-value = 0.1209
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##The assumption of the normality of the residus is not respected, so we will try transformations: 

vtx_wheat_early$sqrt_Glomeraceae_Glomus_sp._VTX00113<-sqrt(vtx_wheat_early$Glomeraceae_Glomus_sp._VTX00113)
names(vtx_wheat_early)

vtx_wheat_early$log_Glomeraceae_Glomus_sp._VTX00113<-log(vtx_wheat_early$Glomeraceae_Glomus_sp._VTX00113)
names(vtx_wheat_early)

###Let's rerun anovas with log-transformed
anov_vtx00113_inoculation_WE2<-lm(log_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_early)
anova(anov_vtx00113_inoculation_WE2)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_WE2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_WE2))
#W = 0.86541, p-value = 0.01207
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_early)
#Bartlett's K-squared = 3.1029, df = 1, p-value = 0.07815
##the variance are equal

###Let's rerun anovas with log-transformed
anov_vtx00113_inoculation_WE3<-lm(sqrt_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_early)
anova(anov_vtx00113_inoculation_WE3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_WE3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_WE3))
#W = 0.95447, p-value = 0.469
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_early)
#Bartlett's K-squared = 2.2652, df = 1, p-value = 0.1323
##the variance are equal

#All the assumptions are respected with the sqrt transformation so the third anova can be interpreted: 
anova(anov_vtx00113_inoculation_WE3)
#Response: sqrt_Glomeraceae_Glomus_sp._VTX00113
#             Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1  0.503  0.5028  0.1546 0.6991
#Residuals   17 55.307  3.2534 

##The relative abundance of VTX00113 is not different in inoculated and non-inoculated in early-wheat fields at F(1,17)=0.1546 and p=0.6691

###Wheat-late
##On vtx 00113
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_wheat_late$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00113"))

plot.design(vtx_wheat_late$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_late)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00113_inoculation_WL<-lm(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_late)
anova(anov_vtx00113_inoculation_WL)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_WL)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_WL))
#W = 0.94299, p-value = 0.3554
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_wheat_late)
#Bartlett's K-squared = 0.27675, df = 1, p-value = 0.5988
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##The assumptions are all meet, so we can analyse the results of the first anova 
anova(anov_vtx00113_inoculation_WL)
#Response: Glomeraceae_Glomus_sp._VTX00113
#             Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1    1.9    1.92  0.0031 0.9567
#Residuals   15 9410.5  627.36 

##The relative abundance of VTX00113 is not different in inoculated and non-inoculated in late-wheat fields at F(1,15)=0.0031 and p=0.9567

###Soy
##On vtx 00113
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_soy$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_soy, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00113"))

plot.design(vtx_soy$Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_soy)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00113_inoculation_S<-lm(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_soy)
anova(anov_vtx00113_inoculation_S)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_S)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_S))
#W = 0.84695, p-value = 0.001242
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_soy)
#Bartlett's K-squared = 6.5137, df = 1, p-value = 0.0107
##the variances are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
vtx_soy$sqrt_Glomeraceae_Glomus_sp._VTX00113<-sqrt(vtx_soy$Glomeraceae_Glomus_sp._VTX00113)
names(vtx_soy)

vtx_soy$log_Glomeraceae_Glomus_sp._VTX00113<-log(vtx_soy$Glomeraceae_Glomus_sp._VTX00113)
names(vtx_soy)

###Let's rerun anovas with log-transformed
anov_vtx00113_inoculation_S2<-lm(log_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_soy)
anova(anov_vtx00113_inoculation_S2)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_S2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_S2))
#W = 0.90445, p-value = 0.01971
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_soy)
#Bartlett's K-squared = 0.1779, df = 1, p-value = 0.6732
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_vtx00113_inoculation_S3<-lm(sqrt_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_soy)
anova(anov_vtx00113_inoculation_S3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00113_inoculation_S3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00113_inoculation_S3))
#W = 0.94198, p-value = 0.1498
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Glomeraceae_Glomus_sp._VTX00113 ~ Inoculation, data=vtx_soy)
#Bartlett's K-squared = 2.3932, df = 1, p-value = 0.1219
##the variance are equal


##The assumptions are all meet, so we can analyse the results of the third anova 
anova(anov_vtx00113_inoculation_S3)
#Response: sqrt_Glomeraceae_Glomus_sp._VTX00113
#             Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1  1.538  1.5376  0.3764 0.5453
#Residuals   24 98.029  4.0845 

##The relative abundance of VTX00113 is not different in inoculated and non-inoculated in soy fields at F(1,24)=0.3764 and p=0.5453

# ANOVAs - VTX00114 - split by crop and growing stage --------
###VTX00114
###Corn-early
##On vtx 00114
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_corn_early$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_corn_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00114"))

plot.design(vtx_corn_early$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_corn_early)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00114_inoculation_CE<-lm(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_corn_early)
anova(anov_vtx00114_inoculation_CE)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_CE)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are  along the line: maybe the residus are normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_CE))
#W = 0.95228, p-value = 0.262
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_corn_early)
#Bartlett's K-squared = 0.11462, df = 1, p-value = 0.7349
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

#All the assumptions are respected, so the first anova can be interpreted: 
anova(anov_vtx00114_inoculation_CE)
#Response: Glomeraceae_Glomus_sp._VTX00114
#             Df   Sum Sq    Mean Sq F value Pr(>F)
#Inoculation  1 0.000003 0.00000288  0.0011 0.9733
#Residuals   24 0.060332 0.00251382      

##The relative abundance of VTX00113 is not different in inoculated and non-inoculated in early-corn fields at F(1,24)=0.0011 and p=0.9733

###Corn-late
##On vtx 00114
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_corn_late$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_corn_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00114"))

plot.design(vtx_corn_late$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_corn_late)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00114_inoculation_CL<-lm(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_corn_late)
anova(anov_vtx00114_inoculation_CL)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_CL)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are along the line: maybe the residus are normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_CL))
#W = 0.96091, p-value = 0.3664
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_corn_late)
#Bartlett's K-squared = 0.48656, df = 1, p-value = 0.4855
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

#All the assumptions are respected, so the first anova can be interpreted: 
anova(anov_vtx00114_inoculation_CL)
#Response: Glomeraceae_Glomus_sp._VTX00114
#             Df   Sum Sq    Mean Sq F value Pr(>F)
#Inoculation  1 0.000441 0.00044081  0.3396 0.5651
#Residuals   26 0.033748 0.00129799

##The relative abundance of VTX00113 is not different in inoculated and non-inoculated in late-corn fields at F(1,26)=0.3396 and p=0.5651

###Wheat-early
##On vtx 00114
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_wheat_early$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00114"))

plot.design(vtx_wheat_early$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_early)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00114

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00114_inoculation_WE<-lm(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_early)
anova(anov_vtx00114_inoculation_WE)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_WE)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_WE))
#W = 0.76374, p-value = 0.0003559
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_early)
#Bartlett's K-squared = 3.3524, df = 1, p-value = 0.06711
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

##The assumption of the normality of the residus is not respected, so we will try transformations: 

vtx_wheat_early$sqrt_Glomeraceae_Glomus_sp._VTX00114<-sqrt(vtx_wheat_early$Glomeraceae_Glomus_sp._VTX00114)
names(vtx_wheat_early)

vtx_wheat_early$log_Glomeraceae_Glomus_sp._VTX00114<-log(vtx_wheat_early$Glomeraceae_Glomus_sp._VTX00114)
names(vtx_wheat_early)

###Let's rerun anovas with log-transformed
##replace -Inf by NA
vtx_wheat_early$log_Glomeraceae_Glomus_sp._VTX00114[vtx_wheat_early$log_Glomeraceae_Glomus_sp._VTX00114=="-Inf"]<-NA

anov_vtx00114_inoculation_WE2<-lm(log_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_early)
anova(anov_vtx00114_inoculation_WE2)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_WE2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_WE2))
#W = 0.93484, p-value = 0.4342
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_early)
#Impossible to test
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_vtx00114_inoculation_WE3<-lm(sqrt_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_early)
anova(anov_vtx00114_inoculation_WE3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_WE3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_WE3))
#W = 0.88385, p-value = 0.02507
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_early)
#Bartlett's K-squared = 0.58441, df = 1, p-value = 0.4446
##the variance are equal

#All the assumptions are respected with the log transformation so the third anova can be interpreted: 
anova(anov_vtx00114_inoculation_WE2)
#Response: log_Glomeraceae_Glomus_sp._VTX00114
#             Df  Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1  0.0638 0.06382  0.0396 0.8462
#Residuals   10 16.0998 1.60998 

##The relative abundance of VTX00114 is not different in inoculated and non-inoculated in early-wheat fields at F(1,10)=0.0396 and p=0.8462

###Wheat-late
##On vtx 00114
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_wheat_late$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00114"))

plot.design(vtx_wheat_late$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_late)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00114

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00114_inoculation_WL<-lm(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_late)
anova(anov_vtx00114_inoculation_WL)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_WL)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_WL))
#W = 0.84489, p-value = 0.008971
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_late)
#Bartlett's K-squared = 1.9881, df = 1, p-value = 0.1585
##the variance are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

#The assumption of normality of the residus is not meet, so we will do the transformations
vtx_wheat_late$sqrt_Glomeraceae_Glomus_sp._VTX00114<-sqrt(vtx_wheat_late$Glomeraceae_Glomus_sp._VTX00114)
names(vtx_wheat_late)

vtx_wheat_late$log_Glomeraceae_Glomus_sp._VTX00114<-log(vtx_wheat_late$Glomeraceae_Glomus_sp._VTX00114)
names(vtx_wheat_late)

###Let's rerun anovas with log-transformed
##replace -Inf by NA
vtx_wheat_late$log_Glomeraceae_Glomus_sp._VTX00114[vtx_wheat_late$log_Glomeraceae_Glomus_sp._VTX00114=="-Inf"]<-NA

anov_vtx00114_inoculation_WL2<-lm(log_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_late)
anova(anov_vtx00114_inoculation_WL2)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_WL2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_WL2))
#W = 0.95499, p-value = 0.54
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_wheat_late)
#Bartlett's K-squared = 0.33742, df = 1, p-value = 0.5613
##the variance are equal

##The assumptions are all meet, so we can analyse the results of the second anova 
anova(anov_vtx00114_inoculation_WL2)
#Response: log_Glomeraceae_Glomus_sp._VTX00114
#             Df Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1  0.236 0.23597  0.2704 0.6107
#Residuals   15 13.090 0.87267

##The relative abundance of VTX00114 is not different in inoculated and non-inoculated in late-wheat fields at F(1,15)=0.2704 and p=0.6107

###Soy
##On vtx 00114
##Vizualizing data
par(mfrow=c(1,1))
boxplot(vtx_soy$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Relative abundance of Glomeraceae Glomus VTX 00114"))

plot.design(vtx_soy$Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_vtx00114_inoculation_S<-lm(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy)
anova(anov_vtx00114_inoculation_S)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_S)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_S))
#W = 0.77931, p-value = 7.99e-05
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy)
#Bartlett's K-squared = 10.978, df = 1, p-value = 0.0009221
##the variances are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
vtx_soy$sqrt_Glomeraceae_Glomus_sp._VTX00114<-sqrt(vtx_soy$Glomeraceae_Glomus_sp._VTX00114)
names(vtx_soy)

vtx_soy$log_Glomeraceae_Glomus_sp._VTX00114<-log(vtx_soy$Glomeraceae_Glomus_sp._VTX00114)
names(vtx_soy)

###Let's rerun anovas with log-transformed
vtx_soy$log_Glomeraceae_Glomus_sp._VTX00114[vtx_soy$log_Glomeraceae_Glomus_sp._VTX00114=="-Inf"]<-NA

anov_vtx00114_inoculation_S2<-lm(log_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy)
anova(anov_vtx00114_inoculation_S2)
#Significantly different at 0.0761 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_S2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_S2))
#W = 0.9861, p-value = 0.9895
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy)
#BBartlett's K-squared = 2.9327, df = 1, p-value = 0.0868
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_vtx00114_inoculation_S3<-lm(sqrt_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy)
anova(anov_vtx00114_inoculation_S3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_vtx00114_inoculation_S3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_vtx00114_inoculation_S3))
#W = 0.93484, p-value = 0.1011
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy)
#Bartlett's K-squared = 3.5695, df = 1, p-value = 0.05885
##the variance are equal


##The assumptions are all meet, so we can analyse the results of the second and thir anova 
anova(anov_vtx00114_inoculation_S2)
#Response: log_Glomeraceae_Glomus_sp._VTX00114
#             Df  Sum Sq Mean Sq F value  Pr(>F)  
#Inoculation  1  2.5583 2.55826  3.5514 0.07671 .
#Residuals   17 12.2460 0.72035                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##The relative abundance of VTX00114 (with log transformation) is different in inoculated and non-inoculated in soy fields at F(1,17)=3,5514 and p=0.07671. However by rejecting H0 with a p<0.1 I should have rejected H0 with the equality of variances so the variances were not equal and so the conditions of applications were not respected

anova(anov_vtx00114_inoculation_S3)
#Response: sqrt_Glomeraceae_Glomus_sp._VTX00114
#           Df   Sum Sq   Mean Sq F value Pr(>F)
#Inoculation  1 0.003553 0.0035529  0.5689  0.458
#Residuals   24 0.149881 0.0062450

####Let's try the non-parametric version: 
anova.1way(Glomeraceae_Glomus_sp._VTX00114 ~ Inoculation, data=vtx_soy, nperm=999)
#$anova.table
#             Df       Sum Sq      Mean Sq  F value Prob(param) Prob(perm)
#Inoculation  1 0.0005238764 0.0005238764 1.330525   0.2600703      0.264
#Residuals   24 0.0094496745 0.0003937364       NA          NA         NA


# ANOVAs - colonization data ----------------------------------------------
PT_colonization<-read.delim("PT_colonisation_R.txt", row.names=1, header=T, na.strings= "NA",  dec=",")
str(PT_colonization)

###Observed_colonization
##Vizualizing data
par(mfrow=c(1,1))
boxplot(PT_colonization$Observed_colonization ~ Inoculation, data=PT_colonization, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab="Level of observed colonization", outline=F)

plot.design(PT_colonization$Observed_colonization ~ Inoculation, data=PT_colonization)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_observed_col_inoculation<-lm(Observed_colonization ~ Inoculation, data=PT_colonization)
anova(anov_observed_col_inoculation)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation))
#W = 0.46063, p-value < 2.2e-16
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Observed_colonization ~ Inoculation, data=PT_colonization)
#Bartlett's K-squared = 14.371, df = 1, p-value = 0.0001501
##the variances are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
PT_colonization$sqrt_Observed_colonization<-sqrt(PT_colonization$Observed_colonization)
names(PT_colonization)

PT_colonization$log_Observed_colonization<-log(PT_colonization$Observed_colonization)
names(PT_colonization)

###Let's rerun anovas with log-transformed
##replace -Inf by NA
PT_colonization$log_Observed_colonization[PT_colonization$log_Observed_colonization=="-Inf"]<-NA

anov_observed_col_inoculation2<-lm(log_Observed_colonization ~ Inoculation, data=PT_colonization)
anova(anov_observed_col_inoculation2)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation2))
#W = 0.81792, p-value = 1.575e-09
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Observed_colonization ~ Inoculation, data=PT_colonization)
#Bartlett's K-squared = 2.525, df = 1, p-value = 0.1121
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_observed_col_inoculation3<-lm(sqrt_Observed_colonization ~ Inoculation, data=PT_colonization)
anova(anov_observed_col_inoculation3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation3))
#W = 0.8282, p-value = 2.657e-10
##Residus are still not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Observed_colonization ~ Inoculation, data=PT_colonization)
#Bartlett's K-squared = 2.3932, df = 1, p-value = 0.1219
##the variance are equal

###The assumptions are still not met (the normal distribution of the residus), so we are going to use a one-way anova

anova.1way(Observed_colonization ~ Inoculation, data=PT_colonization, nperm=999)
#$anova.table
#             Df     Sum Sq  Mean Sq   F value Prob(param) Prob(perm)
#Inoculation   1   3.449838 3.449838 0.6747437   0.4131172      0.359
#Residuals   114 582.860507 5.112811        NA          NA         NA

##So the observed colonization is not different whether the fields are inoculated or not-inoculated. 


###Spores
##Vizualizing data
par(mfrow=c(1,1))
boxplot(PT_colonization$spores ~ Inoculation, data=PT_colonization, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab="Spores (#spores/g)", outline=F)

plot.design(PT_colonization$spores ~ Inoculation, data=PT_colonization)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_spores_inoculation<-lm(spores ~ Inoculation, data=PT_colonization)
anova(anov_spores_inoculation)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation))
#W = 0.34597, p-value < 2.2e-16
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(spores ~ Inoculation, data=PT_colonization)
#Bartlett's K-squared = 2.5386, df = 1, p-value = 0.1111
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
PT_colonization$sqrt_spores<-sqrt(PT_colonization$spores)
names(PT_colonization)

PT_colonization$log_spores<-log(PT_colonization$spores)
names(PT_colonization)

###Let's rerun anovas with log-transformed
##replace -Inf by NA
PT_colonization$log_spores[PT_colonization$log_spores=="-Inf"]<-NA

anov_spores_inoculation2<-lm(log_spores ~ Inoculation, data=PT_colonization)
anova(anov_spores_inoculation2)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation2))
#W = 0.96601, p-value = 0.01678
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_spores ~ Inoculation, data=PT_colonization)
#Bartlett's K-squared = 0.83526, df = 1, p-value = 0.3608
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_spores_inoculation3<-lm(sqrt_spores ~ Inoculation, data=PT_colonization)
anova(anov_spores_inoculation3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation3))
#W = 0.71231, p-value = 1.409e-13
##Residus are still not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_spores ~ Inoculation, data=PT_colonization)
#Bartlett's K-squared = 0.3741, df = 1, p-value = 0.5408
##the variance are equal

###The assumptions are still not met (the normal distribution of the residus), so we are going to use a one-way anova

anova.1way(spores ~ Inoculation, data=PT_colonization, nperm=999)
#$anova.table
#             Df      Sum Sq   Mean Sq   F value Prob(param) Prob(perm)
#Inoculation   1    72176.45  72176.45 0.1886569     0.66488      0.635
#Residuals   111 42466434.73 382580.49        NA          NA         NA

##So the number of spores is not different whether the fields are inoculated or not-inoculated. 


# ANOVAs - colonization data - per crop and growing stage -----------------
##Let's redo-it by slitting the database per growing stage
##Corn
colonization_corn<-PT_colonization[PT_colonization$Crop=="Corn",]
dim(colonization_corn)
str(colonization_corn)

colonization_corn_early<-colonization_corn[colonization_corn$Growing_stage=="Early",]
dim(colonization_corn_early)
str(colonization_corn_early)

colonization_corn_late<-colonization_corn[colonization_corn$Growing_stage=="Late",]
dim(colonization_corn_late)
str(colonization_corn_late)

##Wheat:
colonization_wheat<-PT_colonization[PT_colonization$Crop=="Wheat",]
dim(colonization_wheat)
str(colonization_wheat)

colonization_wheat_early<-colonization_wheat[colonization_wheat$Growing_stage=="Early",]
dim(colonization_wheat_early)
str(colonization_wheat_early)

colonization_wheat_late<-colonization_wheat[colonization_wheat$Growing_stage=="Late",]
dim(colonization_wheat_late)
str(colonization_wheat_late)

colonization_soy<-PT_colonization[PT_colonization$Crop=="Soy",]
dim(colonization_soy)
str(colonization_soy)



###Observed_colonization - corn early
##Vizualizing data
par(mfrow=c(1,1))
boxplot(colonization_corn_early$Observed_colonization ~ Inoculation, data=colonization_corn_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab="Level of observed colonization", outline=F)

plot.design(Observed_colonization ~ Inoculation, data=colonization_corn_early)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_observed_col_inoculation_CE<-lm(Observed_colonization ~ Inoculation, data=colonization_corn_early)
anova(anov_observed_col_inoculation_CE)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_CE)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_CE))
#W = 0.73753, p-value = 1.809e-05
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Observed_colonization ~ Inoculation, data=colonization_corn_early)
#Bartlett's K-squared = Inf, df = 1, p-value < 2.2e-16
##the variances are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_corn_early$sqrt_Observed_colonization<-sqrt(colonization_corn_early$Observed_colonization)
names(colonization_corn_early)

colonization_corn_early$log_Observed_colonization<-log(colonization_corn_early$Observed_colonization)
names(colonization_corn_early)

###Let's rerun anovas with log-transformed
##replace -Inf by NA
colonization_corn_early$log_Observed_colonization[colonization_corn_early$log_Observed_colonization=="-Inf"]<-NA

anov_observed_col_inoculation_CE2<-lm(log_Observed_colonization ~ Inoculation, data=colonization_corn_early)
anova(anov_observed_col_inoculation2)
##Cannot run it cause to many NA and 0 

##Let's verify the assumptions of the ANOVA
#par(mfrow=c(2,2))
#plot(anov_observed_col_inoculation2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
#shapiro.test(resid(anov_observed_col_inoculation2))
#W = 0.81792, p-value = 1.575e-09
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
#bartlett.test(log_Observed_colonization ~ Inoculation, data=PT_colonization)
#Bartlett's K-squared = 2.525, df = 1, p-value = 0.1121
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_observed_col_inoculation_CE3<-lm(sqrt_Observed_colonization ~ Inoculation, data=colonization_corn_early)
anova(anov_observed_col_inoculation_CE3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_CE3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_CE3))
#W = 0.72605, p-value = 1.231e-05
##Residus are still not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Observed_colonization ~ Inoculation, data=colonization_corn_early)
#Bartlett's K-squared = Inf, df = 1, p-value < 2.2e-16
##the variance are equal

###The assumptions are still not met (the normal distribution of the residus), so we are going to use a one-way anova

anova.1way(Observed_colonization ~ Inoculation, data=colonization_corn_early, nperm=999)
#$anova.table
#             Df   Sum Sq   Mean Sq  F value Prob(param) Prob(perm)
#Inoculation  1 4.038462 4.0384615 12.11538 0.001932323      0.003
#Residuals   24 8.000000 0.3333333       NA          NA         NA

##So the observed colonization is different in the inoculated and non-inoculated corn fields at the early stage of growing. Roots growing in soils that have been inoculated have a higher observed colonization level

###Spores
##Vizualizing data
par(mfrow=c(1,1))
boxplot(spores ~ Inoculation, data=colonization_corn_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Spores (#spores/g)"))

plot.design(spores ~ Inoculation, data=colonization_corn_early)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_spores_inoculation_CE<-lm(spores ~ Inoculation, data=colonization_corn_early)
anova(anov_spores_inoculation_CE)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_CE)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_CE))
#W = 0.45291, p-value = 8.714e-09
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(spores ~ Inoculation, data=colonization_corn_early)
#Bartlett's K-squared = 1.0451, df = 1, p-value = 0.3066
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_corn_early$sqrt_spores<-sqrt(colonization_corn_early$spores)
names(colonization_corn_early)

colonization_corn_early$log_spores<-log(colonization_corn_early$spores)
names(colonization_corn_early)

###Let's rerun anovas with log-transformed
##replace -Inf by NA
#PT_colonization$log_spores[PT_colonization$log_spores=="-Inf"]<-NA

anov_spores_inoculation_CE2<-lm(log_spores ~ Inoculation, data=colonization_corn_early)
anova(anov_spores_inoculation_CE2)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_CE2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_CE2))
#W = 0.97392, p-value = 0.7261
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_spores ~ Inoculation, data=colonization_corn_early)
#Bartlett's K-squared = 0.028598, df = 1, p-value = 0.8657
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_spores_inoculation_CE3<-lm(sqrt_spores ~ Inoculation, data=colonization_corn_early)
anova(anov_spores_inoculation_CE3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_CE3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_CE3))
#W = 0.74078, p-value = 2.021e-05
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_spores ~ Inoculation, data=colonization_corn_early)
#Bartlett's K-squared = 0.00099309, df = 1, p-value = 0.9749
##the variance are equal

###The assumptions are met with the log-transformation, so we can interprete the second anova:

anova(anov_spores_inoculation_CE2)
#Response: log_spores
#             Df  Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1   5.584  5.5837  1.2644 0.2719
#Residuals   24 105.989  4.4162 

##So the number of spores is not different whether the fields are inoculated or not-inoculated at F(1,24)=1.2633; p=0.2719


###Observed_colonization - corn late
##Vizualizing data
par(mfrow=c(1,1))
boxplot(colonization_corn_late$Observed_colonization ~ Inoculation, data=colonization_corn_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Level of observed colonization"))

plot.design(Observed_colonization ~ Inoculation, data=colonization_corn_late)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_observed_col_inoculation_CL<-lm(Observed_colonization ~ Inoculation, data=colonization_corn_late)
anova(anov_observed_col_inoculation_CL)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_CL)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_CL))
#W = 0.87865, p-value = 0.00375
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Observed_colonization ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 0.70597, df = 1, p-value = 0.4008
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_corn_late$sqrt_Observed_colonization<-sqrt(colonization_corn_late$Observed_colonization)
names(colonization_corn_late)

colonization_corn_late$log_Observed_colonization<-log(colonization_corn_late$Observed_colonization)
names(colonization_corn_early)

###Let's rerun anovas with log-transformed
##replace -Inf by NA
#colonization_corn_early$log_Observed_colonization[colonization_corn_early$log_Observed_colonization=="-Inf"]<-NA

anov_observed_col_inoculation_CL2<-lm(log_Observed_colonization ~ Inoculation, data=colonization_corn_late)
anova(anov_observed_col_inoculation_CL2)
##Not significantly diff

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_CL2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_CL2))
#W = 0.81822, p-value = 0.0002276
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Observed_colonization ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 1.161, df = 1, p-value = 0.2813
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_observed_col_inoculation_CL3<-lm(sqrt_Observed_colonization ~ Inoculation, data=colonization_corn_late)
anova(anov_observed_col_inoculation_CL3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_CL3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_CL3))
#W = 0.85935, p-value = 0.001452
##Residus are still not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Observed_colonization ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 0.84849, df = 1, p-value = 0.357
##the variance are equal

###The assumptions are still not met (the normal distribution of the residus), so we are going to use a one-way anova

anova.1way(Observed_colonization ~ Inoculation, data=colonization_corn_late, nperm=999)
#$anova.table
#             Df      Sum Sq    Mean Sq   F value Prob(param) Prob(perm)
#Inoculation  1  0.06086957 0.06086957 0.1059371   0.7474228      0.712
#Residuals   26 14.93913043 0.57458194        NA          NA         NA

##So the observed colonization is not different in the inoculated and non-inoculated corn fields at the late stage of growing. 

###Spores
##Vizualizing data
par(mfrow=c(1,1))
boxplot(spores ~ Inoculation, data=colonization_corn_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Spores (#spores/g)"))

plot.design(spores ~ Inoculation, data=colonization_corn_late)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_spores_inoculation_CL<-lm(spores ~ Inoculation, data=colonization_corn_late)
anova(anov_spores_inoculation_CL)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_CL)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_CL))
#W = 0.76302, p-value = 2.541e-05
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(spores ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 0.5945, df = 1, p-value = 0.4407
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_corn_late$sqrt_spores<-sqrt(colonization_corn_late$spores)
names(colonization_corn_late)

colonization_corn_late$log_spores<-log(colonization_corn_late$spores)
names(colonization_corn_early)

###Let's rerun anovas with log-transformed
##replace -Inf by NA
#PT_colonization$log_spores[PT_colonization$log_spores=="-Inf"]<-NA

anov_spores_inoculation_CL2<-lm(log_spores ~ Inoculation, data=colonization_corn_late)
anova(anov_spores_inoculation_CL2)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_CL2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_CL2))
#W = 0.95699, p-value = 0.2952
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_spores ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 0.0088997, df = 1, p-value = 0.9248
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_spores_inoculation_CL3<-lm(sqrt_spores ~ Inoculation, data=colonization_corn_late)
anova(anov_spores_inoculation_CL3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_CL3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_CL3))
#W = 0.87071, p-value = 0.002521
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_spores ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 0.53635, df = 1, p-value = 0.4639
##the variance are equal

###The assumptions are met with the log-transformation, so we can interprete the second anova:

anova(anov_spores_inoculation_CL2)
#Response: log_spores
#             Df  Sum Sq Mean Sq F value Pr(>F)
#Inoculation  1   1.665  1.6653  0.3645 0.5512
#Residuals   26 118.786  4.5687

##So the number of spores is not different whether the fields are inoculated or not-inoculated at F(1,26)=0.3645; p=0.5512



###Observed_colonization - Wheat early
##Vizualizing data
par(mfrow=c(1,1))
boxplot(colonization_wheat_early$Observed_colonization ~ Inoculation, data=colonization_wheat_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Level of observed colonization"))

plot.design(Observed_colonization ~ Inoculation, data=colonization_wheat_early)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_observed_col_inoculation_WE<-lm(Observed_colonization ~ Inoculation, data=colonization_wheat_early)
anova(anov_observed_col_inoculation_WE)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_WE)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_WE))
#W = 0.78643, p-value = 0.0007308
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Observed_colonization ~ Inoculation, data=colonization_wheat_early)
#Bartlett's K-squared = 0.10755, df = 1, p-value = 0.7429
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_wheat_early$sqrt_Observed_colonization<-sqrt(colonization_wheat_early$Observed_colonization)
colonization_wheat_early$sqrt_Observed_colonization

colonization_wheat_early$log_Observed_colonization<-log(colonization_wheat_early$Observed_colonization)
colonization_wheat_early$log_Observed_colonization

###Let's rerun anovas with log-transformed
##replace -Inf by NA
colonization_wheat_early$log_Observed_colonization[colonization_wheat_early$log_Observed_colonization=="-Inf"]<-NA
colonization_wheat_early$log_Observed_colonization

anov_observed_col_inoculation_WE2<-lm(log_Observed_colonization ~ Inoculation, data=colonization_wheat_early)
anova(anov_observed_col_inoculation_WE2)
##Not significantly diff

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_WE2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_WE2))
#W = 0.60935, p-value = 3.325e-05
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Observed_colonization ~ Inoculation, data=colonization_wheat_early)
#Bartlett's K-squared = Inf, df = 1, p-value < 2.2e-16
##the variance are not equal

###Let's rerun anovas with sqrt-transformed
anov_observed_col_inoculation_WE3<-lm(sqrt_Observed_colonization ~ Inoculation, data=colonization_wheat_early)
anova(anov_observed_col_inoculation_WE3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_WE3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_WE3))
#W = 0.76243, p-value = 0.0003417
##Residus are still not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Observed_colonization ~ Inoculation, data=colonization_wheat_early)
#Bartlett's K-squared = 4.256e-05, df = 1, p-value = 0.9948
##the variance are equal

###The assumptions are still not met (the normal distribution of the residus), so we are going to use a one-way anova

anova.1way(Observed_colonization ~ Inoculation, data=colonization_wheat_early, nperm=999)
#$anova.table
#             Df      Sum Sq     Mean Sq     F value Prob(param) Prob(perm)
#Inoculation  1 0.001096491 0.001096491 0.002300095   0.9623076      0.904
#Residuals   17 8.104166667 0.476715686          NA          NA         NA

##So the observed colonization is not different in the inoculated and non-inoculated corn fields at the late stage of growing. 

###Spores
##Vizualizing data
par(mfrow=c(1,1))
boxplot(spores ~ Inoculation, data=colonization_wheat_early, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Spores (#spores/g)"))

plot.design(spores ~ Inoculation, data=colonization_wheat_early)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_spores_inoculation_WE<-lm(spores ~ Inoculation, data=colonization_wheat_early)
anova(anov_spores_inoculation_WE)
##significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_WE)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_WE))
#W = 0.60715, p-value = 1.222e-05
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(spores ~ Inoculation, data=colonization_wheat_early)
#Bartlett's K-squared = 12.858, df = 1, p-value = 0.0003361
##the variances are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_wheat_early$sqrt_spores<-sqrt(colonization_wheat_early$spores)
colonization_wheat_early$sqrt_spores

colonization_wheat_early$log_spores<-log(colonization_wheat_early$spores)
colonization_wheat_early$log_spores

###Let's rerun anovas with log-transformed
##replace -Inf by NA
colonization_wheat_early$sqrt_spores[colonization_wheat_early$sqrt_spores=="Inf"]<-NA
colonization_wheat_early$log_spores[colonization_wheat_early$log_spores=="Inf"]<-NA

anov_spores_inoculation_WE2<-lm(log_spores ~ Inoculation, data=colonization_wheat_early)
#anova(anov_spores_inoculation_CL2)
##to many NA cannot run

##Let's verify the assumptions of the ANOVA
#par(mfrow=c(2,2))
#plot(anov_spores_inoculation_CL2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
#shapiro.test(resid(anov_spores_inoculation_CL2))
#W = 0.95699, p-value = 0.2952
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
#bartlett.test(log_spores ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 0.0088997, df = 1, p-value = 0.9248
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_spores_inoculation_WE3<-lm(sqrt_spores ~ Inoculation, data=colonization_wheat_early)
anova(anov_spores_inoculation_WE3)
##significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_WE3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_WE3))
#W = 0.59249, p-value = 8.814e-06
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_spores ~ Inoculation, data=colonization_wheat_early)
#Bartlett's K-squared = 6.5784, df = 1, p-value = 0.01032
##the variance are equal

###The assumptions are not met, so we will run a non-parametric anova

anova.1way(spores ~ Inoculation, data=colonization_wheat_early, nperm=9999)
#$anova.table
#             Df    Sum Sq    Mean Sq  F value Prob(param) Prob(perm)
#Inoculation  1 0.3953725 0.39537255 9.133349 0.008577006     0.1046
#Residuals   15 0.6493333 0.04328889       NA          NA         NA

##So the number of spores is not different whether the fields are inoculated or not-inoculated at F(1,26)=9.113; p=0.1046



###Observed_colonization - Wheat late
##Vizualizing data
par(mfrow=c(1,1))
boxplot(Observed_colonization ~ Inoculation, data=colonization_wheat_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Level of observed colonization"))

plot.design(Observed_colonization ~ Inoculation, data=colonization_wheat_late)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_observed_col_inoculation_WL<-lm(Observed_colonization ~ Inoculation, data=colonization_wheat_late)
anova(anov_observed_col_inoculation_WL)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_WL)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_WL))
#W = 0.86181, p-value = 0.01632
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Observed_colonization ~ Inoculation, data=colonization_wheat_late)
#Bartlett's K-squared = 0.00023262, df = 1, p-value = 0.9878
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_wheat_late$sqrt_Observed_colonization<-sqrt(colonization_wheat_late$Observed_colonization)
colonization_wheat_late$sqrt_Observed_colonization

colonization_wheat_late$log_Observed_colonization<-log(colonization_wheat_late$Observed_colonization)
colonization_wheat_late$log_Observed_colonization

###Let's rerun anovas with log-transformed
##replace -Inf by NA
colonization_wheat_late$log_Observed_colonization[colonization_wheat_late$log_Observed_colonization=="-Inf"]<-NA
colonization_wheat_late$log_Observed_colonization

anov_observed_col_inoculation_WL2<-lm(log_Observed_colonization ~ Inoculation, data=colonization_wheat_late)
anova(anov_observed_col_inoculation_WL2)
##Not significantly diff

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_WL2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_WL2))
#W = 0.79703, p-value = 0.004587
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Observed_colonization ~ Inoculation, data=colonization_wheat_late)
#Bartlett's K-squared = 0.11234, df = 1, p-value = 0.7375
##the variance are not equal

###Let's rerun anovas with sqrt-transformed
anov_observed_col_inoculation_WL3<-lm(sqrt_Observed_colonization ~ Inoculation, data=colonization_wheat_late)
anova(anov_observed_col_inoculation_WL3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_WL3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_WL3))
#W = 0.76723, p-value = 0.0007461
##Residus are still not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Observed_colonization ~ Inoculation, data=colonization_wheat_late)
#Bartlett's K-squared = 0.062973, df = 1, p-value = 0.8019
##the variance are equal

###The assumptions are still not met (the normal distribution of the residus), so we are going to use a one-way anova

anova.1way(Observed_colonization ~ Inoculation, data=colonization_wheat_late, nperm=999)
#$anova.table
#             Df       Sum Sq      Mean Sq    F value Prob(param) Prob(perm)
#Inoculation  1 9.803922e-04 0.0009803922 0.00121369   0.9726683      0.975
#Residuals   15 1.211667e+01 0.8077777778         NA          NA         NA

##So the observed colonization is not different in the inoculated and non-inoculated wheat fields at the late stage of growing. 

###Spores
##Vizualizing data
par(mfrow=c(1,1))
boxplot(spores ~ Inoculation, data=colonization_wheat_late, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Spores (#spores/g)"))

plot.design(spores ~ Inoculation, data=colonization_wheat_late)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_spores_inoculation_WL<-lm(spores ~ Inoculation, data=colonization_wheat_late)
anova(anov_spores_inoculation_WL)
##significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_WL)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_WL))
#W = 0.48204, p-value = 1.511e-06
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(spores ~ Inoculation, data=colonization_wheat_late)
#Bartlett's K-squared = 11.619, df = 1, p-value = 0.0006528
##the variances are not equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_wheat_late$sqrt_spores<-sqrt(colonization_wheat_late$spores)
colonization_wheat_late$sqrt_spores

colonization_wheat_late$log_spores<-log(colonization_wheat_late$spores)
colonization_wheat_late$log_spores

###Let's rerun anovas with log-transformed
##replace -Inf by NA
colonization_wheat_early$sqrt_spores[colonization_wheat_early$sqrt_spores=="Inf"]<-NA
colonization_wheat_early$log_spores[colonization_wheat_early$log_spores=="Inf"]<-NA

anov_spores_inoculation_WL2<-lm(log_spores ~ Inoculation, data=colonization_wheat_late)
#anova(anov_spores_inoculation_CL2)
##to many NA cannot run

##Let's verify the assumptions of the ANOVA
#par(mfrow=c(2,2))
#plot(anov_spores_inoculation_CL2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
#shapiro.test(resid(anov_spores_inoculation_CL2))
#W = 0.95699, p-value = 0.2952
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
#bartlett.test(log_spores ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 0.0088997, df = 1, p-value = 0.9248
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_spores_inoculation_WL3<-lm(sqrt_spores ~ Inoculation, data=colonization_wheat_late)
anova(anov_spores_inoculation_WL3)
##significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_WL3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_WL3))
#W = 0.74269, p-value = 0.0005195
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_spores ~ Inoculation, data=colonization_wheat_late)
#Bartlett's K-squared = 3.1589, df = 1, p-value = 0.07551
##the variance are equal

###The assumptions are not met, so we will run a non-parametric anova

anova.1way(spores ~ Inoculation, data=colonization_wheat_late, nperm=9999)
#$anova.table
#             Df    Sum Sq  Mean Sq   F value Prob(param) Prob(perm)
#Inoculation  1  13.50114 13.50114 0.3523548   0.5622517     0.8252
#Residuals   14 536.43636 38.31688        NA          NA         NA

##So the number of spores is not different whether the fields are inoculated or not-inoculated at F(1,14)=0.3532; p=0.8252



###Observed_colonization - Soy
##Vizualizing data
par(mfrow=c(1,1))
boxplot(Observed_colonization ~ Inoculation, data=colonization_soy, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Level of observed colonization"))

plot.design(Observed_colonization ~ Inoculation, data=colonization_soy)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_observed_col_inoculation_S<-lm(Observed_colonization ~ Inoculation, data=colonization_soy)
anova(anov_observed_col_inoculation_S)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_S)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_S))
#W = 0.32043, p-value = 6.004e-10
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(Observed_colonization ~ Inoculation, data=colonization_soy)
#Bartlett's K-squared = 15.201, df = 1, p-value = 9.665e-05
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_soy$sqrt_Observed_colonization<-sqrt(colonization_soy$Observed_colonization)
colonization_soy$sqrt_Observed_colonization

colonization_soy$log_Observed_colonization<-log(colonization_soy$Observed_colonization)
colonization_soy$log_Observed_colonization

###Let's rerun anovas with log-transformed
##replace -Inf by NA
#colonization_wheat_late$log_Observed_colonization[colonization_wheat_late$log_Observed_colonization=="-Inf"]<-NA

anov_observed_col_inoculation_S2<-lm(log_Observed_colonization ~ Inoculation, data=colonization_soy)
anova(anov_observed_col_inoculation_S2)
##Not significantly diff

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_S2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_S2))
#W = 0.54087, p-value = 6.582e-08
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(log_Observed_colonization ~ Inoculation, data=colonization_soy)
#Bartlett's K-squared = 4.1901, df = 1, p-value = 0.04066
##the variance are not equal

###Let's rerun anovas with sqrt-transformed
anov_observed_col_inoculation_S3<-lm(sqrt_Observed_colonization ~ Inoculation, data=colonization_soy)
anova(anov_observed_col_inoculation_S3)
##Not significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_observed_col_inoculation_S3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_observed_col_inoculation_S3))
#W = 0.40802, p-value = 3.368e-09
##Residus are still not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_Observed_colonization ~ Inoculation, data=colonization_soy)
#Bartlett's K-squared = 8.9144, df = 1, p-value = 0.002829
##the variance are equal

###The assumptions are still not met (the normal distribution of the residus), so we are going to use a one-way anova

anova.1way(Observed_colonization ~ Inoculation, data=colonization_soy, nperm=999)
#$anova.table
#             Df     Sum Sq   Mean Sq   F value Prob(param) Prob(perm)
#Inoculation  1   4.928205  4.928205 0.2975271    0.590471      0.627
#Residuals   24 397.533333 16.563889        NA          NA         NA

##So the observed colonization is not different in the inoculated and non-inoculated soy fields. 

###Spores
##Vizualizing data
par(mfrow=c(1,1))
boxplot(spores ~ Inoculation, data=colonization_soy, col=c("skyblue1","skyblue4"), xlab="Inoculation", ylab=("Spores (#spores/g)"))

plot.design(spores ~ Inoculation, data=colonization_soy)
##Probably not different, but the plot.design graph showed that plots that were  inoculated seems to have a higher abundance of the vtx00113

##Trying anova
##H0: u1=u2=u3=...=un
##H1: at least one u is different

anov_spores_inoculation_S<-lm(spores ~ Inoculation, data=colonization_soy)
anova(anov_spores_inoculation_WL)
##significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_S)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_S))
#W = 0.68068, p-value = 2.917e-06
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(spores ~ Inoculation, data=colonization_soy)
#Bartlett's K-squared = 0.17062, df = 1, p-value = 0.6796
##the variances are equal

##Assumption 3: Independance of the observations
##Yes, cause each factor is independant of the other and that one can stay constant while the other is varying. 

###The assumptions are not met so we will try transforming the data
colonization_soy$sqrt_spores<-sqrt(colonization_soy$spores)
colonization_soy$sqrt_spores

colonization_soy$log_spores<-log(colonization_soy$spores)
colonization_soy$log_spores

###Let's rerun anovas with log-transformed
##replace -Inf by NA
#colonization_wheat_early$sqrt_spores[colonization_wheat_early$sqrt_spores=="Inf"]<-NA
colonization_soy$log_spores[colonization_soy$log_spores=="Inf"]<-NA

anov_spores_inoculation_S2<-lm(log_spores ~ Inoculation, data=colonization_soy)
#anova(anov_spores_inoculation_CL2)
##to many NA cannot run

##Let's verify the assumptions of the ANOVA
#par(mfrow=c(2,2))
#plot(anov_spores_inoculation_CL2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
#shapiro.test(resid(anov_spores_inoculation_CL2))
#W = 0.95699, p-value = 0.2952
##Residus are distributed normally

##Assumption 2: Homoegeneity of the variances
#bartlett.test(log_spores ~ Inoculation, data=colonization_corn_late)
#Bartlett's K-squared = 0.0088997, df = 1, p-value = 0.9248
##the variance are equal

###Let's rerun anovas with sqrt-transformed
anov_spores_inoculation_S3<-lm(sqrt_spores ~ Inoculation, data=colonization_soy)
anova(anov_spores_inoculation_S3)
##significantly different 

##Let's verify the assumptions of the ANOVA
par(mfrow=c(2,2))
plot(anov_spores_inoculation_S3)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscédacity might be respected
##Graph Normal Q-Q: the points are not along the line: maybe the residus are not normally distributed

#Assupmtion 1: Normal distribution of the residues
shapiro.test(resid(anov_spores_inoculation_S3))
#W = 0.86536, p-value = 0.002867
##Residus are not distributed normally

##Assumption 2: Homoegeneity of the variances
bartlett.test(sqrt_spores ~ Inoculation, data=colonization_soy)
#Bartlett's K-squared = 0.24669, df = 1, p-value = 0.6194
##the variance are equal

###The assumptions are not met, so we will run a non-parametric anova

anova.1way(spores ~ Inoculation, data=colonization_soy, nperm=9999)
#$anova.table
#             Df     Sum Sq   Mean Sq  F value Prob(param) Prob(perm)
#Inoculation  1   2545.232  2545.232 0.195673   0.6621972     0.6425
#Residuals   24 312181.883 13007.578       NA          NA         NA

##So the number of spores is not different whether the fields are inoculated or not-inoculated at F(1,14)=0.3532; p=0.8252



# PCoA - OTUs - raw abundance-------------------------------------------------------------
##Importing database
PT_OTU<-read.delim("PT_AMF_clustered_OTU_R.txt", row.names=1, header=T)
str(PT_OTU)
dim(PT_OTU)

##PCoA on transformed data
PT_OTU_justotu<-PT_OTU[,1:408]
PT_OTU_justotu_hel<-decostand(PT_OTU_justotu, "hel")

##Calculating Bray-Curtis dissimilarity matrix on the hellinger transformed data
PT_OTU_hel_bray<-vegdist(PT_OTU_justotu_hel, method="bray")

##Calculating PCoA

PT_OTU_hel.pcoa<-pcoa(dist(PT_OTU_hel_bray))
PT_OTU_hel.pcoa

###How many axes represent more variability:
(Broken_stick<-PT_OTU_hel.pcoa$values$Broken_stick)
(meanBS<-mean(Broken_stick))
(largeBS<-Broken_stick[Broken_stick>meanBS])
(n=length(largeBS))
##Answer: we could plot up to 42 axis


##Ploting the PCoAs - with inoculation as empty circles
(Scores<-PT_OTU_hel.pcoa$vectors[,c(1,2)])

(Crop<-PT_OTU$Crop)
str(Crop)
crop_s<-as.character(Crop)
str(crop_s)

(Inoculation<-PT_OTU$Inoculation)
Inoculation_s<-as.character(Inoculation)
str(Inoculation_s)

(Scores_1<-cbind(Scores, crop_s, Inoculation_s))
(Scores.df<-as.data.frame(Scores_1))
str(Scores.df)
Scores.df$Axis.1<-as.numeric(Scores.df$Axis.1)
Scores.df$Axis.2<-as.numeric(Scores.df$Axis.2)
str(Scores.df)

##PCoA with Vegetation as colors
with(Scores.df, levels(Inoculation_s))
#colvec <- c("#7b3294", "#008837") #"#88000d"
colvec <- c("#e41a1c", "#377eb8", "#4daf4a")
pchvec<-c(1,16)
col.pcoa<-colvec[Scores.df[,3]]
pch.pcoa<-pchvec[Scores.df[,4]]
#colvec <- c("#d8b365", "#5ab4ac")

#colvec <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
#pchvec<-c(15,16,17, 18)

#### Code taken from https://sites.google.com/site/manabusakamoto/home/r-tutorials/r-tutorial-3

plot(PT_OTU_hel.pcoa$vectors[,1], PT_OTU_hel.pcoa$vectors[,2],  pch=pch.pcoa, cex=1.7, col=col.pcoa,  bty="n", xlab="Component 1", ylab="Component 2")
legend("top", legend = levels(Scores.df[,3]), pch=16, col=colvec, cex=1.2, bty="n" )
#text(PT_OTU_hel.pcoa$vectors[,1], PT_OTU_hel.pcoa$vectors[,2], row.names(PT_OTU_hel.pcoa$vectors), cex=0.6, pos=4, col="blue")

# PCoA - OTUs - relative abundance-------------------------------------------------------------
OTU_relabundance<-read.delim("OTU_table_relAbundance.txt", row.names=1, header=T)

##Transforming data:
PT_OTUrel_justotu<-OTU_relabundance[,1:408]
PT_OTUrel_justotu_hel<-decostand(PT_OTUrel_justotu, "hel")

##Calculating Bray-Curtis dissimilarity matrix on the hellinger transformed data
PT_OTUrel_hel_bray<-vegdist(PT_OTUrel_justotu_hel, method="bray")

##Calculating PCoA

PT_OTUrel_hel.pcoa<-pcoa(dist(PT_OTUrel_hel_bray))
PT_OTUrel_hel.pcoa

###How many axes represent more variability:
(Broken_stick<-PT_OTUrel_hel.pcoa$values$Broken_stick)
(meanBS<-mean(Broken_stick))
(largeBS<-Broken_stick[Broken_stick>meanBS])
(n=length(largeBS))
##Answer: we could plot up to 42 axis


##Ploting the PCoAs - with inoculation as empty circles
(Scores<-PT_OTUrel_hel.pcoa$vectors[,c(1,2)])

(Crop<-PT_OTU$Crop)
str(Crop)
crop_s<-as.character(Crop)
str(crop_s)

(Inoculation<-PT_OTU$Inoculation)
Inoculation_s<-as.character(Inoculation)
str(Inoculation_s)

(Growth_stage<-OTU_relabundance$Growth_stage)
Growth_stage_s<-as.character(Growth_stage)

(Scores_1<-cbind(Scores, crop_s, Inoculation_s, Growth_stage_s))
(Scores.df<-as.data.frame(Scores_1))
str(Scores.df)
Scores.df$Axis.1<-as.numeric(Scores.df$Axis.1)
Scores.df$Axis.2<-as.numeric(Scores.df$Axis.2)
str(Scores.df)

##PCoA with crop as colors and incoluation as fill or not
with(Scores.df, levels(Inoculation_s))
#colvec <- c("#7b3294", "#008837") #"#88000d"
colvec <- c("#e41a1c", "#377eb8", "#4daf4a")
pchvec<-c(1,16)
col.pcoa<-colvec[Scores.df[,3]]
pch.pcoa<-pchvec[Scores.df[,4]]
#colvec <- c("#d8b365", "#5ab4ac")

#colvec <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
#pchvec<-c(15,16,17, 18)

#### Code taken from https://sites.google.com/site/manabusakamoto/home/r-tutorials/r-tutorial-3

plot(PT_OTUrel_hel.pcoa$vectors[,1], PT_OTUrel_hel.pcoa$vectors[,2],  pch=pch.pcoa, cex=1.7, col=col.pcoa,  bty="n", xlab="Component 1", ylab="Component 2")
legend("top", legend = levels(Scores.df[,3]), pch=16, col=colvec, cex=1.2, bty="n" )
#text(PT_OTU_hel.pcoa$vectors[,1], PT_OTU_hel.pcoa$vectors[,2], row.names(PT_OTU_hel.pcoa$vectors), cex=0.6, pos=4, col="blue")

##PCoA with crop as colors and growth_stage as fill or not
with(Scores.df, levels(Inoculation_s))
#colvec <- c("#7b3294", "#008837") #"#88000d"
colvec <- c("#e41a1c", "#377eb8", "#4daf4a")
pchvec<-c(1,16)
col.pcoa<-colvec[Scores.df[,3]]
pch.pcoa<-pchvec[Scores.df[,5]]
#colvec <- c("#d8b365", "#5ab4ac")

#colvec <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
#pchvec<-c(15,16,17, 18)

#### Code taken from https://sites.google.com/site/manabusakamoto/home/r-tutorials/r-tutorial-3

plot(PT_OTUrel_hel.pcoa$vectors[,1], PT_OTUrel_hel.pcoa$vectors[,2],  pch=pch.pcoa, cex=1.7, col=col.pcoa,  bty="n", xlab="Component 1", ylab="Component 2")
legend("top", legend = levels(Scores.df[,3]), pch=16, col=colvec, cex=1.2, bty="n" )
#text(PT_OTU_hel.pcoa$vectors[,1], PT_OTU_hel.pcoa$vectors[,2], row.names(PT_OTU_hel.pcoa$vectors), cex=0.6, pos=4, col="blue")


###Refaire la PcoA - just wheat and corn

OTUrel_WC_justotu<-OTU_relabundance_WC[,1:408]
OTUrel_WC_justotu_hel<-decostand(OTUrel_WC_justotu, "hel")

##Calculating Bray-Curtis dissimilarity matrix on the hellinger transformed data
PT_OTUrel_WC_hel_bray<-vegdist(OTUrel_WC_justotu_hel, method="bray")

##Calculating PCoA

PT_OTUrel_WC_hel.pcoa<-pcoa(dist(PT_OTUrel_WC_hel_bray))
PT_OTUrel_WC_hel.pcoa

###How many axes represent more variability:
(Broken_stick<-PT_OTUrel_WC_hel.pcoa$values$Broken_stick)
(meanBS<-mean(Broken_stick))
(largeBS<-Broken_stick[Broken_stick>meanBS])
(n=length(largeBS))
##Answer: we could plot up to 33 axis


##Ploting the PCoAs - with inoculation as empty circles
(Scores<-PT_OTUrel_WC_hel.pcoa$vectors[,c(1,2)])

(Crop<-OTU_relabundance_WC$Crop)
str(Crop)
crop_s<-as.character(Crop)
str(crop_s)

(Inoculation<-OTU_relabundance_WC$Inoculation)
Inoculation_s<-as.character(Inoculation)
str(Inoculation_s)

(Growth_stage<-OTU_relabundance_WC$Growth_stage)
Growth_stage_s<-as.character(Growth_stage)

(Scores_1<-cbind(Scores, crop_s, Inoculation_s, Growth_stage_s))
(Scores.df<-as.data.frame(Scores_1))
str(Scores.df)
Scores.df$Axis.1<-as.numeric(Scores.df$Axis.1)
Scores.df$Axis.2<-as.numeric(Scores.df$Axis.2)
str(Scores.df)

##PCoA with crop as colors and incoluation as fill or not
with(Scores.df, levels(Inoculation_s))
#colvec <- c("#7b3294", "#008837") # "#88000d"
colvec <- c("#e41a1c", "#377eb8") #  "#4daf4a"
pchvec<-c(1,16)
col.pcoa<-colvec[Scores.df[,3]]
pch.pcoa<-pchvec[Scores.df[,5]]
#colvec <- c("#d8b365", "#5ab4ac")

#colvec <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
#pchvec<-c(15,16,17, 18)

#### Code taken from https://sites.google.com/site/manabusakamoto/home/r-tutorials/r-tutorial-3

plot(PT_OTUrel_WC_hel.pcoa$vectors[,1], PT_OTUrel_WC_hel.pcoa$vectors[,2],  pch=pch.pcoa, cex=1.7, col=col.pcoa,  bty="n", xlab="Component 1", ylab="Component 2")
legend("bottomright", legend = levels(Scores.df[,3]), pch=16, col=colvec, cex=1.2, bty="n" )
#text(PT_OTU_hel.pcoa$vectors[,1], PT_OTU_hel.pcoa$vectors[,2], row.names(PT_OTU_hel.pcoa$vectors), cex=0.6, pos=4, col="blue")

##PCoA with crop as colors and growth_stage as fill or not
with(Scores.df, levels(Inoculation_s))
#colvec <- c("#7b3294", "#008837") #"#88000d"
colvec <- c("#e41a1c", "#377eb8", "#4daf4a")
pchvec<-c(1,16)
col.pcoa<-colvec[Scores.df[,3]]
pch.pcoa<-pchvec[Scores.df[,5]]
#colvec <- c("#d8b365", "#5ab4ac")

#colvec <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
#pchvec<-c(15,16,17, 18)

#### Code taken from https://sites.google.com/site/manabusakamoto/home/r-tutorials/r-tutorial-3

plot(PT_OTUrel_hel.pcoa$vectors[,1], PT_OTUrel_hel.pcoa$vectors[,2],  pch=pch.pcoa, cex=1.7, col=col.pcoa,  bty="n", xlab="Component 1", ylab="Component 2")
legend("top", legend = levels(Scores.df[,3]), pch=16, col=colvec, cex=1.2, bty="n" )
#text(PT_OTU_hel.pcoa$vectors[,1], PT_OTU_hel.pcoa$vectors[,2], row.names(PT_OTU_hel.pcoa$vectors), cex=0.6, pos=4, col="blue")

