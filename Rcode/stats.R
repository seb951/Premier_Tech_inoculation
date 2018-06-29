####setting things up ----
setwd("/Users/jerry/Documents/CSBQ/hijri/Premier_Tech")

#required packages
library(vegan)
library(nlme)
library(lme4)
library(ape)
library(pals)
###stat analyses preparing OTU table and sampling design----
#load OTUs
OTU = read.table("results/OTU.table.VT.merged.ordered.top.tsv",row.names = 1, sep = "\t",header =T,stringsAsFactors = F)

#taxonomy in a separate vector
taxo = OTU[,ncol(OTU)]
OTU = data.frame(t(OTU[,-c(ncol(OTU))]))

#order based on rownames
OTU = OTU[order(rownames(OTU)),]

#normalize based on total number of OTU in sample (given that some seq libraries may have started with more sequences...)
OTU.norm = OTU/rowSums(OTU)

#load sampling design
design = read.table("results/PT_crop_design.txt", sep = "\t", stringsAsFactors = F, row.names = 1,header = T)

#order based on rownames
design = design[order(rownames(design)),]

#remove samples that were discarded in the design....
keep = rownames(design) %in% rownames(OTU.norm)
design.keep = design[keep,]

#####Is R. irregularis (VTX00113) more abundant in inoculated soils? ----

#VTX00113 is the most abundant VTX, and a Glomus spp..
#prepare a matrix with only VTX00113
OTU.norm.VTX00113 = cbind(OTU.norm[,1],design.keep)
colnames(OTU.norm.VTX00113)[1] = colnames(OTU.norm)[1]

#linear mixed effect model (block is random)
lmm1 <- lme(VTX00113~treatment+species+growing_stage,data = OTU.norm.VTX00113,random = ~1|bloc, method = "ML")
anova(lmm1)

#              numDF denDF  F-value p-value
#Intercept)       1   102 316.02848  <.0001
#treatment         1     8   0.25388  0.6279
#species           2     8  32.53176  0.0001 #only significant effect
#growing_stage     1   102   2.78624  0.0981
#

dev.new()
par(mfrow = c(2,2))
#could add a italized expression to the names
boxplot(OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Wheat",1]~OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Wheat",5],main = "VTX00113 - Wheat",ylab = "OTU relative abundance",names = expression(italic(control),italic(inoculated)))
boxplot(OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Corn",1]~OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Corn",5],main = "VTX00113 - Corn",ylab = "OTU relative abundance",names = expression(italic(control),italic(inoculated)))
boxplot(OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Soy",1]~OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Soy",5],main = "VTX00113 - Soy",ylab = "OTU relative abundance",names = expression(italic(control),italic(inoculated)))

dev.print(device=pdf, "figures/figure1_VTX00113.pdf", onefile=FALSE)
dev.off()


#####barplots of VTX ----
#first group things according to incol and species.
OTU.norm.barplot = cbind(colMeans(OTU.norm[design.keep[,1] == "Wheat" & design.keep[,4] == "inoculated", ]), 
      colMeans(OTU.norm[design.keep[,1] == "Wheat" & design.keep[,4] == "control", ]),
      colMeans(OTU.norm[design.keep[,1] == "Corn" & design.keep[,4] == "inoculated", ]), 
      colMeans(OTU.norm[design.keep[,1] == "Corn" & design.keep[,4] == "control", ]),
      colMeans(OTU.norm[design.keep[,1] == "Soy" & design.keep[,4] == "inoculated", ]), 
      colMeans(OTU.norm[design.keep[,1] == "Soy" & design.keep[,4] == "control", ]))
      
colnames(OTU.norm.barplot) = c("wheat.inoc", "wheat.ctl","corn.inoc","corn.ctl","soy.inoc","soy.ctl")      

dev.new()
par(mar=c(16,4,4,2))
barplot(OTU.norm.barplot,beside = F,font = 3, axisnames = T,ylab = "Relative OTU abundance",col = c(cols25(),rep("grey",60)), las = 3, xpd = T,names.arg =  rep(c("Inoculated","Control"),3))
text(y = rep(1.08,6), x = c(1.3,3.6,6.2), labels = c("Wheat","Corn","Soy"), font = 2,cex = 2,xpd = T)

legend(0.2,-0.4,fill = cols25(),legend = paste(rownames(OTU.norm.barplot)[1:10],taxo[1:10],sep = "_"),cex = 0.75,xpd =T)
dev.print(device=pdf, "figures/figure2_OTUabundance.pdf", onefile=FALSE)
dev.off()

##barplot by genera / order / family ----
#TO DO
#TO DO 
#TO DO

###alpha diversity ----

#prepare a matrix with OTU.norm.alpha diversity = "invsimpson"
OTU.norm.alpha = cbind(diversity(OTU.norm, index= "invsimpson"),design.keep)
colnames(OTU.norm.alpha)[1] = "alpha"


#linear mixed effect model on alpha diversity (block is random)
lmm.alpha <- lme(alpha~treatment+species+growing_stage,data = OTU.norm.alpha,random = ~1|bloc, method = "ML")
anova(lmm.alpha)

#              numDF denDF  F-value p-value
#(Intercept)       1   102 609.2532  <.0001
#treatment         1     8   0.0043  0.9495
#species           2     8  14.4030  0.0022 #only signif. effect
#growing_stage     1   102   0.7023  0.4040


#3 boxplots for treatment, species and early/late
dev.new(width=10, height=6,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,3),mar = c(6,5,4,2))

#treatment
boxplot(OTU.norm.alpha$alpha~OTU.norm.alpha$treatment,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Treatment",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(control),italic(inoculated)))

#species
boxplot(OTU.norm.alpha$alpha~OTU.norm.alpha$species,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Species",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(Wheat),italic(Corn),italic(Soy)))

#growing_stage
boxplot(OTU.norm.alpha$alpha~OTU.norm.alpha$growing_stage,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Growing stage",cex.lab = 1.5,cex.main = 1.5,names = expression(italic(early),italic(late)))

dev.print(device=pdf, "figures/figure3_alpha.pdf", onefile=FALSE)
dev.off()

###beta diversity ----
#standardization of raw data OR normalized data?
OTU.norm.hel <-decostand(OTU, "hel")

#permanova
(adonis(formula=OTU.norm.hel ~ species*treatment*growing_stage*bloc, data=design.keep, permutations=999, method="bray"))

#note that results are the same irrespective of prior normalization or not.
#                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#species                   2    8.2220  4.1110  47.493 0.43650  0.001 ***
#treatment                 1    0.0461  0.0461   0.532 0.00244  0.820    
#growing_stage             1    0.3383  0.3383   3.908 0.01796  0.003 ** 
#bloc                      8    1.2314  0.1539   1.778 0.06537  0.001 ***
#species:growing_stage     1    0.2425  0.2425   2.801 0.01287  0.011 *  
#treatment:growing_stage   1    0.0843  0.0843   0.974 0.00448  0.413    
#growing_stage:bloc        3    0.2753  0.0918   1.060 0.01462  0.368    
#Residuals                97    8.3963  0.0866         0.44575           
#Total                   114   18.8362                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Crop generate different microbial communities; growth stage generate different microbial community and those two factors interact together
##However, the pcoa reveals that the impact of growth stage could also only be due to the fact that soybeans have only been sampled at early stage. So we will re-run a permanova with only the wheat and the corn samples


####PcoA ----
#Calculating Bray-Curtis dissimilarity matrix on the hellinger transformed data
OTU.norm.hel.bray <-vegdist(OTU.norm.hel, method="bray")

#Calculating PCoA
OTU.norm.hel.bray.pcoa<-pcoa(dist(OTU.norm.hel.bray))

#How many axes represent more variability (42)
bs = OTU.norm.hel.bray.pcoa$values$Broken_stick
length(bs[bs>mean(bs)])

#PVE of first 2 axes (4.7% & 3.8%)
(OTU.norm.hel.bray.pcoa$values$Broken_stick/sum(OTU.norm.hel.bray.pcoa$values$Broken_stick))[1:2]

#Ploting the PCoAs - with inoculation as empty circles
#crops are "darkred","darkblue","darkorange
col = design.keep$species
col = gsub("Corn","darkred",col);col=gsub("Wheat","darkorange",col);col=gsub("Soy","darkblue",col)
dev.new()
plot(OTU.norm.hel.bray.pcoa$vectors[,1],OTU.norm.hel.bray.pcoa$vectors[,2],col = col, pch = ifelse(design.keep$treatment == "inoculated",19,21),
     ylab = "PC2", xlab = "PC1")
legend(1,1.5,fill = c("darkred","darkorange","darkblue"),legend = c("    Corn","    Wheat","    Soy"),box.lwd = 1)
legend(1.15,1.5,fill = rep("transparent",3), border = c("darkred","darkorange","darkblue"),legend = rep("",3),box.lwd = 0,box.col = "transparent")

#inoculation text (this is a pain...)
text(c(1.37,1.47),c(1.63,1.56),c("inoculated","control"),srt = 45,pos = 3,font =3)
dev.print(device=pdf, "figures/figure4_pcoa.pdf", onefile=FALSE)
dev.off()

###sandbox ----
#homoscédacity
if(1==2){
par(mfrow=c(2,2))
plot(anov_spores_inoculation_CE2)
## Graph Scale-location: the dispersion doesn't seem to increase that much, so homoscedacity might be respected
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
###Let's rerun anovas with log-transformed

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

}

