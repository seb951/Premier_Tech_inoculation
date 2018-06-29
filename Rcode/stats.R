####setting things up ----
setwd("/Users/jerry/Documents/CSBQ/hijri/Premier_Tech")

#required packages
library(vegan)
library(nlme)
library(lme4)
library(ape)
library(pals)

###Preparing OTU table and sampling design----
OTU = read.table("results/OTU.table.VT.merged.ordered.top.tsv",row.names = 1, sep = "\t",header =T,stringsAsFactors = F)

#taxonomy in a separate vector
taxo = OTU[,ncol(OTU)]
OTU = data.frame(t(OTU[,-c(ncol(OTU))]))

#order based on rownames
OTU = OTU[order(rownames(OTU)),]

#normalize based on total number of OTU in sample (given that some seq libraries may have started with more sequences, it's probably good to do...)
OTU.norm = OTU/rowSums(OTU)

#load sampling design
design = read.table("results/PT_crop_design.txt", sep = "\t", stringsAsFactors = F, row.names = 1,header = T)

#order based on rownames
design = design[order(rownames(design)),]

#remove samples that were discarded in the design (would need to ask Jacynthe how this was done?)
keep = rownames(design) %in% rownames(OTU.norm)
design.keep = design[keep,]

###VTX00113 ----
#This is the most abundant VTX, and a Glomus spp,so...
#Is R. irregularis (VTX00113) more abundant in inoculated soils?

#prepare a matrix with only VTX00113 
OTU.norm.VTX00113 = cbind(sqrt(OTU.norm[,1]),design.keep)
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

shapiro.test(lmm1$residuals) #normaly distributed with the square root transform.


#boxplots
dev.new()
par(mfrow = c(2,2))
boxplot(OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Wheat",1]~OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Wheat",5],main = "VTX00113 - Wheat",ylab = "OTU relative abundance",names = expression(italic(control),italic(inoculated)))
boxplot(OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Corn",1]~OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Corn",5],main = "VTX00113 - Corn",ylab = "OTU relative abundance",names = expression(italic(control),italic(inoculated)))
boxplot(OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Soy",1]~OTU.norm.VTX00113[OTU.norm.VTX00113[,2]=="Soy",5],main = "VTX00113 - Soy",ylab = "OTU relative abundance",names = expression(italic(control),italic(inoculated)))
dev.print(device=pdf, "figures/figure1_VTX00113.pdf", onefile=FALSE)
dev.off()

#figure legends
system("echo 'Figure 1: Relative abundance of the most abundant Virtual Taxa (VTX00113 representing Rhizophagus irregularis, 33% of all reads) in all three species (wheat, Corn, Soy) tested.' >figures/legends")

###barplots of all VTX to look at diversity ----
#first group things according to treatment and species.
OTU.norm.barplot = cbind(colMeans(OTU.norm[design.keep[,1] == "Wheat" & design.keep[,4] == "inoculated", ]), 
      colMeans(OTU.norm[design.keep[,1] == "Wheat" & design.keep[,4] == "control", ]),
      colMeans(OTU.norm[design.keep[,1] == "Corn" & design.keep[,4] == "inoculated", ]), 
      colMeans(OTU.norm[design.keep[,1] == "Corn" & design.keep[,4] == "control", ]),
      colMeans(OTU.norm[design.keep[,1] == "Soy" & design.keep[,4] == "inoculated", ]), 
      colMeans(OTU.norm[design.keep[,1] == "Soy" & design.keep[,4] == "control", ]))
      
colnames(OTU.norm.barplot) = c("wheat.inoc", "wheat.ctl","corn.inoc","corn.ctl","soy.inoc","soy.ctl")      

#barplot
dev.new()
par(mar=c(16,4,4,2))
barplot(OTU.norm.barplot,beside = F,font = 3, axisnames = T,ylab = "Relative OTU abundance",col = c(cols25(),rep("grey",60)), las = 3, xpd = T,names.arg =  rep(c("Inoculated","Control"),3))
text(y = rep(1.08,6), x = c(1.3,3.6,6.2), labels = c("Wheat","Corn","Soy"), font = 2,cex = 2,xpd = T)
legend(0.2,-0.4,fill = cols25(),legend = paste(rownames(OTU.norm.barplot)[1:10],taxo[1:10],sep = "_"),cex = 0.75,xpd =T)
dev.print(device=pdf, "figures/figure2_OTUabundance.pdf", onefile=FALSE)
dev.off()

#figure legends
system("echo 'Figure 2: mean Relative abundance of the Virtual Taxa per treatment and species' >>figures/legends")

##barplot by genera / order / family ----
#TO DO, but is this necessary?
#TO DO 
#TO DO

###alpha diversity ----
#prepare a matrix with alpha diversity as "invsimpson" index
OTU.norm.alpha = cbind(diversity(OTU.norm, index= "invsimpson"),design.keep)
OTU.norm.alpha[,1] = log(OTU.norm.alpha[,1])
colnames(OTU.norm.alpha)[1] = "alpha"

#linear mixed effect model on alpha diversity (block is random)
lmm.alpha <- lme(alpha~treatment+species+growing_stage,data = OTU.norm.alpha,random = ~1|bloc, method = "ML")
anova(lmm.alpha)

#              numDF denDF  F-value p-value
#(Intercept)       1   102 609.2532  <.0001
#treatment         1     8   0.0043  0.9495
#species           2     8  14.4030  0.0022 #again, only signif. effect
#growing_stage     1   102   0.7023  0.4040

shapiro.test(lmm.alpha$residuals) #normaly distributed with the log transform

#boxplot: 3 boxplots for treatment, species and early/late
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

#figure legends
system("echo 'Figure 3: alpha diversity per Treatment / Species and Growing stage' >>figures/legends")


###beta diversity ----
#raw data OR normalized data?
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).

#standardization using hellinger transform
OTU.norm.hel <-decostand(OTU, "hel")

#PERMANOVA (I only keep the species*growing_stage interaction as the other interactions were all NS)
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
adonis(formula=OTU.norm.hel ~ species*growing_stage+treatment+bloc, data=design.keep, permutations=999, method="bray")

#                          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  species                 2    8.2220  4.1110  47.421 0.43650  0.001 ***
#  growing_stage           1    0.3395  0.3395   3.916 0.01802  0.003 ** 
#  treatment               1    0.0448  0.0448   0.517 0.00238  0.854    
#  bloc                    8    1.2314  0.1539   1.776 0.06537  0.003 ** 
#  species:growing_stage   1    0.2425  0.2425   2.797 0.01287  0.016 *  
#  Residuals             101    8.7559  0.0867         0.46485           
#  Total                 114   18.8362                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Species generate VERY different microbial communities (r2=0.47)
#Growth stage generate SMALL (r2=0.02) different microbial community and those two factors interact together (r2=0.01)
#Bloc has a SMALL (r2=0.06), but sign. effect too? Not actually sure it should be included or how. Maybe as STRATA, but that feels odd too.
#However, treatment has NO EFFECT, in any scenario tested.
#However, the pcoa reveals that the impact of growth stage could also only be due to the fact that soybeans have only been sampled at early stage,
#So we will re-run a permanova with only the wheat and the corn samples

####beta diversity per species to look at growth stage?----
#TO DO
#TO DO
#TO DO

####PcoA ----
#Calculating Bray-Curtis dissimilarity matrix on the hellinger transformed data
OTU.norm.hel.bray <-vegdist(OTU.norm.hel, method="bray")

#Calculating PCoA
OTU.norm.hel.bray.pcoa<-pcoa(dist(OTU.norm.hel.bray))

#How many axes represent more variability (42)
bs = OTU.norm.hel.bray.pcoa$values$Broken_stick
length(bs[bs>mean(bs)])

#PVE of first 2 axes (4.7% & 3.8%)
axis.1.2 = round((OTU.norm.hel.bray.pcoa$values$Broken_stick/sum(OTU.norm.hel.bray.pcoa$values$Broken_stick))[1:2],4)*100

#Ploting the PCoAs - with inoculation as empty circles
#crops are "darkred","darkblue","darkorange
col = design.keep$species
col = gsub("Corn","darkred",col);col=gsub("Wheat","darkorange",col);col=gsub("Soy","darkblue",col)
dev.new()
plot(OTU.norm.hel.bray.pcoa$vectors[,1],OTU.norm.hel.bray.pcoa$vectors[,2],col = col, pch = ifelse(design.keep$treatment == "inoculated",19,21),
     ylab = paste("PC2 (",axis.1.2[1],"%)",sep = ""), xlab = paste("PC1 (",axis.1.2[1],"%)",sep = ""))
legend(1,1.5,fill = c("darkred","darkorange","darkblue"),legend = c("    Corn","    Wheat","    Soy"),box.lwd = 1)
legend(1.15,1.5,fill = rep("transparent",3), border = c("darkred","darkorange","darkblue"),legend = rep("",3),box.lwd = 0,box.col = "transparent")

#inoculation text (this is a pain...)
text(c(1.37,1.47),c(1.63,1.56),c("inoculated","control"),srt = 45,pos = 3,font =3)
dev.print(device=pdf, "figures/figure4_pcoa.pdf", onefile=FALSE)
dev.off()

#figure legends
system("echo 'Figure 4: PcoA of the all samples colored coded according to species (wheat, corn, soy). Empty / full circles represent the treatment effect' >>figures/legends")


###sandbox ----

#Assumption 1: Normal distribution of the residues
#shapiro.test(resid(anov_spores_inoculation_CE))

##Assumption 2: Homoegeneity of the variances (homoscedacity)
#don't even think this is required in a permanova or a linear mixed effect model...
#bartlett.test(spores ~ Inoculation, data=colonization_corn_early)

#Assumption 3: indepence of samples.


