####setting things up ----
setwd("/Users/jerry/Documents/CSBQ/hijri/Premier_Tech")

#required packages
library(vegan)
library(nlme)
library(lme4)
library(ape)
library(pals)
library(dplyr)
library(ggplot2)
library(RVAideMemoire)

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

#Seven samples were discarded
rownames(design)[keep == F]

###VTX00113 ----
#This is the most abundant VTX, and a Glomus spp,so this is really what we are interested in...
#Is R. irregulare (VTX00113) more abundant in inoculated soils?

#prepare a matrix with only VTX00113 (note that you can also test VTX00114 (OTU.norm[,33]) and it shows the same pattern as VTX00113)
OTU.norm.VTX00113 = cbind(sqrt(OTU.norm[,1]),design.keep)
colnames(OTU.norm.VTX00113)[1] = colnames(OTU.norm)[1]

#linear mixed effect model (block is random)
lmm1 <- lme(VTX00113~treatment+species+growing_stage,data = OTU.norm.VTX00113,random = ~1|bloc, method = "REML")
anova(lmm1)

#              numDF denDF  F-value p-value
#(Intercept)       1   102 728.8265  <.0001
#treatment         1     8   0.2599  0.6240
#species           2     8  39.2675  0.0001***sig effect
#growing_stage     1   102   4.0835  0.0459** sig
#
shapiro.test(lmm1$residuals) #normaly distributed with the square root transform.

#test the interaction in a linear model?
lmm2 <- lm(VTX00113~treatment+species*growing_stage,data = OTU.norm.VTX00113)
anova(lmm2)

#                       Df Sum Sq Mean Sq F value   Pr(>F)    
#treatment               1 0.0139 0.01391  0.4142  0.52119    
#species                 2 2.9888 1.49440 44.4875 7.49e-15 ***
#growing_stage           1 0.1398 0.13978  4.1612  0.04378 *  
#species:growing_stage   1 0.1289 0.12885  3.8358  0.05272 .  *not sign.. !!!!


#boxplots
dev.new()
#boxplot: 3 boxplots for treatment, species and early/late
dev.new(width=10, height=6,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,3),mar = c(6,5,4,2))

#treatment
boxplot(OTU.norm.VTX00113$VTX00113~OTU.norm.VTX00113$treatment,beside = F,font = 3, axisnames = T,ylab = "VTX00113 relative abundance",las = 3, xpd = T,main = "Treatment",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(control),italic(inoculated)))

#species
boxplot(OTU.norm.VTX00113$VTX00113~OTU.norm.VTX00113$species,beside = F,font = 3, axisnames = T,ylab = "VTX00113 relative abundance",las = 3, xpd = T,main = "Species",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(Corn),italic(Soy),italic(Wheat)))

#growing_stage
boxplot(OTU.norm.VTX00113$VTX00113~OTU.norm.VTX00113$growing_stage,beside = F,font = 3, axisnames = T,ylab = "VTX00113 relative abundance",las = 3, xpd = T,main = "Growing stage",cex.lab = 1.5,cex.main = 1.5,names = expression(italic(early),italic(late)))
dev.print(device=pdf, "figures/figure1_VTX00113.pdf", onefile=FALSE)
dev.off()

#figure legends
system("echo 'Figure 1: Relative abundance of the most abundant Virtual Taxa (VTX00113 represents the R. irregulare DAOM-197198 inoculant, 38% of all reads) in all three species (wheat, Corn, Soy) tested. Species and Growing stage had a significant effect on the presence of VTX00113' >figures/legends")



###barplots of all VTX to look at diversity among species and inoculation----
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
barplot(OTU.norm.barplot,beside = F,font = 3, border = 0, axisnames = T,ylab = "Relative OTU abundance",col = c(cols25(),rep("grey",60)), las = 3, xpd = T,names.arg =  rep(c("Inoculated","Control"),3),space = 0.05)
text(y = rep(1.08,6), x = c(1.1,3.2,5.2), labels = c("Wheat","Corn","Soy"), font = 2,cex = 2,xpd = T)
legend(0.1,-0.4,border = "transparent",fill = cols25(),legend =  paste(rownames(OTU.norm.barplot)[1:10],rep(";",10),taxo[1:10],c("*",rep("",9)),sep = ""),cex = 0.75,xpd =T)
dev.print(device=pdf, "figures/figure2_OTUabundance.pdf", onefile=FALSE)
dev.off()

#print figure caption
system("echo 'Figure 2: Mean relative abundance of the Virtual Taxa and taxonomic information according to the MaarjAM (Opik et al. 2014) database per treatment and species. *VTX00113 represents the R. irregulare DAOM-197198 inoculant.' >>figures/legends")

###barplot by genera / order / family ----
#get taxonomy info.
OTU.norm.barplot.taxo = as.data.frame(OTU.norm.barplot)
OTU.norm.barplot.taxo$genera = unlist(strsplit(taxo,";"))[seq(4,(length(taxo)*5),by = 5)]
OTU.norm.barplot.taxo$order = unlist(strsplit(taxo,";"))[seq(3,(length(taxo)*5),by = 5)]
OTU.norm.barplot.taxo$family = unlist(strsplit(taxo,";"))[seq(2,(length(taxo)*5),by = 5)]

#summarize with dplyr.
OTU.norm.barplot.taxo_summary.GENERA = OTU.norm.barplot.taxo[,-c(8,9)] %>% group_by(genera) %>% summarise_all(sum)
OTU.norm.barplot.taxo_summary.ORDER = OTU.norm.barplot.taxo[,-c(7,9)] %>% group_by(order) %>% summarise_all(sum)
OTU.norm.barplot.taxo_summary.FAMILY = OTU.norm.barplot.taxo[,-c(7,8)] %>% group_by(family) %>% summarise_all(sum)

#barplots according to genera/order/family
dev.new()
par(mar= c(6,4,4,2))

#genera
g = barplot(as.matrix(as.data.frame(OTU.norm.barplot.taxo_summary.GENERA)[,2:7]),border = 0,
        ylab = "Relative OTU abundance",beside = F,names.arg = rep("",6),
        col = c(cols25()[1:6],rep("grey",60)),xlim = c(0,18*1.2),space = c(0.1))

#order
o = barplot(as.matrix(as.data.frame(OTU.norm.barplot.taxo_summary.ORDER)[,2:7]),names.arg = rep("",6),border = 0,
        add=T,beside = F,col = c(cols25()[7:11],rep("grey",60)),space = c(1.1*6.5,rep(0.1,5)))

#family
f = barplot(as.matrix(as.data.frame(OTU.norm.barplot.taxo_summary.FAMILY)[,2:7]),names.arg = rep("",6),border = 0,
        add=T,beside = F,col = c(cols25()[12:13],rep("grey",60)),space = c(1.1*13-0.1,rep(0.1,5)))

#add x axix labels
axis(1,at = c(g,o,f), labels = rep(colnames(OTU.norm.barplot.taxo_summary.GENERA)[2:7],3),font = 3, las = 3,cex.axis=0.8)

#add legends
legend(0.5,0.9,legend = OTU.norm.barplot.taxo_summary.GENERA$genera,fill = cols25()[1:6],cex = 0.7,bg = "white")
legend(7.3,0.9,legend = OTU.norm.barplot.taxo_summary.ORDER$order,fill = cols25()[7:11],cex = 0.7,bg = "white")
legend(14.5,0.9,legend = OTU.norm.barplot.taxo_summary.FAMILY$family,fill = cols25()[12:13],cex = 0.7,bg = "white")
text(y = rep(1.04,6), x = c(3.3,10.5,17.5), labels = c("Genera","Order","Family"), font = 2,cex = 2,xpd = T)
dev.print(device=pdf, "figures/figure3_OTUabundance_per_genera_order_family.pdf", onefile=FALSE)
dev.off()

#print figure caption
system("echo 'Figure 3: Mean relative abundance of the Virtual Taxa grouped by genera, order and family according to the MaarjAM (Opik et al. 2014) database per treatment and species.' >>figures/legends")

###alpha diversity ----
#prepare a matrix with alpha diversity as "invsimpson" index
OTU.norm.alpha = cbind(diversity(OTU.norm, index= "invsimpson"),design.keep)
OTU.norm.alpha[,1] = log(OTU.norm.alpha[,1])
colnames(OTU.norm.alpha)[1] = "alpha"

#linear mixed effect model on alpha diversity (block is random)
lmm.alpha <- lme(alpha~treatment+species+growing_stage,data = OTU.norm.alpha,random = ~1|bloc, method = "REML")
anova(lmm.alpha)

#              numDF denDF  F-value p-value
#(Intercept)       1   102 864.2705  <.0001
#treatment         1     8   0.0030  0.9575
#species           2     8  12.7100  0.0033 ***again only sig effect
#growing_stage     1   102   2.4562  0.1202

shapiro.test(lmm.alpha$residuals) #normaly distributed with the log transform

#test the interaction in a linear model?
lmm.alpha2 <- lm(alpha~treatment+species*growing_stage,data = OTU.norm.alpha)
anova(lmm.alpha2)

#                       Df  Sum Sq Mean Sq F value    Pr(>F)    
#treatment               1  0.0006 0.00061  0.0030    0.9561    
#species                 2  5.1225 2.56127 12.7627 1.047e-05 ***
#growing_stage           1  0.4950 0.49496  2.4664    0.1192    
#species:growing_stage   1  0.2922 0.29224  1.4562    0.2301   # not signif!  
# Residuals             109 21.8746 0.20068 

###alpha diversity:boxplot ----
#3 boxplots for treatment, species and early/late
dev.new(width=10, height=6,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,3),mar = c(6,5,4,2))

#treatment
boxplot(OTU.norm.alpha$alpha~OTU.norm.alpha$treatment,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Treatment",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(control),italic(inoculated)))

#species
boxplot(OTU.norm.alpha$alpha~OTU.norm.alpha$species,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Species",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(Corn),italic(Soy),italic(Wheat)))

#growing_stage
boxplot(OTU.norm.alpha$alpha~OTU.norm.alpha$growing_stage,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Growing stage",cex.lab = 1.5,cex.main = 1.5,names = expression(italic(early),italic(late)))

dev.print(device=pdf, "figures/figure4_alpha.pdf", onefile=FALSE)
dev.off()

#figure legends
system("echo 'Figure 4: alpha diversity per treatment, species and growing stage. Species and growing stage had a significant effect on alpha diversity.' >>figures/legends")



###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).

#standardization using hellinger transform
OTU.norm.hel <-decostand(OTU, "hel")

#PERMANOVA (I only keep the species*growing_stage interaction as the other interactions were all NS)
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
adonis(formula=OTU.norm.hel ~ species*growing_stage+treatment+bloc, data=design.keep, permutations=9999, method="bray")

#                          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#species                 2    8.2220  4.1110  47.421 0.43650 0.0001 ***
#growing_stage           1    0.3395  0.3395   3.916 0.01802 0.0014 ** 
#treatment               1    0.0448  0.0448   0.517 0.00238 0.8526    
#bloc                    8    1.2314  0.1539   1.776 0.06537 0.0017 ** 
#species:growing_stage   1    0.2425  0.2425   2.797 0.01287 0.0127 *  
#Residuals             101    8.7559  0.0867         0.46485           
#Total                 114   18.8362                 1.00000           
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Species generate VERY different microbial communities (r2=0.47)
#Growth stage generate SMALL (r2=0.02) different microbial community and those two factors interact together (r2=0.01)
#Bloc has a SMALL (r2=0.06), but sign. effect too? Not actually sure it should be included or how. Maybe as STRATA, but that feels odd too.
#However, treatment has NO EFFECT, in any scenario tested.
#However, the pcoa reveals that the impact of growth stage could also only be due to the fact that soybeans have only been sampled at early stage,
#So we will re-run a permanova with only the wheat and the corn samples

#### posthoc permanova to look at species + growing stage effect----
#beta diversity per species to look at growth stage?----

#use species AND growing stage as factors (9999 perms for more precision given that one p is close to 0.05) 
pairwise.perm.manova(dist(OTU,method = "euclidean"),fact = paste(design.keep[,1],design.keep[,2],sep = "_"),nperm = 9999)

####PcoA: good to illustrate the community structure----
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
pch.pcoa = rep(0,nrow(design.keep))
pch.pcoa[design.keep$treatment == "inoculated"] = 19; pch.pcoa[design.keep$treatment == "control"] = 21

dev.new()
plot(OTU.norm.hel.bray.pcoa$vectors[,1],OTU.norm.hel.bray.pcoa$vectors[,2],col = col, pch = pch.pcoa,
     ylab = paste("PC2 (3.79%)",sep = ""), xlab = paste("PC1 (4.66%)",sep = ""))
abline(h = 0,lty=3,cex =0.7);abline(v=0,lty=3,cex =0.7)
legend(1,1.5,fill = c("darkred","darkorange","darkblue"),legend = c("    Corn","    Wheat","    Soy"),box.lwd = 1)
legend(1.15,1.5,fill = rep("transparent",3), border = c("darkred","darkorange","darkblue"),legend = rep("",3),box.lwd = 0,box.col = "transparent")

#inoculation text (this is a pain...)
text(c(1.37,1.47),c(1.63,1.56),c("inoculated","control"),srt = 45,pos = 3,font =3)
dev.print(device=pdf, "figures/figure5_pcoa.pdf", onefile=FALSE)
dev.off()

#figure legends
system("echo 'Figure 5: Principle Coordinate Analysis (PCoA) of the all samples colored coded according to species (wheat, corn, soy). Empty / full circles represent the treatment effect.' >>figures/legends")


###RDA ----
rownames(OTU.norm.barplot)[1:10]
rda = rda(formula=OTU.norm.hel ~ species*growing_stage, data=design.keep,scale = T)
anova.cca(rda, step=1000) #sign of model 
anova.cca(rda, by="axis", step=1000) #axis
anova.cca(rda, by="terms", step=1000) #variables

#Df Variance      F Pr(>F)    
#Model      9   20.114 4.9005  0.001 ***
#Residual 105   47.886   

#Df Variance       F Pr(>F)    
#RDA1       1    9.792 21.4709  0.001 ***
#RDA2       1    6.120 13.4186  0.001 ***
#RDA3       1    1.344  2.9472  0.006 **   


#Note that I tested also an RDA with species*growing_stage*treatment but treatment is NS, so I removed it.
#Df Variance       F Pr(>F)    
#species                           2   15.801 17.3232  0.001 ***
#growing_stage                     1    0.953  2.0888  0.006 ** 
#species:growing_stage             1    0.882  1.9349  0.018 * 

dev.new()

#empty plot, and keep plotting data structure to reuse with points()
rda = rda(formula=OTU.norm.hel ~ species*growing_stage, data=design.keep,scale = T)
rda.plot = plot(rda, type = "n", scaling = 3,ylab = "RDA1 (9.0%)",xlab ="RDA2 (14.3%)",xlim = c(-1.5,1.5))
col = design.keep$species
col = gsub("Corn","#8B000099",col);col=gsub("Wheat","#FF8C0099",col);col=gsub("Soy","#00008B99",col)

#plot the samples
points(rda.plot$sites,col = col,pch = 21)

#plot the VTX but NOT the top 10.
points(rda.plot$species[-c(1:10),],pch = 4,col = "darkcyan",lwd=1.2)

#plot the centroids
points(rda.plot$centroids,pch = 3, col=c("#8B0000","#00008B","#FF8C00","black","black"),lwd=2)
text(rda.plot$centroids,labels=c("Corn","Soy","Wheat","Late","Early"),col =c("darkred","darkblue","darkorange","black","black"), cex = 1,font=2,pos = 4)

#plot the VTX top 10 (with labels)
text(rda.plot$species[c(1:10),],labels = rownames(rda.plot$species[c(1:10),]),cex = 0.7, col = "darkcyan",font = 2,adj = 0.8)

#figure legends
system("echo 'Redundancy Analysis (RDA). Note that only the top10 most abundant VT (\"X\" signs) were labeled, \"O\" sign are the samples themselves and \"+\" sign represent the centroids of the species (wheat, corn and soy) and growing stage effect (early vs. late). Inoculation treatment is not labeled as samples did not group according to the treatment (e.g. see Figure 5).
       X and Y labels represent the first two Principal Component and the percentage of variance explained")

dev.print(device=pdf, "figures/figure6_rda.pdf", onefile=FALSE)
dev.off()

###sandbox ----

#font = as.numeric(ifelse(regexpr("Glomus",taxo)>0,"2","1"))
#cex = as.numeric(ifelse(regexpr("Glomus",taxo)>0,"0.6","0.4"))

#Assumption 1: Normal distribution of the residues
#shapiro.test(resid(anov_spores_inoculation_CE))

##Assumption 2: Homoegeneity of the variances (homoscedacity)
#don't even think this is required in a permanova or a linear mixed effect model...
#bartlett.test(spores ~ Inoculation, data=colonization_corn_early)

#Assumption 3: indepence of samples.