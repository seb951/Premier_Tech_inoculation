####setting things up ----
setwd("/Users/jerry/Documents/CSBQ/hijri/Premier_Tech")

#file describing the bioinfo analyses done by Jacynthe Masse.
#/Users/costanza/Desktop/04_MohammedHijri/01_MiSeq_PremierTech/01_Jacynthe/Pipeline_usedlast/PremierTech_MH_bioinfo.txt

#fasta file after the OTU were identified, but before things were clustered at 97%.
#/Users/costanza/Desktop/04_MohammedHijri/01_MiSeq_PremierTech/01_Jacynthe/Pipeline_usedlast/Former/PT_combined_denoised_formated2.fasta

#OTU table after OTU were identified, but before clustering (file was identified at L475 in script PT_AMF.otu_table_from_biom.txt labelled above)
#/Users/costanza/Desktop/04_MohammedHijri/01_MiSeq_PremierTech/01_Jacynthe/Pipeline_usedlast/Former/PT_AMF.otu_table_from_biom.txt

#Sampling design
#/Users/costanza/Desktop/04_MohammedHijri/01_PremierTech_analyses/02_Anayses/PT_crop.design.xlsx

#R script from Jacynthe
#/Users/costanza/Desktop/04_MohammedHijri/01_PremierTech_analyses/02_Anayses/PT_AMF_R/Script_PT_AMF_analyses.R


###blast / getting the Virtual Taxa ----
#local Blastn against Maarjam glomeromycetes.fasta database (a little slow on a single CPU, but not bad)
#set a pretty stringent evalue (1e-50), otherwise, will over-cluster...
blast_cmd = paste("blastn -query results/PT_combined_denoised_formated2.fasta -db blast_database/glomeromycetes.fasta -max_target_seqs 10 -evalue 1e-50 -outfmt 6 -out results/PT_combined_denoised_formated2.blastn.out")
system(blast_cmd)

#load blast results
blastn = read.table("results/PT_combined_denoised_formated2.blastn.out",header = F, stringsAsFactors = F)

#load OTU names and remove it
system("grep '>' results/PT_combined_denoised_formated2.fasta >results/OTU.names")
OTU.names = read.table("results/OTU_names",stringsAsFactors = F, header = F)
OTU.names = cbind(gsub(">","",OTU.names[,1]),0)
colnames(OTU.names) = c("OTU","VT")
system("rm results/OTU_names")

#new blast table with only best annotated hit (VTX) per OTU.
for(i in 1:nrow(OTU.names))
{
  temp = 0
  temp = blastn[blastn[,1] == OTU.names[i,1],]
  
  if(length(temp)>1){
  #ordered according to evalue
  temp.ordered = temp[order(temp[,11]),]
  
  #only keep annotations with a VTX
  temp.ordered.vtx = temp.ordered[regexpr("VTX",temp.ordered[,2])>0,]
  
  #keep best hit
  split = strsplit(temp.ordered.vtx[1,2],"_")[[1]]
  OTU.names[i,2] = split[length(split)]
  }
}

#merge the OTU table according to the VT
OTU.table = read.table("results/PT_AMF.otu_table_from_biom.txt", sep = "\t", stringsAsFactors = F, header = T,skip=1,comment.char = "")
OTU.table = OTU.table[,-ncol(OTU.table)]
OTU.table$VT = 0

###add VT info to OTU table (only if present, otherwise, it stays with its 'OTU name')
for(i in 1:nrow(OTU.table))
  {
    temp = OTU.names[OTU.names[,1] == OTU.table[i,1],2]
    if(temp !=0) OTU.table$VT[i] = temp #add VT 
  #  if(temp == 0) OTU.table$VT[i] = OTU.table[i,1] #add OTU name, because no VT.
  }

###create new OTU table with VT merged ----
vt = rle(sort(OTU.table$VT))$values
vt = vt[length(vt):1]
OTU.table.VT.merged = data.frame(matrix(0, nrow= length(vt),ncol = (ncol(OTU.table)-1)))
colnames(OTU.table.VT.merged) = c("vt",colnames(OTU.table)[2:(ncol(OTU.table)-1)])
  
for(i in 1:nrow(OTU.table.VT.merged))
  {
    temp = OTU.table[OTU.table$VT == vt[i],-c(1,ncol(OTU.table))]
    OTU.table.VT.merged[i,1] = vt[i]
    OTU.table.VT.merged[i,2:ncol(OTU.table.VT.merged)] = colSums(temp)
 }

#order dataframe and keep only the top "percent" % expressed
percent = 1
OTU.table.VT.merged.ordered = OTU.table.VT.merged[order(rowSums(OTU.table.VT.merged[,2:117]),decreasing=T),]
OTU.table.VT.merged.ordered.top = OTU.table.VT.merged.ordered[1:round(nrow(OTU.table.VT.merged.ordered)*percent),]

#fraction of the reads it represents... (99.5% for the top 50%)
sum(rowSums(OTU.table.VT.merged.ordered[,2:117])[1:round(nrow(OTU.table.VT.merged.ordered)*percent)]) / sum(rowSums(OTU.table.VT.merged.ordered[,2:117]))

OTU.table.VT.merged.ordered.top$annotation = 0

#add annotation to the last line
for(i in 1:nrow(OTU.table.VT.merged.ordered.top))
{
  grep = paste("grep '",OTU.table.VT.merged.ordered.top[i,1],"' reference_material/vtx_taxonomy.txt >temp",sep ="")
  system(grep)
  
  temp = read.table("temp",stringsAsFactors = F, header = F)
  OTU.table.VT.merged.ordered.top$annotation[i] = paste(temp[,4:8],collapse = ";")
}
system("rm temp")

#get a proper rowname and colnames
colnames(OTU.table.VT.merged.ordered.top) = gsub(".","_",colnames(OTU.table.VT.merged.ordered.top),fixed = T)
rownames(OTU.table.VT.merged.ordered.top) = OTU.table.VT.merged.ordered.top[,1]
OTU.table.VT.merged.ordered.top = OTU.table.VT.merged.ordered.top[,-1]

#save the new merged OTU table
write.table(OTU.table.VT.merged.ordered.top,"results/OTU.table.VT.merged.ordered.top.tsv",sep = "\t",col.names = T,row.names = T)

