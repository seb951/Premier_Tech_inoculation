
### Repo description  
sebastien.renaut@gmail.com  
2018-2019  

This repo will describe the analyses and main results looking at the effect of Rhizoglomus irregulare inoculation on the microbial community of glomeromycetes in corn/soy and wheat.   

Experiment was done in collaboration with Premier Tech and led by Mohamed Hijri.   

Analyses were done by several, mainly Jacynthe Masse. Sebastien Renaut redid most analyses + graphs (getting Virtual taxa, stat analyses) in June 2018  

There are several scripts that describe the analyses:  

### Rcode
*PremierTech_MH_bioinfo.txt* #this describes the pipeline to get to an OTU table (it's kind of messy at the moment)  
*basic_sequencing_stats.R* #A few basic stats were calculated here on the raw, trimmed & merged sequencing data  
*file_description* #Check this to see were the temporary files used in the (*virtual_taxa.R*) and (*statistical_analyses.R*) originated  
*virtual_taxa.R* #Rcode to get from OTUs to VT (based on the MaarJam database)  
*statistical_analyses.R* #statistical analyses to generate results and figures.  

### Figures  
*figure1_VTX00113.pdf*  
*figure2_OTUabundance.pdf*  
*figure3_OTUabundance_per_genera_order_family.pdf*  
*figure4_alpha.pdf*  
*figure5_pcoa.pdf*  
*figure6_rda.pdf*  
*legends #figure legends*  

### Results  
Most results right now (linear mixed effect models + permanova) are actually in Rcode/stats.R  
*OTU.table.VT.merged.ordered.top.tsv*  
*PT_AMF.otu_table_from_biom.txt*  
*PT_combined_denoised_formated2.blastn.out*  
*PT_combined_denoised_formated2.fasta*  
*PT_crop_design.txt*  
*r_irregulares_16S.blastout* #results of the blast of all VT against the R. irregulare genome.  
*TableS1_OTU_VT_abundance.txt # a new table with all OTU names, VTX correspondance, phylogeny, and absolute abundance.
