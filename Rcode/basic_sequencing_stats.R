#!/usr/bin/env Rscript

#Note that for now, data are sitting on my imac desktop. I will transfer it once manuscript is accepted.


#How many raw reads:
system("wc -l ~/Desktop/04_MohammedHijri/01_MiSeq_PremierTech/01_Jacynthe/02_Bioinformatic/*fastq >/Users/jerry/Documents/CSBQ/hijri/Premier_Tech/raw.counts")
raw = read.table("~/Desktop/raw.counts"); 1.96*sd(raw[-nrow(raw),1])/sqrt(nrow(raw)-1)

#How many trimmed
system("wc -l ~/Desktop/04_MohammedHijri/01_MiSeq_PremierTech/01_Jacynthe/02_Bioinformatic/*R?_paired.fq >/Users/jerry/Documents/CSBQ/hijri/Premier_Tech/trimmed.counts")
trim = read.table("~/Desktop/trimmed.counts"); 1.96*sd(trim[-nrow(trim),1])/sqrt(nrow(trim)-1)

#Getting trimmed/merged sequence lengths.
system("grep -v '>' ~/Desktop/04_MohammedHijri/01_MiSeq_PremierTech/01_Jacynthe/02_Bioinformatic/PT_AMF.trim.contigs.good.fasta | awk '{print length}' >/Users/jerry/Documents/CSBQ/hijri/Premier_Tech/sum")
sum = read.table("~/Desktop/sum"); 1.96*sd(sum[nrow(sum),1])/sqrt(nrow(sum))

#which Virtual Taxa in the one that was inoculated? (BLASTN)
#prerequisite: Download the Rhizophagus irregularis DAOM 181602 genome from NCBI. 
#download and format the VT taxa fasta file (https://maarjam.botany.ut.ee/ganja/?do=downloadFile&id=27)

system("blastn -query /Users/jerry/Documents/CSBQ/hijri/Premier_Tech/reference_material/vtx_taxonomy.fasta -db /Users/jerry/Documents/CSBQ/hijri/Premier_Tech/blast_database/GCA_002897155.1_Rir_HGAP_ii_V2_genomic.fa -evalue 1e-100 -outfmt 6 | awk '$3>99' >/Users/jerry/Documents/CSBQ/hijri/Premier_Tech/results/r_irregulares_16S.blastout")
