#!/usr/bin/env Rscript

#How many raw reads:
system("wc -l ~/Desktop/04_MohammedHijri/01_MiSeq_PremierTech/01_Jacynthe/02_Bioinformatic/*fastq >~/Desktop/raw.counts")
raw = read.table("~/Desktop/raw.counts"); 1.96*sd(raw[-nrow(raw),1])/sqrt(nrow(raw)-1)

#How many trimmed
system("wc -l ~/Desktop/04_MohammedHijri/01_MiSeq_PremierTech/01_Jacynthe/02_Bioinformatic/*R?_paired.fq >~/Desktop/trimmed.counts")
trim = read.table("~/Desktop/trimmed.counts"); 1.96*sd(trim[-nrow(trim),1])/sqrt(nrow(trim)-1)

#Getting trimmed/merged sequence lengths.
system("grep -v '>' PT_AMF.trim.contigs.good.fasta | awk '{print length}' >sum")
sum = read.table("~/Desktop/sum"); 1.96*sd(sum[nrow(sum),1])/sqrt(nrow(sum))

#which Virtual Taxa in the one that was inoculated? (BLASTN)
#prerequisite: Download the Rhizophagus irregularis DAOM 181602 genome from NCBI. 
#download and format the VT taxa fasta file (https://maarjam.botany.ut.ee/ganja/?do=downloadFile&id=27)

system("blastn -query ../reference_material/vtx_taxonomy.fasta -db ../blast_database/GCA_002897155.1_Rir_HGAP_ii_V2_genomic.fa -evalue 1e-100 -outfmt 6 | awk '$3>99' >../results/r_irregulares_16S.blastout")
