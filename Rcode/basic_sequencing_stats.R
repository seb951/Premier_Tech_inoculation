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


