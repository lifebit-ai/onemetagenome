#!/usr/bin/env Rscript

library(plyr)

df <- read.csv("queryLca.csv")

unique_species <- df[!duplicated(df$LCA.NCBI.Taxon.ID),]

#count number of reads per species
reads_number <- count(df$LCA.NCBI.Taxon.ID)

#data.frame(reads_number)

colnames(reads_number)[1] <- "LCA.NCBI.Taxon.ID"

combined <- merge(unique_species,reads_number,by="LCA.NCBI.Taxon.ID")

total_reads_number <- read.table("reads_number.txt", quote="\"", comment.char="")
total_reads_number <- total_reads_number[[1]][1]

Percentage.of.Reads <- combined$freq/total_reads_number * 100

combined <- cbind(combined, Percentage.of.Reads)

#re-order columns
combined <- combined[c(4,3,1,5,6,2)]

colnames(combined) <- c("LCA Scientific Name","LCA Rank Name","LCA NCBI Taxon ID","No. of Reads","Percentage of Total Reads","Query Accession")
write.csv(combined,file="table.csv", row.names=FALSE)
