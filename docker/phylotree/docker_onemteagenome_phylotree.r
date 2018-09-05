#!/usr/bin/env Rscript
library(taxize)

queryLca <- read.delim("queryLca.tsv", header=FALSE)

#subset df to only include rows for species
species <- queryLca[(queryLca$V3=="species"),]
#remove duplicates
unique_species <- species[!duplicated(species$V2), ]

species_id <- as.character(unique_species[,2])
species_names <- as.character(unique_species[,4])

taxize_species <- classification(species_id, db = "ncbi")
taxize_tree <- class2tree(taxize_species, check = TRUE)

pdf("phylotree.pdf")
png("phylotree.png")
plot(taxize_tree)
dev.off()
