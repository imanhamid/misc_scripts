#!/usr/bin/env Rscript

chr <- commandArgs(trailingOnly = TRUE)[1]
chr_length <- commandArgs(trailingOnly = TRUE)[2]

infile <- paste0("/work/ih49/1000GP_Phase3/genetic_map_", chr, "_combined_b37.txt")
outfile <- paste0("/work/ih49/1000GP_Phase3/", chr, "_recombinationrates.txt")

chr_map <- read.table(infile, header=TRUE)

chr_map$ends <- c(chr_map$position[-1] - 2,chr_length)

chr_map$rates <- c(chr_map$COMBINED_rate.cM.Mb.[-length(chr_map$COMBINED_rate.cM.Mb.)]*1e-8,0)
adjusted_map <- chr_map[,c("ends","rates")]

write.table(adjusted_map, file=outfile, append=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
