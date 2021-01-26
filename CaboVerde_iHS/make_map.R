#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
  
chr <- commandArgs(trailingOnly=TRUE)[1]

legend_file <- paste0("/work/ih49/CV_localancestry/iHS_calculation/CV_chr", chr, ".legend")
map_file <- paste0("/work/ih49/CV_localancestry/iHS_calculation/genetic_map_chr", chr, "_combined_b37.txt")
outfile <- paste0("/work/ih49/CV_localancestry/iHS_calculation/chr", chr, ".map")

chr_pos <- read.table(legend_file, header=FALSE)
chr_genetic_map <- read.table(map_file, header=TRUE)

get_map <- function(x) {
  genetic_map_pos <- chr_genetic_map[max(which(chr_pos$V2[x]>=chr_genetic_map$position)),3]
  map$Genetic_Map.cM.[x] <<- genetic_map_pos
}

map <- chr_genetic_map[match(chr_pos$V2, chr_genetic_map$position),]
map$position2 <- chr_pos$V2

indices <- which(is.na(map$Genetic_Map.cM.))

sapply(X=indices, FUN=get_map)

map$chr <- chr
map$locus_ID <- chr_pos$V1

map <- map[,c("chr", "locus_ID", "Genetic_Map.cM.", "position2")]

write.table(map, file=outfile, sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
