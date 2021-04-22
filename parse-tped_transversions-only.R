#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))

tped <- read.table("HO_kopila_missing0_transversions-REDO.tped", header=FALSE)
names(tped) <- c("chr", "pos_ID", "cM", "pos", "kpa003_1", "kpa003_2",
                 "kpa004_1", "kpa004_2", "kpa009_1", "kpa009_2", "kpa014_1", 
                 "kpa014_2", "kpa018_1", "kpa018_2")
tped <- tped[which(tped$chr!="X" & tped$chr!="Y"),]
pos_IDs <- read.table("../HumanOrigins/AADR/v44.3_HO_public.snp", header=FALSE)
names(pos_IDs) <- c("rsID", "chr", "cM", "pos", "ref", "alt")
pos_IDs <- pos_IDs[which(pos_IDs$chr!=23 & pos_IDs$chr!=24),]
pos_IDs_transversions <- pos_IDs[which((pos_IDs$ref=="A" & pos_IDs$alt=="G") | (pos_IDs$ref=="G" & pos_IDs$alt=="A") | (pos_IDs$ref=="C" & pos_IDs$alt=="T") | (pos_IDs$ref=="T" & pos_IDs$alt=="C")),]
pos_IDs_transversions$chr_pos <- paste0(pos_IDs_transversions$chr, "_", pos_IDs_transversions$pos)
tped_matched <- tped[match(pos_IDs_transversions$chr_pos, tped$pos_ID),]
tped_matched$rsID <- pos_IDs_transversions$rsID
tped_matched$ref <- pos_IDs_transversions$ref
tped_matched$alt <- pos_IDs_transversions$alt
tped_matched2 <- tped_matched[which(!is.na(tped_matched$pos)),]


tped_matched2$kpa003_1 <- ifelse(tped_matched2$kpa003_1==tped_matched2$ref, tped_matched2$kpa003_1, ifelse(tped_matched2$kpa003_1==tped_matched2$alt, tped_matched2$kpa003_1, 0))
tped_matched2$kpa003_2 <- ifelse(tped_matched2$kpa003_2==tped_matched2$ref, tped_matched2$kpa003_2, ifelse(tped_matched2$kpa003_2==tped_matched2$alt, tped_matched2$kpa003_2, 0))
tped_matched2$kpa004_1 <- ifelse(tped_matched2$kpa004_1==tped_matched2$ref, tped_matched2$kpa004_1, ifelse(tped_matched2$kpa004_1==tped_matched2$alt, tped_matched2$kpa004_1, 0))
tped_matched2$kpa004_2 <- ifelse(tped_matched2$kpa004_2==tped_matched2$ref, tped_matched2$kpa004_2, ifelse(tped_matched2$kpa004_2==tped_matched2$alt, tped_matched2$kpa004_2, 0))
tped_matched2$kpa009_1 <- ifelse(tped_matched2$kpa009_1==tped_matched2$ref, tped_matched2$kpa009_1, ifelse(tped_matched2$kpa009_1==tped_matched2$alt, tped_matched2$kpa009_1, 0))
tped_matched2$kpa009_2 <- ifelse(tped_matched2$kpa009_2==tped_matched2$ref, tped_matched2$kpa009_2, ifelse(tped_matched2$kpa009_2==tped_matched2$alt, tped_matched2$kpa009_2, 0))
tped_matched2$kpa014_1 <- ifelse(tped_matched2$kpa014_1==tped_matched2$ref, tped_matched2$kpa014_1, ifelse(tped_matched2$kpa014_1==tped_matched2$alt, tped_matched2$kpa014_1, 0))
tped_matched2$kpa014_2 <- ifelse(tped_matched2$kpa014_2==tped_matched2$ref, tped_matched2$kpa014_2, ifelse(tped_matched2$kpa014_2==tped_matched2$alt, tped_matched2$kpa014_2, 0))
tped_matched2$kpa018_1 <- ifelse(tped_matched2$kpa018_1==tped_matched2$ref, tped_matched2$kpa018_1, ifelse(tped_matched2$kpa018_1==tped_matched2$alt, tped_matched2$kpa018_1, 0))
tped_matched2$kpa018_2 <- ifelse(tped_matched2$kpa018_2==tped_matched2$ref, tped_matched2$kpa018_2, ifelse(tped_matched2$kpa018_2==tped_matched2$alt, tped_matched2$kpa018_2, 0))


tped_matched3 <- tped_matched2[which(tped_matched2$kpa003_1!=0 | tped_matched2$kpa004_1!=0 | tped_matched2$kpa009_1!=0 | tped_matched2$kpa014_1!=0 | tped_matched2$kpa018_1!=0),]
tped_final <- tped_matched3[,c(1,15, 3:14)]
write.table(tped_final, "HO_kopila_missing0_transversions_REDO.tped", append=FALSE, quote = FALSE, row.names = FALSE, sep = "\t", col.names=FALSE)
