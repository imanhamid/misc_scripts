#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))

tped <- read.table("/Users/iman/Desktop/highfreq_all-indiv_missing0.tped", header=FALSE)

names(tped) <- c("chr", "pos_ID", "cM", "pos", "kpa003_1", "kpa003_2",
                 "kpa004_1", "kpa004_2", "kpa009_1", "kpa009_2", "kpa014_1", 
                 "kpa014_2", "kpa018_1", "kpa018_2")

tped2 <- tped
tped2[fix,]$chr <- tped2[fix,]$chr-1
tped2$pos_ID2 <- paste0(tped2$chr, "_", tped2$pos)

pos_IDs <- read.table("/Users/iman/Desktop/all_variantIDs.txt", header=FALSE)

names(pos_IDs) <- c("chr", "pos", "rsID", "REF", "ALT")

pos_IDs$pos_ID <- paste0(pos_IDs$chr, "_", pos_IDs$pos)
pos_IDs <- pos_IDs[order(pos_IDs$chr, pos_IDs$pos),]

tped_matched <- tped2[match(pos_IDs$pos_ID, tped2$pos_ID2),]

tped_matched$pos_ID <- pos_IDs$rsID

tped_matched2 <- tped_matched

tped_matched2$kpa003_1 <- ifelse(tped_matched2$kpa003_1==pos_IDs$REF, tped_matched2$kpa003_1, ifelse(tped_matched2$kpa003_1==pos_IDs$ALT, tped_matched2$kpa003_1, 0))
tped_matched2$kpa003_2 <- ifelse(tped_matched2$kpa003_2==pos_IDs$REF, tped_matched2$kpa003_2, ifelse(tped_matched2$kpa003_2==pos_IDs$ALT, tped_matched2$kpa003_2, 0))
tped_matched2$kpa004_1 <- ifelse(tped_matched2$kpa004_1==pos_IDs$REF, tped_matched2$kpa004_1, ifelse(tped_matched2$kpa004_1==pos_IDs$ALT, tped_matched2$kpa004_1, 0))
tped_matched2$kpa004_2 <- ifelse(tped_matched2$kpa004_2==pos_IDs$REF, tped_matched2$kpa004_2, ifelse(tped_matched2$kpa004_2==pos_IDs$ALT, tped_matched2$kpa004_2, 0))
tped_matched2$kpa009_1 <- ifelse(tped_matched2$kpa009_1==pos_IDs$REF, tped_matched2$kpa009_1, ifelse(tped_matched2$kpa009_1==pos_IDs$ALT, tped_matched2$kpa009_1, 0))
tped_matched2$kpa009_2 <- ifelse(tped_matched2$kpa009_2==pos_IDs$REF, tped_matched2$kpa009_2, ifelse(tped_matched2$kpa009_2==pos_IDs$ALT, tped_matched2$kpa009_2, 0))
tped_matched2$kpa014_1 <- ifelse(tped_matched2$kpa014_1==pos_IDs$REF, tped_matched2$kpa014_1, ifelse(tped_matched2$kpa014_1==pos_IDs$ALT, tped_matched2$kpa014_1, 0))
tped_matched2$kpa014_2 <- ifelse(tped_matched2$kpa014_2==pos_IDs$REF, tped_matched2$kpa014_2, ifelse(tped_matched2$kpa014_2==pos_IDs$ALT, tped_matched2$kpa014_2, 0))
tped_matched2$kpa018_1 <- ifelse(tped_matched2$kpa018_1==pos_IDs$REF, tped_matched2$kpa018_1, ifelse(tped_matched2$kpa018_1==pos_IDs$ALT, tped_matched2$kpa018_1, 0))
tped_matched2$kpa018_2 <- ifelse(tped_matched2$kpa018_2==pos_IDs$REF, tped_matched2$kpa018_2, ifelse(tped_matched2$kpa018_2==pos_IDs$ALT, tped_matched2$kpa018_2, 0))

tped_matched3 <- tped_matched2[-which(is.na(tped_matched2$pos)), -15]

write.table(tped_matched3, "/Users/iman/Desktop/highfreq_all-indiv_missing0_transversions.tped", quote = FALSE, row.names = FALSE, sep = "\t")

length(which(tped_matched3$kpa014_1!=0 & tped_matched3$kpa018_1!=0))
