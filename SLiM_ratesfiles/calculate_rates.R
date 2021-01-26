AAmap <- read.table("AAmap.chr1.txt", header=TRUE)

#calculating cM/Mb by dividing the difference in cM divided by difference in Mb
#between current and next row

cM <- (AAmap$Genetic_Map.cM.[-1] - AAmap$Genetic_Map.cM.)

Mb <- ((AAmap$Physical_Position_Build36.hg18.[-1] - AAmap$Physical_Position_Build36.hg18.)/1e6)

AAmap$recombinationrate <- cM/Mb