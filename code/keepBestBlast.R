# A. Balard
# 13 April 2021

library(dplyr)

filenames <- list.files("../data/blast_results/IRGmus2otherRodents/row_tblasn", full.names=TRUE)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

Results_blast_best <- data.frame()
Results_blast_full <- data.frame()

for (f in filenames){
  blastdata <- read.delim(file = f, header = F)
  
  # remove identical sequences search b5/b5*, b6/b6*, and   
  # keep only "alternative" sequences for Irgb10 and Irgq (corresponds to NCBI)
  blastdata <- blastdata[!blastdata$V1 %in% c("Irgb5*", "Irgb6*", "Irgb10", "Irgq"),]
  blastdata$V1[blastdata$V1 %in% "Irgb10alternative(longer)"] <- "Irgb10"
  blastdata$V1[blastdata$V1 %in% "Irgqalternative(longer)"] <- "Irgq"
  
  # remove "alternative" tandems (= concatenated protein) TO CHECK
  blastdata <- blastdata[!grepl("alternative", blastdata$V1),]
  
  blastHeaders <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sseq")
  colnames(blastdata) <- blastHeaders
  
  # forward or reverse
  blastdata$forward_reverse <- "F"
  blastdata[blastdata$sstart > blastdata$send,"forward_reverse"] <- "R"
  
  # write in good direction
  blastdata$start <- blastdata$sstart
  blastdata$start[blastdata$forward_reverse %in% "R"] <- blastdata$send[blastdata$forward_reverse %in% "R"]
  blastdata$end <- blastdata$send
  blastdata$end[blastdata$forward_reverse %in% "R"] <- blastdata$sstart[blastdata$forward_reverse %in% "R"]
  
  # find all candidates: find overlapping sequences
  require(GenomicRanges)
  fullblastRange <- GRanges(blastdata$sseqid, IRanges(blastdata$start, blastdata$end))
  candidateDF <- reduce(GRanges(blastdata$sseqid, IRanges(blastdata$start, blastdata$end)))
  whichCandidate <- data.frame(findOverlaps(fullblastRange, candidateDF, type="any"))
  blastdata$candidate <- whichCandidate$subjectHits
  
  blastdata$Species <- substrRight(f, 4)
  
  
  # for each candidate IRG, keep the row with the lowest p-value, then highest bitscore, then highest percentage identity
  blastdataFinal <- blastdata %>% group_by(candidate) %>% 
    dplyr::slice_max(order_by = bitscore) %>% # NB if by eval, keep the first randomly if eval = 0
    data.frame()
  
  Results_blast_full <- rbind(Results_blast_full, blastdata)
  Results_blast_best <- rbind(Results_blast_best, blastdataFinal)
}

# homogenize for synteny map: 
## (1) chr, start, end, forward_reverse colums
names(Results_blast_best)[names(Results_blast_best)%in% "sseqid"] <- "chr"
## (2)IrgGroup in "pseudogene", "IrgbTandem", "Irgb6-OG1", "Irgb6-OG2", "Irgb6-tbn","Irgb10",
#  "Irga2-6", "Irga4", "Irga3", "Irgm1", "Irgc", "Irgm2", "Irgd", "Irgq", "Irgm3
Results_blast_best$IrgGroup <- Results_blast_best$qseqid

unique(Results_blast_best$qseqid)
Results_blast_best$IrgGroup[Results_blast_best$IrgGroup %in% c("Irga1","Irga5","Irga8","Irga7", "Irgb7")] <- "pseudogene"
Results_blast_best$IrgGroup[Results_blast_best$IrgGroup %in% c("Irga2","Irga6")] <- "Irga2-6"
Results_blast_best$IrgGroup[
  Results_blast_best$IrgGroup %in% 
    c("Irgb9","Irgb8", "Irgb5*b3alternative", "Irgb8b9alternative", "Irgb1",
      "Irgb5*","Irgb5", "Irgb4", "Irgb3", "Irgb2", "Irgb2b1alternative")] <- "IrgbTandem"
Results_blast_best$IrgGroup[Results_blast_best$IrgGroup %in% "Irgb10alternative(longer)"] <- "Irgb10"
Results_blast_best$IrgGroup[Results_blast_best$IrgGroup %in% "Irgb6*"] <- "Irgb6"
Results_blast_best$IrgGroup[Results_blast_best$IrgGroup %in% "Irgqalternative(longer)"] <- "Irgq"

## (3) IrgName
Results_blast_best <- data.frame(Results_blast_best %>% group_by(Species, qseqid) %>% mutate(id = row_number()))
Results_blast_best$IrgName <- paste(Results_blast_best$qseqid, Results_blast_best$id, sep = ".")

## (4) Species names
Results_blast_best$Species[Results_blast_best$Species %in% "Aamp"] <- "A.amphibius"
Results_blast_best$Species[Results_blast_best$Species %in% "Asyl"] <- "A.sylvaticus"
Results_blast_best$Species[Results_blast_best$Species %in% "Marv"] <- "M.arvalis"
Results_blast_best$Species[Results_blast_best$Species %in% "Maur"] <- "M.auratus"
Results_blast_best$Species[Results_blast_best$Species %in% "Mgla"] <- "M.glareolus"
Results_blast_best$Species[Results_blast_best$Species %in% "Mgl2"] <- "M.glareolus genome Bank_vole2_10x"
Results_blast_best$Species[Results_blast_best$Species %in% "Mmus"] <- "M.musculus"
Results_blast_best$Species[Results_blast_best$Species %in% "Moch"] <- "M.ochrogaster"
Results_blast_best$Species[Results_blast_best$Species %in% "Mpah"] <- "M.pahari"
Results_blast_best$Species[Results_blast_best$Species %in% "Pleu"] <- "P.leucopus"
Results_blast_best$Species[Results_blast_best$Species %in% "Rnor"] <- "R.norvegicus"
Results_blast_best$Species[Results_blast_best$Species %in% "Rrat"] <- "R.rattus"

## WRITE OUT
write.csv(Results_blast_best, "../data/data4syntheny_tblastn.csv", row.names = F)

#test <- Results_blast_full[Results_blast_full$qseqid %in% c("Irgb8", "Irgb4"),]
