# A. Balard
# 13 April 2021

library(dplyr)
library(seqinr)

# pb in files

filenames1 <- list.files("../data/blast_results/IRGmus2otherRodents/", full.names=TRUE, pattern = "\\.outfmt6$")
filenames2 <- list.files("../data/blast_results/IRGmus2transcriptomes/", full.names=TRUE, pattern = "\\.outfmt6$")

makedata4syntheny <- function(filenames){
  
  Results_blast_best <- data.frame()
  Results_blast_full <- data.frame()
  
  for (f in filenames){
    blastdata <- read.delim(file = f, header = F)
    
    # remove identical sequences search b6/b6*, and   
    blastdata <- blastdata[!blastdata$V1 %in%  "Irgb6*",]
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
    
    blastdata$Species <- substr(sub(".*tBlastnIRGmusGRCm39_vs_", "", f), 1,4)
    
    # for each candidate IRG, keep the row with the highest bitscore (if by eval, keep the first randomly if eval = 0)
    blastdataFinal <- blastdata %>% group_by(candidate) %>% 
      dplyr::slice_max(order_by = bitscore) %>% 
      data.frame()
    
    Results_blast_full <- rbind(Results_blast_full, blastdata)
    Results_blast_best <- rbind(Results_blast_best, blastdataFinal)
  }
  
  # homogenize for synteny map:
  ## (1) chr, start, end, forward_reverse colums
  names(Results_blast_best)[names(Results_blast_best)%in% "sseqid"] <- "chr"
  ## (2)IrgGroup in "pseudogene", "IrgbTandem", "Irgb6-OG1", "Irgb6-OG2", "Irgb6-tbn","Irgb10",
  #  "Irga2-6", "Irga4", "Irga3", "Irgm1", "Irgc", "Irgm2", "Irgd", "Irgq", "Irgm3
  Results_blast_best$IrgGroup <- gsub("_refMmus", "", gsub("_refMus", "", Results_blast_best$qseqid))
  Results_blast_best$IrgGroup[grepl("Irga1|Irga5|Irga8|Irga7|Irgb7", Results_blast_best$IrgGroup)] <- "pseudogene"
  Results_blast_best$IrgGroup[grepl("Irga2|Irga6", Results_blast_best$IrgGroup)] <- "Irga2-6"
  Results_blast_best$IrgGroup[Results_blast_best$IrgGroup %in%
                                c("Irgb2","Irgb1","Irgb4","Irgb5","Irgb3","Irgb9","Irgb8")] <- "IrgbTandem"
  Results_blast_best$IrgGroup[Results_blast_best$IrgGroup %in% "Irgb6*"] <- "Irgb6"

  ## (3) IrgName
  Results_blast_best <- data.frame(Results_blast_best %>% group_by(Species, qseqid) %>% mutate(id = row_number()))
  Results_blast_best$IrgName <- gsub("_refMus", "", paste(Results_blast_best$qseqid, Results_blast_best$id, sep = "."))
  
  # (4) Species names
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
  # for transcriptomics
  Results_blast_best$Species[Results_blast_best$Species %in% "AAL-"] <- "A.agrarius"
  Results_blast_best$Species[Results_blast_best$Species %in% "BVK-"] <- "M.glareolus"
  Results_blast_best$Species[Results_blast_best$Species %in% "FMN-"] <- "M.arvalis"
  
  return(Results_blast_best)
}

## WRITE OUT both fasta and file for synteny
write.csv(makedata4syntheny(filenames1), "../data/data4syntheny_tblastn.csv", row.names = F)
## make fasta
d <- makedata4syntheny(filenames1) ; dlist <- as.list(d$sseq)
names(dlist) <- paste0(d$IrgName, "_", d$Species, " (tblastn Mmus vs genome)")
write.fasta(sequences = dlist, file.out = "../data/blast_results/IRGmus2otherRodents/tBlastnIRGmusGRCm39_vs_allRodents_besthits.fasta", 
            names = names(dlist), nbchar = 10000)

write.csv(makedata4syntheny(filenames2), "../data/data4syntheny_tblastn-transcriptome.csv", row.names = F)
## make fasta
d <- makedata4syntheny(filenames2) ; dlist <- as.list(d$sseq)
names(dlist) <- paste0(d$IrgName, "_", d$Species, " (tblastn Mmus vs de novo transcriptome)")
write.fasta(sequences = dlist, file.out = "../data/blast_results/IRGmus2transcriptomes/tBlastnIRGmusGRCm39_vs_denovotranscriptomes_besthits.fasta", 
            names = names(dlist), nbchar = 10000)
  