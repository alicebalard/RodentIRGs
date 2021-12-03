# A. Balard
# October 2021
# Analyses done on Harriet server from Heitlinger team

library(tidyr)
library(rentrez) # to retrieve NCBI sequences
library(tidyverse)

##########################
## OrthoFinder analysis ##
##########################

## We match the OrthoFinder results with the known IRG names in Mmus
tabRef <- read.csv("../data/Table1_IRGs_naming.csv")

tabOF <- read.csv("../data/orthofinder_results/N0_HOG_IRG.tsv", sep="\t")

tabRef <- tabRef[c("Name_Bekpen_2005", "proteins_GRCm39")]

## rm empty rows: pseudo genes, or not found (Irgb4)
tabRef <- tabRef[!tabRef$proteins_GRCm39 %in% "",]
tabRef <- tabRef[!tabRef$proteins_GRCm39 %in% "?",]

tabRef$proteins_GRCm39  <- gsub(" ", "", tabRef$proteins_GRCm39)

## one protein name per line with amazing separate_rows()
tabRef <- data.frame(separate_rows(tabRef, proteins_GRCm39, sep = ";"))

tabOF <- data.frame(separate_rows(tabOF, Mmus_GCF_000001635.27_GRCm39_protein, sep=", "))

# rename a column before merging
tabOF$proteins_GRCm39 <- tabOF$Mmus_GCF_000001635.27_GRCm39_protein

# merge
tabMERGED <- merge(tabRef, tabOF)

nrow(tabMERGED)

## wide to long
gathercols <- names(tabMERGED)[6:19]
data_long <- gather(tabMERGED, species, protein_name, gathercols)

## separate rows (one line per protein)
data_long  <-  data.frame(separate_rows(data_long, protein_name, sep = ", "))

head(data_long,20)

# no orthologue found
data_long[data_long$protein_name %in% "",]

# retrieve fasta from NCBI with rentrez (cut in 3 cause bug)
fastaSequences <- lapply(data_long[1:499,"protein_name"], function(x)
  entrez_fetch(db = "protein", rettype = 'fasta', id = x))

fastaSequences2 <- lapply(data_long[500:999,"protein_name"], function(x)
  entrez_fetch(db = "protein", rettype = 'fasta', id = x))

fastaSequences3 <- lapply(data_long[1000:nrow(data_long),"protein_name"], function(x)  entrez_fetch(db = "protein", rettype = 'fasta', id = x))

fastaSequencesAll <- unlist(c(fastaSequences, fastaSequences2, fastaSequences3))

# check: all goood
length(fastaSequencesAll)
nrow(data_long)

which(data_long$protein_name %in% "")
which(fastaSequencesAll %in% "Supplied id parameter is empty.\n")

data_long$protein_seq  <- fastaSequencesAll

## Get gene ID from protein ID (to remove duplicates, i.e. proteins from the same gene)
genes  <- lapply(data_long[,"protein_name"], function(x){
  if (x==""){
    geneID <- ""
  } else {    
    gene_links <- entrez_link(dbfrom='protein', id=x, db='gene')
    geneID <- gene_links$links[["protein_gene"]]
  }
  return(geneID)
})

## add to data frame
data_long$geneID  <- genes

# load("/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/data_long.RData")

# check if all lines have a geneID
data_long[sapply(data_long$geneID, is.null),]
##NP_001129214.2
## has been updated for NP_001129214.3

## MANUAL CORRECTION
data_long$protein_name[data_long$protein_name %in% "NP_001129214.2"]  <- "NP_001129214.3"
gene_links <- entrez_link(dbfrom='protein', id="NP_001129214.3", db='gene')
geneID <- gene_links$links[["protein_gene"]]
data_long$geneID[data_long$protein_name %in% "NP_001129214.3"] <- geneID

## remove duplicates based on (1) similar gene ID and (2) similar protein sequence

# Extract characters after pattern to obtain only the protein sequence without header
data_long$protein_seq_only <- gsub("\n", "", gsub(".*]","", data_long$protein_seq))

## RM duplicates based on sequence and gene ID
data_final <- data_long[!duplicated(data_long[c("geneID","protein_seq_only")]),]

nrow(data_final) #269 candidates

## obtain chromosomic position
getPos <- function(geneid){
  x = entrez_fetch(db="gene", rettype="docsum", geneid)
  ## get chromosome accession number
  a = gsub("</ChrAccVer>.*", "", gsub(".*<ChrAccVer>","", x))
  ## get chr start
  b = gsub("</ChrStart>.*", "", gsub(".*<ChrStart>","", x))
  ## get chr stop
  c = gsub("</ChrStop>.*", "", gsub(".*<ChrStop>","", x))
  return(list(chr=a,start=b,stop=c))
}

## fetch position info
posDF <- data.frame(t(sapply(data_final$geneID, getPos)))

data_FINAL <- cbind(data_final,posDF)

data_FINAL[1,]

## Remove isoforms: exclude similar geneID
data_FINAL <- data_FINAL[!duplicated(data_FINAL$geneID),]

nrow(data_FINAL) # 238

## Add a counter for the IRG name by species
data_FINAL$number = 1
data_FINAL <- data_FINAL %>%
  group_by(species, Name_Bekpen_2005) %>%
  mutate(ticker = cumsum(number))

data_FINAL$IRGname <- paste0(substr(data_FINAL$species, 1, 4), "_", data_FINAL$Name_Bekpen_2005, ".", data_FINAL$ticker)

data_FINAL  <- data_FINAL[!data_FINAL$protein_seq_only %in% "Supplied id parameter is empty.",]

nrow(data_FINAL)# 237

## Rename NameBekpen for tandems (1 protein, 2 names, let's merge)
data_FINAL$Name_Bekpen_2005[data_FINAL$Name_Bekpen_2005 %in% "Irgb1"] <-  "Irgb1b2"
data_FINAL$Name_Bekpen_2005[data_FINAL$Name_Bekpen_2005 %in% "Irgb3"] <-  "Irgb5*b3"
data_FINAL$Name_Bekpen_2005[data_FINAL$Name_Bekpen_2005 %in% "Irgb9"] <-  "Irgb9b8"

##############################
## SAVING POINT 
save(data_FINAL, file = "../data/data_FINAL_OF.RData")

## to load
load("../data/data_FINAL_OF.RData")
##############################

## Check orthogroups
table(data_FINAL$OG, data_FINAL$Name_Bekpen_2005)


### Make fasta file from orthofinder results
fastaOF  <- paste0(">", data_FINAL$IRGname," ", data_FINAL$protein_name, "\n", data_FINAL$protein_seq_only)

write.table(fastaOF, file = "/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/fasta_sequences/candidatesFromOF.fasta", sep = "\n", col.names=FALSE, row.names = FALSE, quote=FALSE )

################################################
## Reciprocal Best Blast Hit: retrieve sequences
################################################

musIRG <- read.csv("../data/GRCm39IRGprotseq.txt")

## Rename NameBekpen for tandems (1 protein, 2 names, let's merge)
musIRG$Name_Bekpen_2005[musIRG$Name_Bekpen_2005 %in% "Irgb1"] <-  "Irgb1b2"
musIRG$Name_Bekpen_2005[musIRG$Name_Bekpen_2005 %in% "Irgb3"] <-  "Irgb5*b3"
musIRG$Name_Bekpen_2005[musIRG$Name_Bekpen_2005 %in% "Irgb9"] <-  "Irgb9b8"

## Retrieve the sequences from NCBI
musIRGfasta <- lapply(musIRG[,"proteins_GRCm39"], function(x)
  entrez_fetch(db = "protein", rettype = 'fasta', id = x))

## Extract the IRG names
namesseq <- paste0(">Mmus_", musIRG$Name_Bekpen_2005, " ",
                   gsub(">", "", gsub("\n", "", gsub("].*","", unlist(musIRGfasta)))),
                   "]")

# Extract characters after pattern to obtain only the protein sequence without header
seq <- gsub("\n", "", gsub(".*]","", unlist(musIRGfasta)))

head(as.vector(t(data.frame(namesseq, seq))))

write.table(as.vector(t(data.frame(namesseq, seq))), "/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/fasta_sequences/Mmus_IRG_protseq.fasta", quote = F, row.names = F, col.names = F)
##############################


## The RBBH itself is done in bash (mainScript.sh)


##############################
## Retrieve RBBH outfmt6 files
PATH = "../data/blast_results/RBBS_Irgmus2set2"

myfiles  <- list.files(PATH, full.names=T)
myfiles <- myfiles[grep("outfmt6.2", myfiles)]

myfiles1 <- myfiles[grep("vs_IRGmusGRCm39", myfiles)]
myfiles2 <- myfiles[grep("tBlastnIRGmus", myfiles)]

rodentsAsQuery  <- lapply(myfiles1, read.table)
musAsQuery <- lapply(myfiles2, read.table)

## rename list elements
names(rodentsAsQuery)  <- gsub("_G", "", gsub("tBlastn", "",stringr::str_extract(myfiles1, "tBlastn.{0,6}")))

names(musAsQuery)  <- gsub("_G", "", gsub("_vs_", "",stringr::str_extract(myfiles2, "_vs_.{0,6}")))

## rename all columns within list
namesBlastO6  <- c("qseqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sseq")

musAsQuery <-  lapply(musAsQuery, setNames, namesBlastO6)
rodentsAsQuery <-  lapply(rodentsAsQuery, setNames, namesBlastO6)

#############################
#### TO RM AFTER REAL RUN
rodentsAsQuery <- lapply(rodentsAsQuery, function(y){
  sub <- strsplit(y$sseqid, "_")
  y$sseqid <- paste0("Mmus_", sapply(sub, function(x) {
    name = x[[1]]
    ## Rename NameBekpen for tandems (1 protein, 2 names, let's merge)
    name[name %in% "Irgb1"] <-  "Irgb1b2"
    name[name %in% "Irgb3"] <-  "Irgb5*b3"
    name[name %in% "Irgb9"] <-  "Irgb9b8"
    return(name)
  }))
  return(y)
})

#############################

# Step 1: Keep those reads with RBBH for IRGname + genomic sequence, and with similar sequence end and start +/- 5bp
library(dplyr)

# if multiple IRG at the same place, keep the more likely (low e-value, high bitscore)
getRBBHdf <- function(i, range = -150:150){
  R  <- rodentsAsQuery[[i]]
  M <- musAsQuery[[i]]
  
  dfR <- data.frame(A=R$sseqid, B=R$qseqid, C=R$qstart, D=R$qend,
                    E=R$evalue, F=R$bitscore, seq = R$sseq)
  dfM <- data.frame(A=M$qseqid, B=M$sseqid, C=M$sstart, D=M$send,
                    E=M$evalue, F=M$bitscore, seq = M$sseq)
  
  ## amplify the first dataframe (to extend the possible start and end of 100bp, still it's the same position)
  amplifiedR <- data.frame(A=unlist(lapply(dfR$A, function(x) rep(x,length(range)))),
                           B=unlist(lapply(dfR$B, function(x) rep(x,length(range)))),
                           C=unlist(lapply(dfR$C, function(x) x+range)),
                           D=unlist(lapply(dfR$D, function(x) rep(x, length(range)))),
                           EaR=unlist(lapply(dfR$E, function(x) rep(x, length(range)))),
                           FaR=unlist(lapply(dfR$F, function(x) rep(x, length(range)))))
  amplifiedR <- data.frame(A=unlist(lapply(amplifiedR$A, function(x) rep(x,length(range)))),
                           B=unlist(lapply(amplifiedR$B, function(x) rep(x,length(range)))),
                           C=unlist(lapply(amplifiedR$C, function(x) rep(x, length(range)))),
                           D=unlist(lapply(amplifiedR$D, function(x) x+range)),
                           EaR=unlist(lapply(amplifiedR$EaR, function(x) rep(x,length(range)))),
                           FaR=unlist(lapply(amplifiedR$FaR, function(x) rep(x,length(range)))))
  
  ## keep the corresponding sequences with matching start and end between both blast
  df  <-  inner_join(amplifiedR , dfM, by=c("A", "B", "C", "D"))
  return(df)
}

dfList <- lapply(1:length(rodentsAsQuery), getRBBHdf)

## Step 2: remove candidates for several IRGs
## For each hit, if there is another hit with better score save that, 
## if not save the line on focus

# for each rodent species
newdfList <- lapply(dfList, function(x){
  # for each line
  newdf <- data.frame()
  for (i in 1:nrow(x)){
    t <- x[i,]
    subtab <- x[(x$B %in% t$B & 
                   x$C %in% (t$C + c(-150:150)) &
                   x$D %in% (t$D + c(-150:150))),]
    newdf <- rbind(newdf, 
                   subtab[subtab$EaR %in% min(subtab$EaR) & subtab$FaR %in% max(subtab$FaR),])
  }
  newdf <- unique(newdf)
  return(newdf)
})

# Check if exact same hit for 2 different IRG (e.g. Irgb6 and b6*) and merge
newdfList <- lapply(newdfList, function(x){
  x <- x[!duplicated(paste(x$B, x$C, x$D)),]
})

vecFas <- vector()

## prepare the fasta file format (including a counter for each IRG)
test <- newdfList[[2]]
## rename the list:
names(newdfList) <- names(rodentsAsQuery)

names(rodentsAsQuery)[]

for (i in 1:length(newdfList)){
  newdfList[[i]]$number = rep(1, nrow(newdfList[[i]]))
  newdfList[[i]]$species = rep(names(newdfList)[i], nrow(newdfList[[i]]))
  N <- newdfList[[i]] %>%
    group_by(A) %>%
    mutate(counter=cumsum(number))
  N$irgName <- paste0(N$species, "_", gsub("Mmus_", "", N$A), ".", N$number)
  vecFas <- c(vecFas, paste0(">", N$irgName, "\n", N$seq))
}

prefasta <- gsub("-", "", vecFas)
  
# Save fasta file
write.table(prefasta, "../data/fasta_sequences/candidatesFromRBBH.fasta", quote = F, col.names = F, row.names = F)
              