# A. Balard
# October 2021
# Analyses done on Harriet server from Heitlinger team

library(tidyr)
library(rentrez) # to retrieve NCBI sequences
library(tidyverse)

## We match the OrthoFinder results with the known IRG names in Mmus
tabRef <- read.csv("/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/Table1_IRGs_naming.csv")

tabOF <- read.csv("/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/orthofinder_results/N0_HOG_IRG.tsv", sep="\t")

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

##############################
## SAVING POINT 
save(data_FINAL, file = "/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/data_FINAL.RData")
## to load
#load("/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/data_FINAL.RData")
##############################

data_FINAL  <- data_FINAL[!data_FINAL$protein_seq_only %in% "Supplied id parameter is empty.",]

nrow(data_FINAL)# 237

### Make fasta file from orthofinder results
fastaOF  <- paste0(">", data_FINAL$IRGname," ", data_FINAL$protein_name, "\n", data_FINAL$protein_seq_only)

write.table(fastaOF, file = "/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/candidatesFromOF.fasta", sep = "\n", col.names=FALSE, row.names = FALSE, quote=FALSE )
