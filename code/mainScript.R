# A. Balard
# October 2021
# Analyses done on Harriet server from Heitlinger team

library(tidyr)
library(rentrez) # to retrieve NCBI sequences

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
# fastaSequences <- lapply(data_long[1:499,"protein_name"], function(x)
#     entrez_fetch(db = "protein", rettype = 'fasta', id = x))
# fastaSequences2 <- lapply(data_long[500:999,"protein_name"], function(x)
#     entrez_fetch(db = "protein", rettype = 'fasta', id = x))
# fastaSequences3 <- lapply(data_long[1000:nrow(data_long),"protein_name"], function(x)  entrez_fetch(db = "protein", rettype = 'fasta', id = x))

# fastaSequencesAll <- unlist(c(fastaSequences, fastaSequences2, fastaSequences3))

# check: all goood
length(fastaSequencesAll)
nrow(data_long)

which(data_long$protein_name %in% "")
which(fastaSequencesAll %in% "Supplied id parameter is empty.\n")

data_long$protein_seq  <- fastaSequencesAll

## Get gene ID from protein ID (to remove duplicates, i.e. proteins from the same gene)
# genes  <- lapply(data_long[,"protein_name"], function(x){
#     if (x==""){
#         geneID <- ""
#     } else {    
#         gene_links <- entrez_link(dbfrom='protein', id=x, db='gene')
#         geneID <- gene_links$links[["protein_gene"]]
#     }
#     return(geneID)
# })
                                        #

## add to data frame

data_long$geneID  <- genes

## remove duplicates based on (1) similar gene ID and (2) similar protein sequence

# Extract characters after pattern to obtain only the protein sequence without header
data_long$protein_seq_only <- gsub("\n", "", gsub(".*]","", data_long$protein_seq))

## RM duplicates based on sequence and gene ID
data_final <- data_long[!duplicated(data_long[c("geneID","protein_seq_only")]),]

nrow(data_final) #269 candidates

## obtain chromosomic position
#     entrez_fetch(db = "protein", rettype = 'fasta', id = x))
# fastaSequences3 <- lapply(data_long[1000:nrow(data_long),"protein_name"], function(x)  entrez_fetch(db = "protein", rettype = 'fasta', id = x))

data_final$geneID[1]

test=entrez_fetch(db="gene", rettype="XML", data_final$geneID[1])

test=XML::xmlToList(test)

str(test)

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

t1 <- data.frame(t(sapply(data_final$geneID[c(1,2)], getPos)))

t2 <- data_final[c(1,2),]

t3 <- cbind(t1,t2)

t3$chr

## Add a counter for the IRG name by species
table(data_final$Name_Bekpen_2005,data_final$species)

data_final[data_final$species %in% "Mmus_GCF_000001635.27_GRCm39_protein" & data_final$Name_Bekpen_2005 %in% "Irgm1",]
# 2 isoforms

