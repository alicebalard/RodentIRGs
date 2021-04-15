# A. Balard
# 14 April 2021

library(dplyr)

dataOF <- read.csv("../data/data4syntheny_OF.csv") 

# rm pseudogenes (no protein product)
dataOF <- dataOF[!dataOF$IrgGroup %in% "pseudogene",]

# for each IRGgroup and Species, number the candidates
dataOF <- dataOF %>% group_by(Species, IrgName) %>% mutate(id = row_number())

# rm locations without protein name on NCBI (e.g. M.glareolus)
dataOF <- dataOF[!is.na(dataOF$protein_name),]

# make smaller species name
dataOF$species <- gsub("\\.", "", substring(dataOF$Species, 1,5))

# add myProtId to each protein name
protIdTable <- data.frame(myProtId = paste(dataOF$species, 
                                           paste(dataOF$IrgName, dataOF$id, sep = "."), sep = "_"),
                          NCBIprotName = dataOF$protein_name)

# export table
write.table(protIdTable, "../data/protIdTable_OF.csv", row.names = F, sep=",",  col.names=FALSE, quote = F)

# Follow up on bash to fetch sequences on NCBI:
# extractFastaSequences_part2.sh