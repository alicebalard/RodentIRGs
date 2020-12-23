library(ggplot2)
setwd("../")

my23posdf <- read.table("23BlastPosition_noseq.txt")
my23posdf

names(my23posdf) <- c("Pos", "Chr", "Start", "End")
my23posdf$Chr <- factor(my23posdf$Chr, levels = c("LR761981.1", "LR761995.1", "LR761979.1","LR761921.1"))

ggplot(my23posdf) +
  geom_rect(aes(xmin=Start, xmax=End, ymin = 0, ymax=1)) +
  facet_grid(.~Chr, scales = "free_x") +
  geom_label(aes(label=Pos, x = Start, y = rep(c(0.2, 0.5, 0.8),8))) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y= element_blank(),
        axis.title =  element_blank())




#################

  
blastn <- read.table("finalBlastn_IRGmus.vs.refmyodesglarolus.outfmt6")
tblastn <- read.table("finalTblastn_IRGmus.vs.refmyodesglarolus.outfmt6")

# add counter per hit
library(tidyverse)
blastn <- blastn %>% 
  group_by(V1) %>%
  mutate(hitnum = row_number()) %>% data.frame()
tblastn <- tblastn %>% 
  group_by(V1) %>%
  mutate(hitnum = row_number()) %>% data.frame()

##################### 
## Convert to fasta##
#####################

fasta_blastn <- paste(paste0(">", paste(blastn$V1,blastn$hitnum, sep = "_hit"),
                             paste0(" ",blastn$V2, ":", blastn$V9,"-", blastn$V10)), 
                      blastn$V13, sep = "\n")
write.table(fasta_blastn, file = "finalBlastn_IRGmus.vs.refmyodesglarolus.fasta", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

fasta_tblastn <- paste(paste0(">", paste(tblastn$V1,tblastn$hitnum, sep = "_hit"),
                             paste0(" ",tblastn$V2, ":", tblastn$V9,"-", tblastn$V10)), 
                       tblastn$V13, sep = "\n")
write.table(fasta_tblastn, file = "finalTblastn_IRGmus.vs.refmyodesglarolus.fasta", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##########
## Visu ##
##########

#### choose
blast <- blastn
# blast <- tblastn

# does it start by ATG?
blast$isATG = FALSE
blast[substr(blast$V13, 1,3) == "ATG","isATG"] <- TRUE

# forward or reverse?
blast$sens <- "F"
blast[blast$V9 >blast$V10,"sens"] <- "R"

### KEEP BEST 5 HITS
# blast <- blast[blast$hitnum %in% 1:5,]

######
## all
d <- blast

#####################
## first chromosome:
d <- blast[blast$V2 %in% "emb|LR761921.1|",]

#####################
## second chromosome:
d <- blast[blast$V2 %in% "emb|LR761981.1|",]

##################### Start
d <- data.frame(hitnum=d$hitnum, sens=d$sens, isATG=d$isATG, posy=1:nrow(d), chr = d$V2, Irg = d$V1,
                value=d$V9)

pas = 50000
ggplot(d, aes(x=value, y=posy, col = isATG, label = paste(Irg, hitnum))) + 
  geom_point(aes(pch = sens), size = 10) +
  ggtitle(label = as.character(d$chr)) +
  theme_linedraw() +
  scale_x_continuous(breaks = seq(ceiling(min(d$value)/pas) *pas,
                                  ceiling(max(d$value)/pas) *pas, pas),
                     labels = paste0(seq(ceiling(min(d$value)/pas) *pas,
                                         ceiling(max(d$value)/pas) *pas, pas)/1000, "K"))+ # every 50K
  scale_y_continuous(breaks = seq(0, nrow(d)))+
  geom_label()+
  theme(axis.text.y = element_blank(), axis.ticks.y= element_blank(),
        axis.title =  element_blank()) #+
    # facet_grid(chr~., scales = "free") # to uncomment for "all"

# list 10 candidates
giveMeHit <- function(x, y) blast[blast$V1 %in% x & blast$hitnum %in% y,]  

giveMeHit("Irga7",8)

# Hit1:a6.1
# Hit2:a7.5
# Hit3:a8.2
# Hit4:a7.2
# Hit5:a4.5 REV
# Hit6: 8.6
# Hit7: 7.1
# Hit8: 7.8 REV 
# Hit9: 7.4 INV!!
# Hit10: 3.1 REV

# plot the hits
dfIrga <- data.frame(labels = c("Irga1/2/6_mgl","Irga4/7.1_mgl","Irga4/7.2_mgl","Irga4/7.3_mgl",
                                "Irga4/7.4_mgl","Irga4/7.5_mgl", "Irga4/7.6_mgl","Irga4/7.7_mgl",
                                "Irga8/3.1_mgl","Irga8/3.2_mgl"),
                     start= c(21849150,21880821,21910811,21935793,22018870,22068242,
                              22100022,22167960,21895336,22002670),
                     end=c(21850387, 21882081, 21912058, 21937049, 22020117, 22069495,
                           22098755, 22169171, 21896600, 22003903),
                     chr="LR761921.1",
                     supportOF = c(T,T,T,F,F,F,T,T,F,T))

dfIrga$direction = "F"
dfIrga[dfIrga$start > dfIrga$end, "direction"] <- "R"

ggplot(dfIrga) + 
  ggtitle(label = as.character(dfIrga$chr)) +
  geom_rect(aes(xmin=start, xmax=end, ymin=1,ymax=2, col = supportOF))+
  geom_text(angle = 90, aes(x=start, y=3, label = labels)) +
  geom_text(aes(label = direction, x=start, y = 2.1))+
  theme_linedraw()+
  scale_y_continuous(limits = c(1,4)) +
  theme(axis.text.y = element_blank(), axis.ticks.y= element_blank(),
        axis.title =  element_blank()) +
  # every 50K
  scale_x_continuous(breaks = seq(ceiling(min(dfIrga$start)/1000) *1000, 
                                  ceiling(max(dfIrga$end)/1000) *1000, 50000),
                     labels = paste0(seq(ceiling(min(dfIrga$start)/1000) *1000, 
                                         ceiling(max(dfIrga$end)/1000) *1000, 50000)/1000, "K")) 

d <- blast[blast$V2 %in% "emb|LR761921.1|",]
d <- data.frame(hitnum=d$hitnum, sens=d$sens, isATG=d$isATG, posy=1:nrow(d), chr = d$V2, Irg = d$V1, 
                rbind(data.frame(value=d$V9, pos = "start"), data.frame(value=d$V10, pos = "end")))

library(ggplot2)
ggplot(d, aes(x=value, y=posy, col = isATG, label = paste(Irg, hitnum))) + 
  geom_point(aes(pch = sens), size = 10) +
  theme_linedraw() +
  scale_x_continuous(breaks = seq(ceiling(min(d$value)/1000) *1000, 
                                  ceiling(max(d$value)/1000) *1000, 50000),
                     labels = paste0(seq(ceiling(min(d$value)/1000) *1000, 
                                         ceiling(max(d$value)/1000) *1000, 50000)/1000, "K"))+ # every 50K
  scale_y_continuous(breaks = seq(0, nrow(d)))+
  geom_label()

# list 10 candidates
giveMeHit <- function(x, y) blast[blast$V1 %in% x & blast$hitnum %in% y,]  

giveMeHit("Irga7",8)

# Hit1:a6.1
# Hit2:a7.5
# Hit3:a8.2
# Hit4:a7.2
# Hit5:a4.5 REV
# Hit6: 8.6
# Hit7: 7.1
# Hit8: 7.8 REV 
# Hit9: 7.4 INV!!
# Hit10: 3.1 REV

# plot the hits
dfIrga <- data.frame(labels = c("Irga1/2/6_mgl","Irga4/7.1_mgl","Irga4/7.2_mgl","Irga4/7.3_mgl",
                                "Irga4/7.4_mgl","Irga4/7.5_mgl", "Irga4/7.6_mgl","Irga4/7.7_mgl",
                                "Irga8/3.1_mgl","Irga8/3.2_mgl"),
                     start= c(21849150,21880821,21910811,21935793,22018870,22068242,
                              22100022,22167960,21895336,22002670),
                     end=c(21850387, 21882081, 21912058, 21937049, 22020117, 22069495,
                           22098755, 22169171, 21896600, 22003903),
                     chr="LR761921.1",
                     supportOF = c(T,T,T,F,F,F,T,T,F,T))

dfIrga$direction = "F"
dfIrga[dfIrga$start > dfIrga$end, "direction"] <- "R"

ggplot(dfIrga) + 
  ggtitle(label = as.character(dfIrga$chr)) +
  geom_rect(aes(xmin=start, xmax=end, ymin=1,ymax=2, col = supportOF))+
  geom_text(angle = 90, aes(x=start, y=3, label = labels)) +
  geom_text(aes(label = direction, x=start, y = 2.1))+
  theme_linedraw()+
  scale_y_continuous(limits = c(1,4)) +
  theme(axis.text.y = element_blank(), axis.ticks.y= element_blank(),
        axis.title =  element_blank()) +
  # every 50K
  scale_x_continuous(breaks = seq(ceiling(min(dfIrga$start)/1000) *1000, 
                                  ceiling(max(dfIrga$end)/1000) *1000, 50000),
                     labels = paste0(seq(ceiling(min(dfIrga$start)/1000) *1000, 
                                         ceiling(max(dfIrga$end)/1000) *1000, 50000)/1000, "K")) 


##### Solving the Irga problem...

# First, let's match the hits with the orthofinder orthologues

Irgachro <- blast[blast$V2 %in% "emb|LR761921.1|",]

giveMeHit <- function(x, y, pos) {
  blastn = blast[blast$V1 %in% x & blast$hitnum %in% y,]  
  paste(paste0(">", pos, " ", paste(blastn$V1,blastn$hitnum, sep = "_hit"),
               paste0(" ",blastn$V2, ":", blastn$V9,"-", blastn$V10)), 
        blastn$V13, sep = "\n")
}
  
Pa.1 <- giveMeHit("Irga3",1, "P1")
Pa.2 <- giveMeHit("Irga4",3, "P2")
Pa.3 <- giveMeHit("Irga8",2, "P3")
Pa.4 <- giveMeHit("Irga4",1, "P4")
Pa.5 <- giveMeHit("Irga7",3, "P5")
Pa.6 <- giveMeHit("Irga8",6, "P6")
Pa.7 <- giveMeHit("Irga3",3, "P7")
Pa.8 <- giveMeHit("Irga3",5, "P8")
Pa.9 <- giveMeHit("Irga7",8, "P9")
Pa.10 <- giveMeHit("Irga1",2, "P10")
Pa.11 <- giveMeHit("Irga4",6, "P11")
Pa.12 <- giveMeHit("Irga3",1, "P12")

allAcand <- c(Pa.1,Pa.2,Pa.3,Pa.4,Pa.5,Pa.6,Pa.7,Pa.8,Pa.9,Pa.10,Pa.11, Pa.12)

write.table(allAcand, file = "allAcandidates.fasta", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


