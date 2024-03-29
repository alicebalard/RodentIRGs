library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)

## load  intermediate data from OF preprocessed in R mainScript.R
load("../data/data_FINAL_OF.RData")

# IRGname, start, stop, chr, species
data_FINAL <- data.frame(apply(data_FINAL, 2 , unlist))
data_FINAL$start <- as.numeric(data_FINAL$start)
data_FINAL$stop <- as.numeric(data_FINAL$stop)
data_FINAL$chr <- as.factor(data_FINAL$chr)

# short name species (NB check 2 versions of Mgla for RBBH)
data_FINAL$sp <- substr(data_FINAL$species, 1, 4)

## Plot for all chromosomes starting on the left, at the same scale
chromoDF <- data_FINAL %>% group_by(chr, sp) %>%
  summarise(startChr = min(c(start, stop)),
            stopChr = max(c(start, stop)), .groups = "keep")

chromoDF$chrRight <- chromoDF$stopChr - chromoDF$startChr
chromoDF$chrLeft <- 0

## Find long breaks in chromosomes to put [...] break symbol 
# which chr have to be broken? 
hist(chromoDF$chrRight)

# per chromosome, what are the min and max of start and end?
rangeDF <- data_FINAL %>% group_by(chr) %>% summarise(minstart = min(start),
                                                      maxstart = max(start),
                                                      range = abs(maxstart - minstart)) %>%  data.frame()

ggplot(rangeDF, aes(y=chr, x=range)) + geom_point(alpha=.3) +
  geom_vline(xintercept = 1e6)

chromoDF <- merge(chromoDF, rangeDF[c("chr", "range")])

# is chr broken?
chromoDF$broken <- "no"
chromoDF$broken[chromoDF$range > 1e6] <- "yes"

# what is the position for genes on first part/unbroken chromosomes?
## prepare genes positions
geneDF <- data_FINAL
geneDF <- merge(geneDF, chromoDF)

# original position
geneDF$geneLeft <- as.numeric(apply(geneDF, 1, function(x) min(x[["start"]], x[["stop"]])))
geneDF$geneRight <- as.numeric(apply(geneDF, 1, function(x) max(x[["start"]], x[["stop"]])))

# new position to start all on the left
geneDF$geneLeft <- geneDF$geneLeft - geneDF$startChr
geneDF$geneRight <- geneDF$geneRight - geneDF$startChr 

# and for genes on part 2?
geneDF$whichpart <- "part1"
geneDF$whichpart[geneDF$geneLeft > 1e6] <- "part2"

# where do we cut?
breakStart <- geneDF[geneDF$whichpart %in% "part1",] %>% group_by(chr) %>% 
  summarise(stopPart1 = max(geneLeft, geneRight) + 1000)
breakEnd <- geneDF[geneDF$whichpart %in% "part2",] %>% group_by(chr) %>%
  summarise(startPart2 = min(geneLeft, geneRight) - 1000)
chrEnd <- geneDF[geneDF$whichpart %in% "part2",] %>% group_by(chr) %>% summarise(stopChr2 = max(geneLeft) + 1000)

geneDF <- merge(merge(merge(geneDF, breakStart, all.x = T), breakEnd, all.x = T), chrEnd, all.x = T)

geneDF$breaksize <- geneDF$startPart2 - geneDF$stopPart1 - 2000 # 2000 for borders

## Update all positions: chr, breaks, genes
geneDF$chrLeft <- 1000 # + 1000 on the left
geneDF$chrRight[geneDF$broken %in% "no"] <- geneDF$chrRight[geneDF$broken %in% "no"] + 1000
geneDF$chrRight[geneDF$broken %in% "yes"] <- geneDF$chrRight[geneDF$broken %in% "yes"] - geneDF$breaksize[geneDF$broken %in% "yes"]

# update: break
geneDF$breakstart <- NA
geneDF$breakstart[geneDF$broken %in% "yes"] <- geneDF$stopPart1[geneDF$broken %in% "yes"]  
geneDF$breakend[geneDF$broken %in% "yes"] <- geneDF$startPart2[geneDF$broken %in% "yes"] - 
  geneDF$breaksize[geneDF$broken %in% "yes"]

# update position: genes
geneDF$geneLeft <- geneDF$geneLeft + 1000
geneDF$geneRight<- geneDF$geneRight + 1000

geneDF$geneLeft[geneDF$whichpart %in% "part2"] <- geneDF$geneLeft[geneDF$whichpart %in% "part2"] - 
  geneDF$breaksize[geneDF$whichpart %in% "part2"]  
geneDF$geneRight[geneDF$whichpart %in% "part2"] <- geneDF$geneRight[geneDF$whichpart %in% "part2"] - 
  geneDF$breaksize[geneDF$whichpart %in% "part2"]

## Simplify names for legend
geneDF$IRGgroup <- geneDF$Name_Bekpen_2005
geneDF$IRGgroup[geneDF$Name_Bekpen_2005 %in% c("Irgb9b8", "Irgb1b2", "Irgb5" , "Irgb5*b3")] <- "Irgb tandems"
geneDF$IRGgroup[geneDF$Name_Bekpen_2005 %in% c("Irga2", "Irga4", "Irga3", "Irga6")] <- "Irga"
geneDF$IRGgroup[geneDF$Name_Bekpen_2005 %in% c("Irgb6*", "Irgb6") ] <- "Irgb6"
geneDF$IRGgroup[geneDF$Name_Bekpen_2005 %in% c("Irgm1", "Irgm2", "Irgm3") ] <- "Irgm"

# order species according to phylogeny:
geneDF$sp <- factor(geneDF$sp, levels = c("Itri", "Anil", "Mcou", "Rrat", "Rnor",
                                          "Mpah", "Mmus", "Mcar", "Aamp", "Moch", 
                                          "Otor", "Pman", "Pleu", "Cgri"))

# colors for plot:
mycolors <- c("yellow","red",6,"violet","darkblue","lightblue","green1","brown")

# to add a scale instead of x axis
scaleDF = data.frame(sp=geneDF$sp, x=300000, xend=310000, chr="NW_020822456.1", chr2="NW_020822466.1")
scaleDF <- scaleDF[scaleDF$sp %in% "Cgri",] # to keep the levels for factor
scaleDF$label <- NA
scaleDF$label[1] = "10kb"

# for gene label
geneDF$shortName <- gsub("Irg", "", geneDF$Name_Bekpen_2005)

# gene direction
geneDF$direction <- "forward"
geneDF$direction[geneDF$start > geneDF$stop] <- "reverse"

## and plot
syntenyPlot <- ggplot(geneDF) +
  # chr
  geom_segment(aes(x=chrLeft, xend = chrRight, y=chr, yend=chr), col = "grey", size = 1) +
  # break in long chr
  geom_segment(aes(x=breakstart+1000, xend=breakend, y=chr, yend=chr), col = "black", size = 7) +
  # GENES:
  # arrows for directions
  geom_segment(data=geneDF[geneDF$direction %in% "forward",],
               aes(x = geneLeft, y = chr, xend = geneRight, yend = chr, col = IRGgroup),
               arrow = arrow(length = unit(0.4, "cm"))) +
  geom_segment(data=geneDF[geneDF$direction %in% "reverse",],
               aes(x = geneRight, y = chr, xend = geneLeft, yend = chr, col = IRGgroup),
               arrow = arrow(length = unit(0.4, "cm"))) +
  # blocks for position
  geom_segment(aes(x=geneLeft, xend=geneRight, y=chr, yend=chr, col = IRGgroup), size = 3) +
  # genes labels
  geom_text(aes(x=geneLeft+(geneRight-geneLeft)/2, y=chr, label=shortName), size = 2) + 

  theme_bw() +
  # add scale
  geom_errorbarh(data=scaleDF,
                 aes(xmin=x,y=chr,xmax=xend), size = 1) +
  geom_text(data=scaleDF, aes(y=chr2, label=label), x=305000) + 
  facet_wrap(sp~., scales = "free_y", nrow = 14) +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(expand=c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text.x  = element_blank(), axis.ticks.x  = element_blank(), axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        panel.border=element_blank()) 

syntenyPlot

### SAVING POINT:
pdf(file = "../figures/SyntenyPlot.pdf", width = 13, height = 15)
syntenyPlot
dev.off()
