library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(cowplot)

data <- read.csv("../data/data4syntheny.csv", sep = "\t")

# manual addition for faceting (for 2 species with only facet 1)
data <- rbind(data,
              data.frame(IrgName = "TO REMOVE", Species="M.glareolus", chr="LR761921.1", start_NCBI=24392229, start=24392229,
                         end=24392259, forward_reverse = NA, gene_name=NA, protein_name=NA, exon_nbr = NA, found_by=NA,
                         rankSpecies=7),
              data.frame(IrgName = "TO REMOVE", Species="M.auratus", chr="NW_004801614.1", start_NCBI=71000000, start=71000000,
                         end=71000500, forward_reverse = "F", gene_name=NA, protein_name=NA, exon_nbr = NA, found_by=NA,
                         rankSpecies=9))

# order chromosomes
data$chr <- as.factor(data$chr)
data <- data[order(data$chr),]
d1 <- data[grep("Chr", data$chr),]
d1 <- d1[order(as.numeric(gsub("b", "", gsub("Chr","",d1$chr)))),]
d2 <- data[!grepl("Chr", data$chr),]
data <- rbind(d1, d2)
# add counter per species
data <- data %>% group_by(Species) %>% mutate(rankChr = dense_rank(chr))
#add full represented chr length
d1 <- data %>% group_by(Species, chr) %>%
  summarise(chrStart=min(start), chrEnd=max(end), chrLength=max(end)-min(start))
data <- merge(data, d1)
# set x fact: is gene far from chr start? if yes, break
facetBreak <- 2500000
data$facetX <- "part1"
data$facetX[data$start - data$chrStart > facetBreak] <- "part2"

d1 <- data[data$facetX %in% "part1",] %>% group_by(Species, chr) %>%
  summarise(part1chrStart=min(start), part1chrEnd=max(end), part1chrLength=max(end)-min(start))
d2 <- data[data$facetX %in% "part2",] %>% group_by(Species, chr) %>%
  summarise(part2chrStart=min(start), part2chrEnd=max(end), part2chrLength=max(end)-min(start))
d3 <- merge(d1, d2, all.x = T)
data <- merge(data, d3, by=c("Species", "chr"), all.x = T)
rm(d1, d2, d3)

# set up start and end for chromosomes
data$facetX_chrStart <- data$part1chrStart
data$facetX_chrStart[data$facetX %in% "part2"] <- data$part2chrStart[data$facetX %in% "part2"] 
data$facetX_chrEnd <- data$part1chrEnd
data$facetX_chrEnd[data$facetX %in% "part2"] <- data$part2chrEnd[data$facetX %in% "part2"] 
data$facetX_chrLength <- data$facetX_chrEnd - data$facetX_chrStart

# set up end and start of arrow
data$arrowStart <- data$start - data$facetX_chrStart
data$arrowEnd <- data$start - data$facetX_chrStart + 10000
data$arrowStart[data$forward_reverse %in% "R"] <- data$end[data$forward_reverse %in% "R"] - data$facetX_chrStart[data$forward_reverse %in% "R"]
data$arrowEnd[data$forward_reverse %in% "R"] <- data$end[data$forward_reverse %in% "R"] - data$facetX_chrStart[data$forward_reverse %in% "R"] - 10000

# choose colors
group.colors <- c(IrgbTandem = "blue", "Irgb6-OG1" = "red", "Irgb6-OG2" = "orange", Irgb10 ="green", "Irga2-6" = "yellow",
                  Irga4=1, Irga3=1, Irgm1=1, Irgc=1, Irgm2=1, Irgd=1, Irgq=1, Irgm3=1, "TO REMOVE"="white")

makePlotPerSpecies <- function(df){
  ann_text_chr <- data.frame(x = - 80000, y = unique(df$rankChr)+0.2,
                             lab = unique(df$chr), fontface="bold.italic",
                             facetX = factor("part1",levels = c("part1","part2")))
  
  ggplot(df)+
    # frame
    geom_rect(aes(xmin=-100000, xmax= max(data$facetX_chrLength), 
                  ymin=0,ymax=max(rankChr)+0.5), fill = "white", col = "white")+
    # chr
    geom_rect(aes(xmin=0, xmax= facetX_chrLength, ymin=rankChr-0.2,ymax=rankChr+0.2), fill = "lightgrey", col = 1, alpha=.4)+
    scale_y_reverse(breaks=1:length(unique(df$chr)), labels = unique(df$chr))+
    # annotations: chromosome
    geom_label(data = ann_text_chr, aes(x=x+10000, y=y-0.25, label =lab), size = 5, fill = "black", col = "white") +
    # genes
    geom_rect(aes(xmin=start - facetX_chrStart, xmax=end - facetX_chrStart, ymin=rankChr-0.2,ymax=rankChr+0.2, fill=IrgName),
              col="black", alpha = 0.7) +
    scale_fill_manual(values = group.colors)+
    # arrow in genes
    geom_segment(aes(x=arrowStart, xend=arrowEnd, y=rankChr, yend=rankChr),
                 arrow= arrow(length=unit(0.30,"cm"),type = "closed"), color ="black") +
    # genes labels 
    geom_label_repel(aes(x = (end + start)/2 - facetX_chrStart, y = rankChr, label=IrgName), size = 4, nudge_y = 0.1)+
    # genomic position
    geom_text_repel(aes(x = start - facetX_chrStart, y = rankChr+0.2, label=start), 
                    force = 0.5, nudge_y = -0.2, direction = "x", vjust = 0, segment.size = 0.2)+
    # facet x
    facet_grid(.~facetX, ) +
    # theme
    theme_nothing() +
    theme(legend.position = "none", panel.border=element_blank(),
          plot.background = element_rect(color = "black", size = 10))
}

muridae <- c("M.musculus", "M.pahari","R.norvegicus","R.rattus")
cricetidae <- c("M.ochrogaster","A.amphibius", "M.glareolus", "P.leucopus", "M.auratus")

plot.list_mus <- list()
ratio_mus <- vector()
for (s in muridae){
  p <- makePlotPerSpecies(data[data$Species %in% s,])
  ratio_mus[s] <- length(unique(data[data$Species %in% s,"chr"])) +1
  plot.list_mus[[s]] <- p
}
plot_mus <- plot_grid(plotlist = plot.list_mus, nrow=length(muridae), labels = muridae, label_fontface = "bold.italic", 
                      label_x = 0, label_y = 0.9, label_size = 20, rel_heights = ratio_mus)
pdf(file = "../figures/SyntenyPlot_muridae.pdf", width = 25, height =20)
plot_mus
dev.off()

plot.list_cri <- list()
ratio_cri <- vector()
for (s in cricetidae){
  p <- makePlotPerSpecies(data[data$Species %in% s,])
  ratio_cri[s] <- length(unique(data[data$Species %in% s,"chr"])) +1
  plot.list_cri[[s]] <- p
}
plot_cri <- plot_grid(plotlist = plot.list_cri, nrow=length(cricetidae), labels = cricetidae, label_fontface = "bold.italic",
                      label_x = 0, label_y = 0.9, label_size = 20, rel_heights = ratio_cri)
pdf(file = "../figures/SyntenyPlot_cricetidae.pdf", width = 25, height =25)
plot_cri
dev.off()
