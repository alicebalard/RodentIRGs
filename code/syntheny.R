  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(gridExtra)
  library(ggrepel)
  library(cowplot)
  
  data <- read.csv("../data/data4syntheny.csv", sep = "\t")
  
  # order chromosomes
  data$chr <- as.factor(data$chr)
  data <- data[order(data$chr),]
  d1 <- data[grep("Chr", data$chr),]
  d1 <- d1[order(as.numeric(gsub("b", "", gsub("Chr","",d1$chr)))),]
  d1$chr <- factor(d1$chr, levels = unique(d1$chr)) 
  d2 <- data[!grepl("Chr", data$chr),]
  data <- rbind(d1, d2)
  # add counter per species
  data <- data %>% group_by(Species) %>% mutate(rankChr = dense_rank(chr))
  
  #add full represented chr length
  d1 <- data %>% group_by(Species, chr) %>%
    summarise(chrStart=min(start), chrEnd=max(end), chrLength=max(end)-min(start))
  data <- merge(data, d1)
  # set x break: is gene far from chr start? if yes, break
  breakStart <- 2500000
  data$part <- "part1"
  data$part[data$start - data$chrStart > breakStart] <- "part2"
  
  d1 <- data[data$part %in% "part1",] %>% group_by(Species, chr) %>%
    summarise(part1chrStart=min(start), part1chrEnd=max(end), part1chrLength=max(end)-min(start))
  d2 <- data[data$part %in% "part2",] %>% group_by(Species, chr) %>%
    summarise(part2chrStart=min(start), part2chrEnd=max(end), part2chrLength=max(end)-min(start))
  data <- merge(data,merge(d1, d2, all = T))
  data$part2chrLength[is.na(data$part2chrLength)] <- 0
  
  data$isBreak <- F
  data$isBreak[data$part2chrLength > 1] <- T
  
  # coord for plot
  ## Part 1
  data$geneStartPlot <- data$start - data$chrStart
  data$geneEndPlot <- data$end - data$chrStart
  ## Part 2
  breakSizePlot <- 50000
  # reducton: real end - break plot end
  reduct <- data$part2chrStart[data$part %in% "part2"] - (
    data$part1chrEnd[data$part %in% "part2"] - data$chrStart[data$part %in% "part2"] + breakSizePlot)
  # data$start - data$chrStart - reduction
  data$geneStartPlot[data$part %in% "part2"] <- data$start[data$part %in% "part2"]  - reduct
  data$geneEndPlot[data$part %in% "part2"] <- data$end[data$part %in% "part2"]  - reduct
  
  # # set up end and start of arrow
  data$arrowStart <- data$geneStartPlot
  data$arrowEnd <- data$geneEndPlot
  data$arrowStart[data$forward_reverse %in% "R"] <- data$geneEndPlot[data$forward_reverse %in% "R"]
  data$arrowEnd[data$forward_reverse %in% "R"] <- data$geneStartPlot[data$forward_reverse %in% "R"] 
  
  # choose colors
  group.colors <- c(IrgbTandem = "blue", "Irgb6-OG1" = "red", "Irgb6-OG2" = "orange", "Irgb6-tbn" = "pink",
                    Irgb10 ="green", "Irga2-6" = "yellow",
                    Irga4=1, Irga3=1, Irgm1=1, Irgc=1, Irgm2=1, Irgd=1, Irgq=1, Irgm3=1)
  
  makePlotPerSpecies <- function(df){
    ggplot(df)+
      # space for species name on top as label in plot_grid
      geom_rect(aes(xmin=0, xmax = 100000, ymin=0, ymax=0), col = "white") +
      # chr part1 & part2
      geom_rect(aes(xmin=0, xmax= part1chrLength, 
                    ymin=rankChr-0.2,ymax=rankChr+0.2), fill = "lightgrey", col = "white", alpha=.4)+
      geom_rect(aes(xmin=part1chrLength+breakSizePlot, xmax= part1chrLength + breakSizePlot + part2chrLength, 
                      ymin=rankChr-0.2,ymax=rankChr+0.2), fill = "lightgrey", col = "white", alpha=.4)+
      scale_y_reverse(breaks=1:length(unique(df$chr)), labels = unique(df$chr))+
      # make a cut for long chromosomes
      geom_segment(aes(x=part1chrEnd-part1chrStart, xend= part1chrEnd-part1chrStart+breakSizePlot,
                    y=rankChr,yend=rankChr, col = isBreak), linetype = 2, size =2)+
      scale_color_manual(values = c("white", "black"))+
      # genes
      geom_rect(aes(xmin=geneStartPlot, xmax=geneEndPlot, ymin=rankChr-0.2,ymax=rankChr+0.2, fill=IrgName),
                col="black", alpha = 0.7) +
      scale_fill_manual(values = group.colors)+
      # arrow in genes
      geom_segment(aes(x=arrowStart, xend=arrowEnd, y=rankChr, yend=rankChr),
                   arrow= arrow(length=unit(0.40,"cm"),type = "closed"), color ="white") +
      # genes labels
      geom_label_repel(aes(x = (geneEndPlot + geneStartPlot)/2, y = rankChr, label=IrgName), size = 4, nudge_y = 0.1)+
    # genomic position
    geom_text_repel(aes(x = geneStartPlot, y = rankChr+0.2, label=start),
                    force = 0.5, nudge_y = -0.2, direction = "x", vjust = 0, segment.size = 0.2)+
    # theme
    theme_classic() +
      scale_x_continuous(expand = c(0,0)) +
      theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line =  element_blank(),
            legend.position = "none", panel.border=element_blank(), 
            axis.text.y = element_text(size = 10), axis.text.x = element_blank(),
            plot.background = element_rect(color = "black", size = 2))
  }
  
  muridae <- c("M.musculus", "M.pahari","R.norvegicus","R.rattus")
  cricetidae <- c("M.ochrogaster","A.amphibius", "M.glareolus", "P.leucopus", "M.auratus")
  
  makeSubPlot <- function(gp){
    plot.list <- list()
    ratio <- vector()
    for (s in gp){
      p <- makePlotPerSpecies(data[data$Species %in% s,])
      ratio[s] <- length(unique(data[data$Species %in% s,"chr"])) +1
      plot.list[[s]] <- p
    }
    plot_grid(plotlist = plot.list, nrow=length(gp), rel_heights = ratio, align = "hv",  
              labels = gp, label_fontface = "bold.italic", label_size = 20)
  }
  
  pdf(file = "../figures/SyntenyPlot_muridae.pdf", width = 20, height =20)
  makeSubPlot(muridae)
  dev.off()
  
  pdf(file = "../figures/SyntenyPlot_cricetidae.pdf", width = 20, height =25)
  makeSubPlot(cricetidae)
  dev.off()