library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(cowplot)

data <- read.csv("data4syntheny.tsv", sep = "\t")

data$middleGene <- data$end - (data$end - data$start)/2

# to obtain chromosome ratio (useful at the end)
df1 <- data %>% group_by(Species,chr) %>% summarise(Chrlength=max(end)-min(start))
## keep tiny bits visible
df1$Chrlength[df1$Chrlength < 50000] <- 100000
df2 <- df1 %>% group_by(Species) %>% summarise(AlChrlength=sum(Chrlength))
df3 <- merge(df1, df2)
df3$ChrRatio <- df3$Chrlength / df3$AlChrlength
data <- merge(data, df3)

# order chromosomes 
data <- data[order(data$chr),]
d1 <- data[grep("Chr", data$chr),]
d1 <- d1[order(as.numeric(gsub("_part2", "", gsub("Chr","",d1$chr)))),]
d2 <- data[!grepl("Chr", data$chr),]
data <- rbind(d1, d2)

speciesList <- unique(data$Species)
chrList <- unique(data$chr)
data$combi <- paste0(data$Species, data$chr)

# divisor of length of the ... extra part of chromosomes
div = 5

all_plot_list = list()
for (s in speciesList){
  
  plot_list = list()
  for (c in chrList){
    # check if the species-chromosome combination exists
    if(paste0(s, c) %in% data$combi){
      # if yes, plot
      subdata <- data[data$Species %in% s & data$chr %in% c,]
      
      # choose colors
      subdata$color <- "black"
      subdata$color[subdata$IrgName %in% "IrgbTandems"] <- "red"
      subdata$color[subdata$IrgName %in% "Irgb6-6*"] <- "blue"
      subdata$color[subdata$IrgName %in% "Irgb10"] <- "green"
      subdata$color[subdata$IrgName %in% "Irga2-6"] <- "yellow"
      
      # prepare graduated axis along chr
      step=100000
      xl = seq(0,max(subdata$end), step)
      yl = rep(-0.55, length(seq(0,max(subdata$end), step)))
      labels=paste0(seq(0,max(subdata$end), step)/1000, "kb")
      plot.labels<-data.frame(xl,yl,labels)
      # remove all the useless ones before start
      plot.labels<-plot.labels[plot.labels$xl>min(subdata$start),]
      
      p <- ggplot(subdata)+
        # chr
        geom_rect(aes(xmin=min(start)-(max(end) - min(start))/div, xmax=max(end)+ (max(end) - min(start))/div, 
                      ymin=-1.5,ymax=1), fill = "white")+
        geom_text(aes(x = max(end) - (max(end) - min(start))/2, y = -0.9, label=chr), size = 4)+
        geom_rect(aes(xmin=min(start)-(max(end) - min(start))/div, xmax=max(end)+ (max(end) - min(start))/div, 
                      ymin=-0.5,ymax=+0.5), fill = "grey") +
        # # genomic position
        geom_segment(data = plot.labels, aes(x=xl, xend = xl, y = yl, yend = rep(-0.5, length(yl)))) +
        geom_text(data = plot.labels, aes(x = xl, y = yl-0.07, label=labels), size = 3, angle = 45)+
        ## hline
        geom_segment(aes(x=min(start),xend = max(end), y = -0.75, yend = -0.75)) +
        geom_segment(aes(x=min(start), xend =min(start)-(max(end) - min(start))/div,
                         y = -0.75, yend = -0.75), linetype = 2) +
        geom_segment(aes(x=max(end), xend =max(end)+(max(end) - min(start))/div,
                         y = -0.75, yend = -0.75), linetype = 2) +
        ## ticks
        geom_segment(aes(x=min(start)-(max(end) - min(start))/div, xend = min(start)- (max(end) - min(start))/div,
                         y = -0.75,  yend = -0.65)) +
        geom_segment(aes(x=max(end)+(max(end) - min(start))/div, xend = max(end)+ (max(end) - min(start))/div,
                         y = -0.75,  yend = -0.65)) +
        # gene
        geom_rect(aes(xmin=start, xmax=end, ymin=-0.5,ymax=+0.5, fill=color), col="black", alpha = 0.7) +
        scale_fill_manual(values = unique(subdata$color))+
        geom_label_repel(aes(x = middleGene, y = 0.5, label=IrgName), size = 4, nudge_y = 0.1)+
        
        # theme
        theme_void() +
        theme(legend.position = "none")
      plot_list[[c]] = p
    }
  }
  all_plot_list[[s]] = plot_list
}

# for chr ratio
drat <- unique(data.frame(data$Species, data$chr,data$ChrRatio))

nspe <- length(all_plot_list)
all_plots <- list()
for(s in 1:nspe){
  n <- length(all_plot_list[[s]])
  all_plots[[s]] <- plot_grid(plotlist = all_plot_list[[s]], ncol=n, 
                              rel_widths = drat[drat$data.Species %in% speciesList[s],"data.ChrRatio"], labels = speciesList[s])
}
do.call("plot_grid", c(all_plots, nrow=s))