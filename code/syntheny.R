install.packages("genoPlotR") 

library(genoPlotR)
# 
# cli/DESCRIPTION', probable reason 'No such file or directory'cli/DESCRIPTION', probable reason 'No such file or directory'
# 2: In readLines(file, skipNul = TRUE) :
#   cannot open compressed file '/usr/lib/R/site-library/digest/DESCRIPTION', probable reason 'No such file or directory'
# 3: In readLines(file, skipNul = TRUE) :
#   cannot open compressed file '/usr/lib/R/site-library/rmarkdown/DESCRIPTION', probable reason 'No such file or directory'
# 4: In readLines(file, skipNul = TRUE) :
#   cannot open compressed file '/usr/lib/R/site-library/tidyselect/


data(three_genes)
comparisons[[1]]$col <- apply_color_scheme(c(0.6, 0.4, 0.5), "grey")
names <- c("Huey", "Dewey", "Louie")
names(dna_segs) <- names
tree <- newick2phylog("(((Huey:4.2,Dewey:3.9):3.1,Louie:7.3):1);")
mid_pos <- middle(dna_segs[[1]])
xlims <- list(c(Inf, -Inf), c(-Inf, Inf), c(1850, 2800))
annot <- annotation(x1=c(mid_pos[1], dna_segs[[1]]$end[2]),
                       x2=c(NA, dna_segs[[1]]$end[3]),
                       text=c(dna_segs[[1]]$name[1], "region1"),
                       rot=c(30, 0), col=c("blue", "black"))
#Now plotting these three segments:
 plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
                   annotations=annot, annotation_height=1.3,
                   tree=tree, tree_width=2,
                   xlims=xlims,
                   main="Comparison of Huey, Dewey and Louie")
 

my_dna_segs <- read.csv("segs.csv")
my_dna_segs <- split(my_dna_segs, f = my_dna_segs$species)


 