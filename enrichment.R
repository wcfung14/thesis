# To perform KOG/KEGG/GO enrichment by 'compareCluster' of 'clusterProfiler' package
# reference to https://github.com/xieyichun50/Myriapod-genomes/blob/v1.0.0/script/4function_anno2tree/04gain_loss_specific_enrichment.R and 05dotplot.R
# input:
#     Tco_S7F_head_10-6_6h_20E_DMSO_all_rep3_swapped.csv # results from Degust (http://degust.erc.monash.edu/)
#     Tco.GO.1v1.txt
#     Tco.KEGG.1v1.txt
#     Tco.ko.1v1.txt
#     Tco.KOG.1v1.txt

library(clusterProfiler)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../")

degust_res <- read.csv("./degust/Tco_S7F_head_10-6_6h_20E_DMSO_all_rep3_swapped.csv") # Results from Degust (http://degust.erc.monash.edu/)
degust_res$direction <- ifelse(degust_res$exp < 0, "DOWN", "UP")
head(degust_res)

# Compare Cluster function
compare_cluster <- function(term2gene, term2name, maxgssize) {
  test <- compareCluster(gene_id ~ direction, data = degust_res, 
                         fun = "enricher", TERM2GENE = term2gene, TERM2NAME = term2name,
                         pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                         minGSSize = 1, 
                         maxGSSize = maxgssize)
  
  plotin <- as.data.frame(test)
  return(plotin)
}

# KOG
KOG <- read.delim("./Tco.KOG.1v1.txt")
kog2name <- read.delim("./kog2name.txt")

KOG.all <- compare_cluster(term2gene = KOG, term2name = kog2name, maxgssize = 200000)

KOG.all<-KOG.all[which(is.na(KOG.all$Description)==FALSE),]
names(KOG.all)[3]<-"kogClass"
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)
KOG.all<-merge(KOG.all, kog2name, by = c("kogClass"), all.x = TRUE)
KOG.all<-unique(KOG.all)

# KEGG
KEGG <- read.delim("./Tco.KEGG.1v1.txt")
kegg2name <- read.delim("./kegg2name.txt")

KEGG.all <- compare_cluster(term2gene = KEGG, term2name = kegg2name, maxgssize = 200000)

# GO
GO <- read.delim("./Tco.GO.1v1.txt")
go2name <- read.delim("./go2name.txt")

GO.all <- compare_cluster(term2gene = GO, term2name = go2name, maxgssize = 10000000)

# need to be done by GO levels

# Dotplot function
plotdata_func <- function(plotin) {
  plotinsep <- separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep <- separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotinsep$ratio1 = plotinsep$Genenumerator/plotinsep$BGnumerator
  plotinsep$ratio2 = plotinsep$Genenumerator/plotinsep$Genedenominator
  #names(plotinsep)[names(plotinsep)=="direction"]<-"label"
  
  plotdata <- subset(plotinsep)
  plotdata$ratio1 <- ifelse(plotdata$direction == "DOWN", -plotdata$ratio1, plotdata$ratio1)
  plotdata$ratio2 <- ifelse(plotdata$direction == "DOWN", -plotdata$ratio2, plotdata$ratio2)
  
  labelorderGE <- unique(subset(plotdata, select = c("direction","Genedenominator")))
  #labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
  
  return(plotdata)
}

plot_func <- function(plotdata, title_text, text_size=11) {
  topgroups <- plotdata %>% group_by(Cluster)
  
  plot <- ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = ratio2, y = Description, size = Genenumerator, colour = p.adjust)) +
    labs(title = title_text, x = "Gene ratio", size = "Gene counts", colour = "p.adjust") +
    geom_point(shape = 19)+scale_size_area() +
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2)) +
    scale_x_continuous(limits = c(-0.6, 0.6)) +
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "black"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
          axis.text.y = element_text(size = text_size),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12)) +
    #facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
    theme(strip.text = element_text(size = 10))
  return(plot)
}

KOG_plotdata <- plotdata_func(plotin = KOG.all)
KEGG_plotdata <- plotdata_func(plotin = KEGG.all)

KOG_plot <- plot_func(KOG_plotdata, title_text = "KOG function enrichment")
KEGG_plot <- plot_func(KEGG_plotdata, title_text = "KEGG enrichment", text_size=6)

ggsave("KOG.tiff", plot = KOG_plot, width = 12, height = 12, units = "in", dpi = 300, compression="lzw")
ggsave("KOG.png", plot = KOG_plot, width = 12, height = 12, units = "in", dpi = 300)
ggsave("KEGG.tiff", plot = KEGG_plot, width = 12, height = 12, units = "in", dpi = 300, compression="lzw")
ggsave("KEGG.png", plot = KEGG_plot, width = 12, height = 12, units = "in", dpi = 300)

write.csv(x = KEGG_plotdata, "KEGG_plotdata.csv",row.names = FALSE)
