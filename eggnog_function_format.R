# To format annotation for each protein in the eggNOG annotation
# reference to https://github.com/xieyichun50/Myriapod-genomes/blob/v1.0.0/script/4function_anno2tree/Orthofinder_eggnog_orthogroup_function.R
# input: 
#     Trigoniulus_corallinus_hic.proteins.fa.v2.emapper.annotations.xls
#     go2name.txt
#     kegg2name.txt
#     ko2name.txt
#     kog2name.txt
# output:
#     Tco.GO.1v1.txt
#     Tco.KEGG.1v1.txt
#     Tco.ko.1v1.txt
#     Tco.KOG.1v1.txt

library(dplyr)
library(tidyr)
library(stringr)
#library(ggplot2)
#library(clusterProfiler)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../")

#Whole genome Annotation file
prefix = "Trigoniulus_corallinus_hic.proteins.fa.v2"
eggnog <- read.delim(paste0(prefix,".emapper.annotations.xls"), header = TRUE, skip = 3)
eggnog <- separate(eggnog, X.query_name, c("Genes", "TS"), sep = "-T", remove = FALSE)

##Gene name
gene2name <- eggnog[,c("Genes","Preferred_name", "X.4")]
colnames(gene2name)[colnames(gene2name) == "X.4"] <- "Description"
gene2name <- gene2name[gene2name$Description != "",]
gene2name <- unique(gene2name)

write.table(gene2name, file = "gene2name.txt", sep = '\t', row.names = FALSE, quote = FALSE)

# function for formatting GenesKEGGpair/Geneskopair/GenesKOGpair/GenesGOpair
df_format_func <- function(df, target_column) {
  df_new <- data.frame()
  for (row in 1:nrow(df)) {
    rcnames <- list(df[row, "Genes"], c(strsplit(df[row, target_column], ',')[[1]]))
    rcnames[[1]] <- rep(rcnames[[1]], length(rcnames[[2]]))
    
    df_new <- rbind(df_new, rcnames)
  }
  colnames(df_new) <- c("Genes", target_column)
  df_new <- unique(df_new)
  return(df_new)
}

##KEGG from eggnog
{
  #read in kegg2name
  kegg2name <- read.delim("./kegg2name.txt", sep = "\t", colClasses = "character")
  
  pathways <- eggnog[, c("Genes","KEGG_Pathway")]
  pathways$KEGG_Pathway <- gsub(",map.*$", "", pathways$KEGG_Pathway)
  colnames(pathways) = c("Genes", "KEGG")
  pathways <- subset(pathways, KEGG != "-" & KEGG != "" & is.na(KEGG)==FALSE, select = c("Genes", "KEGG"))
  pathways <- unique(pathways)
 
  #Format GenesKEGGpair
  GenesKEGGpair.1v1 <- df_format_func(pathways, target_column = "KEGG")
  GenesKEGGpair.1v1 <- subset(GenesKEGGpair.1v1, select = c("KEGG", "Genes"))

  write.table(GenesKEGGpair.1v1, file = "./Tco.KEGG.1v1.txt", sep = '\t', row.names = FALSE, quote = FALSE)
  rm(pathways)
}

##KO from eggnog
{
  #read in kegg2name
  ko2name <- read.delim("./ko2name.txt", sep = "\t", colClasses = "character")
  
  pathways <- eggnog[,c("Genes","KEGG_ko")]
  colnames(pathways) = c("Genes", "ko")
  pathways <- subset(pathways, ko != "-" & ko != "" & is.na(ko)==FALSE, select = c("Genes", "ko"))
  pathways <- unique(pathways)
  
  #Format Geneskopair
  Geneskopair.1v1 <- df_format_func(pathways, target_column = "ko")
  Geneskopair.1v1 <- subset(Geneskopair.1v1, select = c("ko", "Genes"))
  
  write.table(Geneskopair.1v1, file = "./Tco.ko.1v1.txt", sep = '\t', row.names = FALSE, quote = FALSE)
  rm(pathways)
}

##KOG from eggnog
{
  #read in kog2name
  kog2name <- read.delim("./kog2name.txt", sep = "\t", colClasses = "character")
  
  pathways <- eggnog[,c("Genes","X.3")] 
  colnames(pathways) = c("Genes", "KOG") 
  pathways <- subset(pathways, KOG != "" & KOG != "-" & is.na(KOG)==FALSE)
  pathways <- unique(pathways)
  
  #Format GenesKOGpair
  GenesKOGpair.1v1 <- df_format_func(pathways, target_column = "KOG")
  GenesKOGpair.1v1 <- subset(GenesKOGpair.1v1, select = c("KOG", "Genes"))
  
  write.table(GenesKOGpair.1v1, file = "./Tco.KOG.1v1.txt", sep = '\t', row.names = FALSE, quote = FALSE)
  rm(pathways)
}

##GO from eggnog
{
  #read in go2name
  go2name <- read.delim("./go2name.txt", sep = "\t", colClasses = "character", header = FALSE)
  colnames(go2name) = c("goClass", "goName", "ONTOLOGY")
   
  pathways <- eggnog[,c("Genes", "GOs")]
  colnames(pathways) = c("Genes", "GO")
  pathways <- subset(pathways, GO != "" & GO != "-" & is.na(GO)==FALSE)
  pathways <- unique(pathways)
  
  #Format GenesGOGpair
  GenesGOpair.1v1 <- df_format_func(pathways, target_column = "GO")
  GenesGOpair.1v1 <- subset(GenesGOpair.1v1, select = c("GO", "Genes"))
  
  write.table(GenesGOpair.1v1, file = "./Tco.GO.1v1.txt", sep = '\t', row.names = FALSE, quote = FALSE)
  rm(pathways)
}

# microbenchmark:: microbenchmark(
#   expression_1 = {
#     df_format_func(pathways, target_column = "GO")
#   },
#   expression_2 = {
#     GenesGOpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
#     GenesGOpair.1v1<-as.data.frame(GenesGOpair.1v1)
#     names(GenesGOpair.1v1)[1]="Genes"
#     names(GenesGOpair.1v1)[2]="GO"
#     
#     for (i in 1:nrow(pathways)) {
#       subtable<-pathways[i,]
#       rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$GO[1], ',')[[1]]))
#       pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
#       pairtable<-as.data.frame(pairtable)
#       pairtable$Genes<-rownames(pairtable)
#       rownames(pairtable)<-1:nrow(pairtable)
#       pairtable<-as.data.frame(pairtable)
#       pairtable.new<-pairtable %>% gather(GO, pair, c(1:ncol(pairtable)-1))
#       pairtable.new<-pairtable.new[,c(1:2)]
#       GenesGOpair.1v1<-rbind(GenesGOpair.1v1, pairtable.new)
#     }
#     GenesGOpair.1v1<-subset(GenesGOpair.1v1, 
#                             is.na(GenesGOpair.1v1$Genes)==FALSE, 
#                             select = c("GO", "Genes"))
#     GenesGOpair.1v1<-unique(GenesGOpair.1v1)
#   },
#   times = 3
# )
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# expression_1 101.8941 103.8538 105.3055 105.8135 107.0112 108.2089     3
# expression_2 132.5566 132.6501 137.5675 132.7435 140.0729 147.4022     3
