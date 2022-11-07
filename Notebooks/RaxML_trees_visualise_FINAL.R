# Script by Edeline Gagnon, 31.10.2022
# Visualisation of RaxML trees for Rees et al. (in prep) Paubrasilia population genomic paper.

library(ggtree)
library(ggtreeExtra)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ape)
library(treeio)

#Set working directory
setwd("C://Users/edeli/OneDrive/Projets_Recherche/2018_Paubrasilia/Manuscripts/Scripts/03_RaxML/analysis-raxml_data2_1000bp/")

#Load trees
tree<-read.raxml("RAxML_bipartitionsBranchLabels.data4_1000bp")
tree@phylo$tip.label


#search and replace the tip names
TipNames <- read.table("master_list_names.txt",header=T,sep="\t")
str(TipNames)
recoderFunc <- function(data, oldvalue, newvalue) {
  # convert any factors to characters
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  # create the return vector
  newvec <- data
  # put recoded values into the correct position in the return vector
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec
}
NewNames <- recoderFunc(tree@phylo$tip.label, TipNames$OldName, TipNames$New.Name)
NewNames
tree@phylo$tip.label <- NewNames

tree2<-drop.tip(tree, c("Cenostigma_pluviosum_77_26553","Cenostigma_pluviosum_84_12678","Cenostigma_laxiflorum_2012","Cenostigma_microphyllum_86_1596","Cenostigma_gaumeri_92_1762",
         "Cenostigma_marginatum_26561","Cenostigma_pluviosum_76_3105"))
tree2
            
# visualise tree without the outgroup

b <- ggtree(tree)+#,branch.length="none") + 
#  geom_tippoint()+
  # Add point when support for node is greater than 0.95 posterior probability)
    geom_point2(aes(subset=bootstrap>95),size=2, shape=23, fill="chartreuse")+
  # Adds tip labels
    geom_tiplab(color="black",size=3)+
  geom_nodelab(aes(label=round(bootstrap, 2), subset=bootstrap<95), vjust=-.3, hjust=1.5,size=3)+
  #Adds an axis
  theme_tree2() +
#  xlim(NA,25)
    xlim(NA,0.125)
  
b

b2 <- ggtree(tree2,branch.length="none") + 
  #  geom_tippoint()+
  # Add point when support for node is greater than 0.95 posterior probability)
  geom_point2(aes(subset=bootstrap>95),size=2, shape=23, fill="chartreuse")+
  # Adds tip labels
  geom_tiplab(color="black",size=3)+
  geom_nodelab(aes(label=round(bootstrap, 2), subset=bootstrap<95), vjust=-.3, hjust=1.5,size=3)+
  #Adds an axis
  theme_tree2() +
  xlim(NA,25)

b2

