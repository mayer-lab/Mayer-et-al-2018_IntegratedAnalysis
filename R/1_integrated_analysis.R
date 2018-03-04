# this script performs the following steps:
# 1) alignment of embryonic (E13, E18) and postnatal (P10) precursor cells with adult cortical interneurons 
# 2) run t-distributed stochastic neighbor embedding (t-SNE) dimensionality reduction on the integrated datasets. 
# 3) visualization of aligned datasets 
# load
install_github("mayer-lab/SeuratForMayer2018") #version of Seurat used in Mayer et al 2018
library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
library(reshape2, ggplot2)
library(gplots)
library(FNN)

# Pre-processed Seurat objects can be downloaded here: https://www.dropbox.com/s/qe2carqnf9eu4sd/Filtered_Mayer-et-al.Rda.zip?dl=0
load("~/Datashare/Filtered_Mayer-et-al.Rda") 
# color Palette URL: http://paletton.com/#uid=7531q0kw0w0jyC+oRxVy4oIDfjr
c.AT1 <- c("#0F4DA8","#FFA100","#BE008A","#B4F200","gray10","gray20","gray30","lightgray","gray95") #AT1 "Pvalb","Sst","Vip","Id2","Igfbp6","Th","Nos1","P10"
c.AT2 <- c("#547FBE","#0B3B82","#FFA100","#FFB639","#9B6200","#BE008A","#C22B99","#740054","#B4F200","#C3F336","#6D9300","gray10","gray30","gray50","lightgrey") #AT2 "Pvalb-1","Pvalb-2","Sst-1","Sst-2","Sst-3","Vip-1","Vip-2","Vip-3","Id2-1","Id2-2","Id2-3", "Igfbp6" "Th","Nos1","E13")
c.E1 <- c("snow2","lightsteelblue3","#0F4DA8","#BE008A","#FFA100","#B4F200","gray10","gray20","gray30")#MGE_E13" "CGE_E13" "Pvalb"   "Sst"     "Vip"     "Id2"     "TH"    "Nos1"    
c.E3 <- c("snow2","lightsteelblue3","lightgrey")
c.E1 <- c("snow2","lightsteelblue3","#0F4DA8","#BE008A","#FFA100","#B4F200","gray10","gray20","gray30")
c.br = c("#e41a1c","#4daf4a","#377eb8","lightgray")

# run a canonical correlation analysis to identify common sources of variation between two datasets
ccaE13=RunCCA(E13,P56)
ccaE18=RunCCA(E18,P56)
ccaP10=RunCCA(P10,P56)

# align the CCA subspaces. This returns a new dimensional reduction called cca.aligned. 
ccaE13=AlignSubspace(ccaE13,reduction.type = "cca",grouping.var = "DevStage",dims.align = 1:10) #cutoff pass: 1,2,4
ccaE18=AlignSubspace(ccaE18,reduction.type = "cca",grouping.var = "DevStage",dims.align = 1:10) #cutoff pass: 1:5
ccaP10=AlignSubspace(ccaP10,reduction.type = "cca",grouping.var = "DevStage",dims.align = 1:15) #cutoff pass: 1:9 

# run t-distributed stochastic neighbour embedding (t-SNE) dimensionality reduction on the integrated dataset
ccaE13 <- RunTSNE(ccaE13, reduction.use = "cca.aligned", dims.use = c(1,2,4)) #we selected CCVs for which at least 30 genes exhibited a minimum bicor of 0.15 in both datasets; see output of AlignSubspace
ccaE18 <- RunTSNE(ccaE18, reduction.use = "cca.aligned", dims.use = c(1:5))
ccaP10 <- RunTSNE(ccaP10, reduction.use = "cca.aligned", dims.use = c(1:9)) 

# visualization
#AdultTypes1 ccaP10
ccaP10@meta.data[, "AdultTypes1"] <- factor(ccaP10@meta.data[, "AdultTypes1"], levels = c("Pvalb","Sst","Vip","Id2","Igfbp6","Th","Nos1","P10")) 
TSNEPlot(ccaP10,do.return = F, pt.size = 0.001,label.size=4, do.label = T, group.by="AdultTypes1",colors.use=c.AT1,no.legend = T )

# AdultTypes1 ccaE18
ccaE18@meta.data[, "AdultTypes1"] <- factor(ccaE18@meta.data[, "AdultTypes1"], levels = c("Pvalb","Sst","Vip","Id2","Igfbp6","Th","Nos1","E18")) 
TSNEPlot(ccaE18,do.return = F, pt.size = 0.001,label.size=4,do.label = T, group.by="AdultTypes1",colors.use=c.AT1,no.legend = T,do.identify = F)

# AdultTypes1 ccaE13
ccaE13@meta.data[, "AdultTypes1"] <- factor(ccaE13@meta.data[, "AdultTypes1"], levels = c("Pvalb","Sst","Vip","Id2","Igfbp6","Th","Nos1","E13")) 
TSNEPlot(ccaE13,do.return = F, pt.size = 0.001,label.size=4,do.label = T, group.by="AdultTypes1",no.legend = T,colors.use=c.AT1,do.identify = FALSE)

# AdultTypes2 ccaP10
ccaP10@meta.data[, "AdultTypes2b"] <- factor(ccaP10@meta.data[, "AdultTypes2b"], levels = c("P56 Pvalb-1","P56 Pvalb-2","P56 Sst-1","P56 Sst-2","P56 Sst-3","P56 Vip-1","P56 Vip-2","P56 Vip-3","P56 Id2-1","P56 Id2-2","P56 Id2-3","P56 Igfbp6", "P56 Th","P56 Nos1","P10")) 
TSNEPlot(ccaP10,do.return = F, pt.size = 0.001,label.size=4, do.label = T, group.by="AdultTypes2b",colors.use=c.AT2,no.legend = F)

# AdultTypes2 ccaE18
ccaE18@meta.data[, "AdultTypes2b"] <- factor(ccaE18@meta.data[, "AdultTypes2b"], levels = c("P56 Pvalb-1","P56 Pvalb-2","P56 Sst-1","P56 Sst-2","P56 Sst-3","P56 Vip-1","P56 Vip-2","P56 Vip-3","P56 Id2-1","P56 Id2-2","P56 Id2-3","P56 Igfbp6", "P56 Th","P56 Nos1","E18")) 
TSNEPlot(ccaE18,do.return = F, pt.size = 0.001,label.size=4,do.label = T, group.by="AdultTypes2b",colors.use=c.AT2,no.legend = F,do.identify = F)

# AdultTypes2 ccaE13
ccaE13@meta.data[, "AdultTypes2b"] <- factor(ccaE13@meta.data[, "AdultTypes2b"], levels = c("P56 Pvalb-1","P56 Pvalb-2","P56 Sst-1","P56 Sst-2","P56 Sst-3","P56 Vip-1","P56 Vip-2","P56 Vip-3","P56 Id2-1","P56 Id2-2","P56 Id2-3","P56 Igfbp6", "P56 Th","P56 Nos1","E13")) 
TSNEPlot(ccaE13,do.return = F, pt.size = 0.001,label.size=4,do.label = T, group.by="AdultTypes2b",no.legend = F,colors.use=c.AT2,do.identify = FALSE)

# eminence 
temp1 = rownames(table(ccaE13@meta.data$Eminence))
ccaE13@meta.data[, "Eminence"] <- factor(ccaE13@meta.data[, "Eminence"], levels = temp1[c(3,1,4,2)]) 
TSNEPlot(ccaE13,do.return = F, pt.size = 0.001,label.size=4,do.label = T, group.by="Eminence",no.legend = F,colors.use=c.AT1,do.identify = FALSE)
#
temp1 = rownames(table(ccaE13@meta.data$Eminence1))
ccaE13@meta.data[, "Eminence1"] <- factor(ccaE13@meta.data[, "Eminence1"], levels = c("E13 MGE","E13 CGE","Pvalb","Sst","Vip","Id2","Igfbp6","Th","Nos1"))
TSNEPlot(ccaE13,do.return = F, pt.size = 0.001,label.size=4,do.label = F, group.by="Eminence1",no.legend = F,colors.use=c.E1,do.identify = FALSE)


