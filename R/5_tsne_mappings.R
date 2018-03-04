# this script performs the following steps: 
# 1) visualize the assignment of precursor cells with an independent t-SNE analysis of each time point. 
# 2) visualize mapping of E18 precursor cells to E13 precursor states
# 3) P56 dataset (P56 raw data downloaded from the Allen Institute)

# tSNE
P10@meta.data$Map1=ccaP10@meta.data[P10@cell.names,"Map1"]
P10@meta.data$Map3=ccaP10@meta.data[P10@cell.names,"Map3"]
#
E18@meta.data$Map1=ccaE18@meta.data[E18@cell.names,"Map1"]
E18@meta.data$Map3=ccaE18@meta.data[E18@cell.names,"Map3"]
#
E13@meta.data$Map1=ccaE13@meta.data[E13@cell.names,"Map1"]
E13@meta.data$Map3=ccaE13@meta.data[E13@cell.names,"Map3"]

# visualization of assignment 
#P10
P10@meta.data[, "Map3"] <- factor(P10@meta.data[, "Map3"], levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned")) 
TSNEPlot(P10,pt.size = 0.1,label.size=4,group.by="Map3",colors.use=c.AT2)
#E18
E18@meta.data[, "Map3"] <- factor(E18@meta.data[, "Map3"], levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned")) 
TSNEPlot(E18,pt.size = 0.1,label.size=4,group.by="Map3",colors.use=c.AT2)
#E13
E13@meta.data$Map1 <- as.character(E13@meta.data$Map1)
E13@meta.data[, "Map1"] <- plyr::mapvalues(E13@meta.data$Map1, from = NA, to = "Unassigned")
E13@meta.data[, "Map1"] <- factor(E13@meta.data[, "Map1"], levels = c("Pvalb","Sst","Vip","Id2","Igfbp6","Nos1","Th","Unassigned")) 
TSNEPlot(E13,pt.size = 0.1,label.size=4, group.by="Map1",colors.use=c("#0F4DA8" ,"#FFA100","#BE008A","#B4F200","gray10" ,"gray20" ,"lightgray","gray95"))

# branch mappings from branch analysis on E18 
CH1 = readRDS("~/CH/Fig4a_10XE18_data_frame.Rds") #(Data see GitHub of Christop Hafemeister)
which.cells = rownames(CH1[grep("CX",CH1$sample),])
E18@meta.data[which.cells,"Branch"] <- CH1[which.cells, "Branch"]
E18@meta.data[which.cells,"BranchEminence"] <- CH1[which.cells, "label"]
CH1[, "label2"] <- as.character(CH1[, "label2"])
CH1[, "label2"]  <- plyr::mapvalues(CH1[, "label2"], from = c(NA), to = c("Unassigned"))
E18@meta.data[which.cells,"BranchEminence2"] <- CH1[which.cells, "label2"]
E18@meta.data[which.cells,"Branch2"] <- plyr::mapvalues(E18@meta.data[which.cells,"BranchEminence2"], from = c("Lhx6+ 1","Lhx6- 1","Lhx6+ 2","Lhx6- 2","Lhx6+ 3","Lhx6- 3"), to = c("1","1","2","2","3","3"))
E18@meta.data[which.cells,"Eminence3"] <- plyr::mapvalues(E18@meta.data[which.cells,"BranchEminence2"], from = c("Lhx6+ 1","Lhx6- 1","Lhx6+ 2","Lhx6- 2","Lhx6+ 3","Lhx6- 3"), to = c("MGE","CGE","MGE","CGE","MGE","CGE"))
E18@meta.data[, "Eminence4"] <- factor(E18@meta.data[, "Eminence3"], levels = c("MGE","CGE","Unassigned")) 

#branch
TSNEPlot(E18,do.return = F, pt.size = 0.001,label.size=4,do.label = F, group.by="Branch2",colors.use= c.br, no.legend = T,do.identify = F)

#eminence
TSNEPlot(E18,do.return = F, pt.size = 0.001,label.size=4,do.label = F, group.by="Eminence4",colors.use= c("#9a5ea1","#98823c","lightgray"), no.legend = T,do.identify = F)

#tSNE Meis and Gad
FeaturePlot(E18,c("Meis2"),cols.use=c("lightgrey","brown1"),no.legend = T)
FeaturePlot(E18,c("Gad1"),cols.use=c("lightgrey","brown1"),no.legend = T)
FeaturePlot(P10,c("Gad1"),cols.use=c("lightgrey","brown1"),no.legend = T)
FeaturePlot(P10,c("Meis2"),cols.use=c("lightgrey","brown1"),no.legend = T)
FeaturePlot(E18,c("Meis2"),cols.use=c("lightgrey","brown1"),no.legend = T,do.hover = T)
FeaturePlot(P10,c("Gad1"),cols.use=c("lightgrey","brown1"),no.legend = F)

# P56 adult dataset
P56@meta.data[, "AdultTypes2"] <- factor(P56@meta.data[, "AdultTypes2"], levels = c("Pvalb-1","Pvalb-2","Sst-1","Sst-2","Sst-3","Vip-1","Vip-2","Vip-3","Id2-1","Id2-2","Id2-3","Igfbp6", "Th","Nos1")) 
TSNEPlot(P56,do.return = F, pt.size = 0.001,label.size=4,do.label = T,group.by = "AdultTypes2",colors.use = c.AT2,no.legend = T)
TSNEPlot(P56,do.return = F, pt.size = 0.001,label.size=4,do.label = T,group.by = "AdultTypes3",colors.use =c.AT2,no.legend = F)

P56=SetAllIdent(P56,"AdultTypes2")
P56@ident <- factor(P56@ident, levels = c("Pvalb-1","Pvalb-2","Sst-1","Sst-2","Sst-3","Vip-1","Vip-2","Vip-3","Id2-1","Id2-2","Id2-3","Igfbp6", "Th","Nos1")) 
P56.g = c(rownames(FindMarkers(P56,c("Sst-1"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Sst-2"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Sst-3"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Vip-1"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Vip-2"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Vip-3"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Id2-1"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Id2-2"),only.pos = TRUE))[1:5], 
          rownames(FindMarkers(P56,c("Id2-3"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Igfbp6"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Th"),only.pos = TRUE))[1:5],
          rownames(FindMarkers(P56,c("Nos1"),only.pos = TRUE))[1:5])
P56.gP1a = rownames(FindMarkers(P56,"Pvalb-1",only.pos = TRUE))
P56.gP1=rownames(FindMarkers(P56,"Pvalb-1","Pvalb-2",only.pos = TRUE))
P56.gP2a = rownames(FindMarkers(P56,"Pvalb-2",only.pos = TRUE))
P56.gP2=rownames(FindMarkers(P56,"Pvalb-2","Pvalb-1",only.pos = TRUE))
DoHeatmap(P56,genes.use = c(P56.gP1a[1:5],P56.gP2a[1:5],P56.gP2[1:5],P56.gP1[1:5],P56.g),slim.col.label = TRUE,cex.row=4)#+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm")) 
FeaturePlot(P56,c("Pvalb","Sst","Vip","Id2"),cols.use=c("lightgrey","brown1"),no.legend = T)
FeaturePlot(P56,c("Igfbp6","Th","Nos1"),cols.use=c("lightgrey","brown1"),no.legend = T)



