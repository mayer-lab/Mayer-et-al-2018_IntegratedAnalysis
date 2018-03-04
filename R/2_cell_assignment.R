# this script performs the following steps:
# 1) assignment of precursor cells to adult subtypes
# 3) visualization of assignment - cardinal types 
# 4) visualization of assignment - subtypes
# related to Fig.4 B-D) lower part and Extended Data Figure 9

library(pbapply)
# assign precursor cells to adult subtypes
MapCells=function(object, timepoint.dev="P10", timepoint.end="P56", num.k=10, thresh.require=0.9,map.col="AdultTypes1",new.col="Map1") {
  tsne.dims=object@calc.params$RunTSNE$dims.use
  input.dims=GetCellEmbeddings(object,reduction.type = "cca.aligned",dims.use = tsne.dims)
  input.dist=as.matrix(dist(input.dims))
  cells.end=FastWhichCells(object,"DevStage",timepoint.end)
  cells.map=FastWhichCells(object,"DevStage",timepoint.dev)
  map.id=pbsapply(cells.map,function(x)
    names(which(sort(table(object@meta.data[names(sort(input.dist[x,cells.end])[1:num.k]),map.col]),decreasing = T)>=num.k*thresh.require))[1])
  nn.same=pbsapply(cells.map, function(x) length(intersect(names(sort(input.dist[x,])[1:num.k]),cells.map)))
  print(table(map.id))
  map.id[is.na(map.id)]="Unassigned"
  map.id[names(which(nn.same>=num.k*1))]="Unassigned"
  print(table(map.id))
  object@meta.data[cells.end,new.col]=object@meta.data[cells.end,map.col]
  object@meta.data[cells.map,new.col]=map.id
  return(object)
}
ccaP10=MapCells(ccaP10,timepoint.dev = "P10",map.col = "AdultTypes1",new.col = "Map1",thresh.require = 0.9,num.k = 10)
ccaP10=MapCells(ccaP10,timepoint.dev = "P10",map.col = "AdultTypes3",new.col = "Map3",thresh.require = 0.9,num.k = 10)
ccaE18=MapCells(ccaE18,timepoint.dev = "E18",map.col = "AdultTypes1",new.col = "Map1",thresh.require = 0.9,num.k = 10)
ccaE18=MapCells(ccaE18,timepoint.dev = "E18",map.col = "AdultTypes3",new.col = "Map3",thresh.require = 0.9,num.k = 10)
ccaE13=MapCells(ccaE13,timepoint.dev = "E13",map.col = "AdultTypes1",new.col = "Map1",thresh.require = 0.9,num.k = 10)
ccaE13=MapCells(ccaE13,timepoint.dev = "E13",map.col = "AdultTypes3",new.col = "Map3",thresh.require = 0.9,num.k = 10)

# visualization of assignment - (A) cardinal types 
#p10
ccaP10@meta.data$MapStage=paste(ccaP10@meta.data$DevStage,ccaP10@meta.data$Map1,sep=" ")
temp1 = rownames(table(ccaP10@meta.data$MapStage))
ccaP10@meta.data[, "MapStage1"] <- plyr::mapvalues(ccaP10@meta.data$MapStage, from = c(temp1), to = c(temp1[1:8],"P56","P56","P56","P56","P56","P56","P56"))
temp2 = rownames(table(ccaP10@meta.data$MapStage1))
ccaP10@meta.data[, "MapStage1"] <- factor(ccaP10@meta.data[, "MapStage1"], levels = temp2[c(4,5,8,1,2,3,6,7,9)])
ccaP10_Mapped=SetAllIdent(ccaP10,"MapStage1")
ccaP10_Mapped=SubsetData(ccaP10_Mapped,ident.remove = c("P56","P10 Unassigned"))
TSNEPlot(ccaP10_Mapped, do.return = F, pt.size = 0.1,label.size=4, do.label = T, group.by="MapStage1",colors.use=c.AT1,no.legend = F )
#e18
ccaE18@meta.data$MapStage=paste(ccaE18@meta.data$DevStage,ccaE18@meta.data$Map1,sep=" ")
ccaE18@meta.data$test2=paste(ccaE18@meta.data$DevStage,ccaE18@meta.data$Map3,sep=" ")
#
temp1 = rownames(table(ccaE18@meta.data$MapStage))
ccaE18@meta.data[, "MapStage1"] <- plyr::mapvalues(ccaE18@meta.data$MapStage, from = c(temp1), to = c(temp1[1:8],"P56","P56","P56","P56","P56","P56","P56"))
temp2 = rownames(table(ccaE18@meta.data$MapStage1))
ccaE18@meta.data[, "MapStage1"] <- factor(ccaE18@meta.data[, "MapStage1"], levels = temp2[c(4,5,8,1,2,3,6,7,9)])
ccaE18_Mapped=SetAllIdent(ccaE18,"MapStage1")
ccaE18_Mapped=SubsetData(ccaE18_Mapped,ident.remove = c("P56","E18 Unassigned"))
TSNEPlot(ccaE18_Mapped, do.return = F, pt.size = 0.1,label.size=4, do.label = T, group.by="MapStage1",colors.use=c.AT1,no.legend = F )
#e13
ccaE13@meta.data$MapStage=paste(ccaE13@meta.data$DevStage,ccaE13@meta.data$Map1,sep=" ")
temp1 = rownames(table(ccaE13@meta.data$MapStage))
ccaE13@meta.data[, "MapStage1"] <- plyr::mapvalues(ccaE13@meta.data$MapStage, from = c(temp1), to = c(temp1[1:7],"P56","P56","P56","P56","P56","P56","P56"))
temp2 = rownames(table(ccaE13@meta.data$MapStage1))
ccaE13@meta.data[, "MapStage1"] <- factor(ccaE13@meta.data[, "MapStage1"], levels = temp2[c(3,4,7,1,2,5,6,8)])
ccaE13_Mapped=SetAllIdent(ccaE13,"MapStage1")
ccaE13_Mapped=SubsetData(ccaE13_Mapped,ident.remove = c("P56","E13 Unassigned"))
TSNEPlot(ccaE13_Mapped, do.return = F, pt.size = 0.1,label.size=4, do.label = T, group.by="MapStage1",colors.use=c.AT1,no.legend = F )

# visualization of assignment - (B) subtypes 
#p10
p10=SubsetData(ccaP10,cells.use = FastWhichCells(ccaP10,group.by = "DevStage",subset.value = "P10"))
map3.ids=names(table(p10@meta.data$Map3))
p10=SetAllIdent(p10,"Map3")
ccaP10_Mapped=SubsetData(p10,ident.remove = c("Unassigned"))
ccaP10_Mapped@ident <- factor(ccaP10_Mapped@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
table(ccaP10_Mapped@ident)
#
TSNEPlot(ccaP10_Mapped, do.return = F, pt.size = 1,label.size=4, do.label = T,no.legend = F,colors.use=c.AT2)
TSNEPlot(SubsetData(ccaP10_Mapped,ident.use = c("Pvalb-1 Basket","Pvalb-2 Chandelier")), do.return = F, pt.size = 1,label.size=4, do.label = F,no.legend = F )
TSNEPlot(SubsetData(ccaP10_Mapped,ident.use = c("Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti")), do.return = F, pt.size = 1,label.size=4, do.label = F, no.legend = F )
TSNEPlot(SubsetData(ccaP10_Mapped,ident.use = c("Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)")), do.return = F, pt.size = 1,label.size=4, do.label = F, no.legend = F )
TSNEPlot(SubsetData(ccaP10_Mapped,ident.use = c("Id2-3 Sncg","Id2-2 NGFC","Id2-1 NDNF")), do.return = F, pt.size = 1,label.size=4, do.label = F, no.legend = F )
#e18
e18=SubsetData(ccaE18,cells.use = FastWhichCells(ccaE18,group.by = "DevStage",subset.value = "E18"))
map3.ids=names(table(e18@meta.data$Map3))
e18=SetAllIdent(e18,"Map3")
ccaE18_Mapped=SubsetData(e18,ident.remove = c("Unassigned"))
ccaE18_Mapped@ident <- factor(ccaE18_Mapped@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
TSNEPlot(ccaE18_Mapped, do.return = F, pt.size = 1,label.size=4, do.label = F,no.legend = F,colors.use=c.AT2)
TSNEPlot(SubsetData(ccaE18_Mapped,ident.use = c("Pvalb-1 Basket","Pvalb-2 Chandelier")), do.return = F, pt.size = 1,label.size=4, do.label = F,no.legend = F )
TSNEPlot(SubsetData(ccaE18_Mapped,ident.use = c("Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti")), do.return = F, pt.size = 1,label.size=4, do.label = F, no.legend = F )
TSNEPlot(SubsetData(ccaE18_Mapped,ident.use = c("Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)")), do.return = F, pt.size = 1,label.size=4, do.label = F, no.legend = F )
TSNEPlot(SubsetData(ccaE18_Mapped,ident.use = c("Id2-3 Sncg","Id2-2 NGFC","Id2-1 NDNF")), do.return = F, pt.size = 1,label.size=4, do.label = F, no.legend = F )

