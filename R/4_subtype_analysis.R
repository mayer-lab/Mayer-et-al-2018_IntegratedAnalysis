# this script performs the subtype analysis related to Extended Data Figure 9 

library(pbapply)
returnTopX=function(data, sort_col, n_return) {
  sort_ids=as.character(unique(data[,sort_col]))
  to_return=data.frame()
  for(i in sort_ids) {
    data_sub=data[data[,sort_col]==i,]
    to_return=rbind(to_return,head(data_sub,n_return))
  }
  return(to_return)
}

#P10
pval.thresh=0.05/nrow(ccaP10@data)
#sst
p10.sst=SubsetData(p10,ident.use = grep("Sst",map3.ids,value = T))
p10.sst@ident <- factor(p10.sst@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
p10.markers.sst=FindAllMarkers(p10.sst,thresh.use = log(2),only.pos = T)
p10.markers.sst.pass=subset(p10.markers.sst,p_val<pval.thresh)
p10.markers.sst.plot=returnTopX(data=p10.markers.sst.pass,sort_col = "cluster",20)
DoHeatmap(p10.sst,genes.use = unique(p10.markers.sst.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(p10.sst@cell.names),remove.key=T,cex.row=5)
#pv
p10.pv=SubsetData(p10,ident.use = grep("Pvalb",map3.ids,value = T))
p10.pv@ident <- factor(p10.pv@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
p10.markers.pv=FindAllMarkers(p10.pv,thresh.use = log(2),only.pos = T)
p10.markers.pv.pass=subset(p10.markers.pv,p_val<pval.thresh)
p10.markers.pv.plot=returnTopX(data=p10.markers.pv.pass,sort_col = "cluster",20)
DoHeatmap(p10.pv,genes.use = unique(p10.markers.pv.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(p10.pv@cell.names),remove.key=T,cex.row=5)
#Id
p10.id=SubsetData(p10,ident.use = grep("Id2",map3.ids,value = T))
p10.id@ident <- factor(p10.id@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
p10.markers.id=FindAllMarkers(p10.id,thresh.use = log(2),only.pos = T)
p10.markers.id.plot=returnTopX(data=p10.markers.id,sort_col = "cluster",20)
DoHeatmap(p10.id,genes.use = unique(p10.markers.id.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(p10.id@cell.names),remove.key=T,cex.row=5)
#vip
p10.vip=SubsetData(p10,ident.use = grep("Vip",map3.ids,value = T))
p10.vip@ident <- factor(p10.vip@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
p10.markers.vip=FindAllMarkers(p10.vip,thresh.use = log(2),only.pos = T)
p10.markers.vip.plot=returnTopX(data=p10.markers.vip,sort_col = "cluster",20)
DoHeatmap(p10.vip,genes.use = unique(p10.markers.vip.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(p10.vip@cell.names),remove.key=T,cex.row=5)

#E18
pval.thresh=0.05/nrow(ccaE18@data)
#sst
e18.sst=SubsetData(e18,ident.use = grep("Sst",map3.ids,value = T))
e18.sst@ident <- factor(e18.sst@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
e18.markers.sst=FindAllMarkers(e18.sst,thresh.use = log(2),only.pos = T)
e18.markers.sst.pass=subset(e18.markers.sst,p_val<pval.thresh)
e18.markers.sst.plot=returnTopX(data=e18.markers.sst.pass,sort_col = "cluster",20)
DoHeatmap(e18.sst,genes.use = unique(e18.markers.sst.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(e18.sst@cell.names),remove.key=T,cex.row=5)
#pv
e18.pv=SubsetData(e18,ident.use = grep("Pvalb",map3.ids,value = T))
e18.pv@ident <- factor(e18.pv@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
e18.markers.pv=FindAllMarkers(e18.pv,thresh.use = log(2),only.pos = T)
e18.markers.pv.pass=subset(e18.markers.pv,p_val<pval.thresh)
e18.markers.pv.plot=returnTopX(data=e18.markers.pv.pass,sort_col = "cluster",20)
DoHeatmap(e18.pv,genes.use = unique(e18.markers.pv.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(e18.pv@cell.names),remove.key=T,cex.row=5)
#id2
e18.id=SubsetData(e18,ident.use = grep("Id2",map3.ids,value = T))
e18.id@ident <- factor(e18.id@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
e18.markers.id=FindAllMarkers(e18.id,thresh.use = log(2),only.pos = T)
e18.markers.id.plot=returnTopX(data=e18.markers.id,sort_col = "cluster",20)
DoHeatmap(e18.id,genes.use = unique(e18.markers.id.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(e18.id@cell.names),remove.key=T,cex.row=5)

#vip
e18.vip=SubsetData(e18,ident.use = grep("Vip",map3.ids,value = T))
e18.vip@ident <- factor(e18.vip@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
e18.markers.vip=FindAllMarkers(e18.vip,thresh.use = log(1.5),only.pos = T)
e18.markers.vip.plot=returnTopX(data=e18.markers.vip,sort_col = "cluster",20)
DoHeatmap(e18.vip,genes.use = unique(e18.markers.vip.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(e18.vip@cell.names),remove.key=T,cex.row=5)


#E13  
pval.thresh=0.05/nrow(ccaE13@data)
e13=SubsetData(ccaE13,cells.use = FastWhichCells(ccaE13,group.by = "DevStage",subset.value = "E13"))
map3.ids=names(table(e13@meta.data$Map3))
e13=SetAllIdent(e13,"Map3")
#sst
e13.sst=SubsetData(e13,ident.use = grep("Sst",map3.ids,value = T))
e13.sst@ident <- factor(e13.sst@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
e13.markers.sst=FindAllMarkers(e13.sst,thresh.use = log(2),only.pos = T)
e13.markers.sst.pass=subset(e13.markers.sst,p_val<pval.thresh)
e13.markers.sst.plot=returnTopX(data=e13.markers.sst.pass,sort_col = "cluster",20)
DoHeatmap(e13.sst,genes.use = unique(e13.markers.sst.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(e13.sst@cell.names),remove.key=T,cex.row=5)
#pv
e13.pv=SubsetData(e13,ident.use = grep("Pvalb",map3.ids,value = T))
e13.pv@ident <- factor(e13.pv@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
e13.markers.pv=FindAllMarkers(e13.pv,thresh.use = log(2),only.pos = T)
e13.markers.pv.pass=subset(e13.markers.pv,p_val<pval.thresh)
e13.markers.pv.plot=returnTopX(data=e13.markers.pv.pass,sort_col = "cluster",20)
DoHeatmap(e13.pv,genes.use = unique(e13.markers.pv.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(e13.pv@cell.names),remove.key=T,cex.row=5)
#id2
e13.id=SubsetData(e13,ident.use = grep("Id2",map3.ids,value = T))
e13.id@ident <- factor(e13.id@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
e13.markers.id=FindAllMarkers(e13.id,thresh.use = log(2),only.pos = T)
e13.markers.id.plot=returnTopX(data=e13.markers.id,sort_col = "cluster",20)
DoHeatmap(e13.id,genes.use = unique(e13.markers.id.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(e13.id@cell.names),remove.key=T,cex.row=5)
#vip
e13.vip=SubsetData(e13,ident.use = grep("Vip",map3.ids,value = T))
e13.vip@ident <- factor(e13.vip@ident, levels = c("Pvalb-1 Basket","Pvalb-2 Chandelier","Sst-1 Martinotti","Sst-2 Martinotti","Sst-3 Non-Martinotti","Vip-1 Bipolar","Vip-2 Bipolar","Vip-3 Multipolar (Cck)","Id2-1 NDNF","Id2-2 NGFC","Id2-3 Sncg","Igfbp6","Nos1","Th","Unassigned"))
e13.markers.vip=FindAllMarkers(e13.vip,thresh.use = log(2),only.pos = T)
e13.markers.vip.plot=returnTopX(data=e13.markers.vip,sort_col = "cluster",20)
DoHeatmap(e13.vip,genes.use = unique(e13.markers.vip.plot$gene),slim.col.label = T,group.label.rot = F,cells.use = sample(e13.vip@cell.names),remove.key=T,cex.row=5)
