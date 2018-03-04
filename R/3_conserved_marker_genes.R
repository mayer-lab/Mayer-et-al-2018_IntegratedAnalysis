# this script performs the following steps:
# 1) identification of cardinal type interneuron marker genes that are conserved between embryo and adult at each developmental stage
# 2) determines occurrence of conserved marker genes along a developmental time-course 
# 3) averaging of datasets 
# 3) visualization; in form of "bar-graph-heat-maps"   
# related to Fig. 4E and Extended Figure 7

#find conserved marker genes
ccaE13=SetAllIdent(ccaE13,"Map1")
ccaE18=SetAllIdent(ccaE18,"Map1")
ccaP10=SetAllIdent(ccaP10,"Map1")
#
pval.thresh=0.05/14162
#
E13.Pvalb_Sst.conserved=FindConservedMarkers(ccaE13,"Pvalb","Sst",grouping.var = "DevStage",thresh.use=log(1.25))
E13.Pvalb_Sst.conserved$sign=apply(SubsetColumn(E13.Pvalb_Sst.conserved,"avg_diff"),1,prod)
E13.Pvalb_Sst.conserved.use=subset(E13.Pvalb_Sst.conserved, (max_pval<pval.thresh & sign>0))
E13.Pvalb_Sst = rownames(subset(E13.Pvalb_Sst.conserved.use, E13_avg_diff>log(1.25)))
E13.Sst_Pvalb = rownames(subset(E13.Pvalb_Sst.conserved.use, E13_avg_diff<(-1*log(1.25))))
#
E13.Vip_Id2.conserved=FindConservedMarkers(ccaE13,"Vip","Id2",grouping.var = "DevStage",thresh.use=log(1.25))
E13.Vip_Id2.conserved$sign=apply(SubsetColumn(E13.Vip_Id2.conserved,"avg_diff"),1,prod)
E13.Vip_Id2.conserved.use=subset(E13.Vip_Id2.conserved, (max_pval<pval.thresh & sign>0))
E13.Vip_Id2 = rownames(subset(E13.Vip_Id2.conserved.use, E13_avg_diff>log(1.25)))
E13.Id2_Vip = rownames(subset(E13.Vip_Id2.conserved.use, E13_avg_diff<(-1*log(1.25))))
#
E13.MGE_CGE.conserved=FindConservedMarkers(ccaE13,c("Pvalb","Sst"),c("Vip","Id2"),grouping.var = "DevStage",thresh.use=log(1.25))
E13.MGE_CGE.conserved$sign=apply(SubsetColumn(E13.MGE_CGE.conserved,"avg_diff"),1,prod)
E13.MGE_CGE.conserved.use=subset(E13.MGE_CGE.conserved, (max_pval<pval.thresh & sign>0))
E13.MGE_CGE = rownames(subset(E13.MGE_CGE.conserved.use, E13_avg_diff>log(1.25)))
E13.CGE_MGE = rownames(subset(E13.MGE_CGE.conserved.use, E13_avg_diff<(-1*log(1.25))))
#
pval.thresh=0.05/14162
#
E18.Pvalb_Sst.conserved=FindConservedMarkers(ccaE18,"Pvalb","Sst",grouping.var = "DevStage",thresh.use=log(1.25))
E18.Pvalb_Sst.conserved$sign=apply(SubsetColumn(E18.Pvalb_Sst.conserved,"avg_diff"),1,prod)
E18.Pvalb_Sst.conserved.use=subset(E18.Pvalb_Sst.conserved, (max_pval<pval.thresh & sign>0))
E18.Pvalb_Sst = rownames(subset(E18.Pvalb_Sst.conserved.use, E18_avg_diff>log(1.25)))
E18.Sst_Pvalb = rownames(subset(E18.Pvalb_Sst.conserved.use, E18_avg_diff<(-1*log(1.25))))
#
E18.Vip_Id2.conserved=FindConservedMarkers(ccaE18,"Vip","Id2",grouping.var = "DevStage",thresh.use=log(1.25))
E18.Vip_Id2.conserved$sign=apply(SubsetColumn(E18.Vip_Id2.conserved,"avg_diff"),1,prod)
E18.Vip_Id2.conserved.use=subset(E18.Vip_Id2.conserved, (max_pval<pval.thresh & sign>0))
E18.Vip_Id2 = rownames(subset(E18.Vip_Id2.conserved.use, E18_avg_diff>log(1.25)))
E18.Id2_Vip = rownames(subset(E18.Vip_Id2.conserved.use, E18_avg_diff<(-1*log(1.25))))
#
E18.MGE_CGE.conserved=FindConservedMarkers(ccaE18,c("Pvalb","Sst"),c("Vip","Id2"),grouping.var = "DevStage",thresh.use=log(1.25))
E18.MGE_CGE.conserved$sign=apply(SubsetColumn(E18.MGE_CGE.conserved,"avg_diff"),1,prod)
E18.MGE_CGE.conserved.use=subset(E18.MGE_CGE.conserved, (max_pval<pval.thresh & sign>0))
E18.MGE_CGE = rownames(subset(E18.MGE_CGE.conserved.use, E18_avg_diff>log(1.25)))
E18.CGE_MGE = rownames(subset(E18.MGE_CGE.conserved.use, E18_avg_diff<(-1*log(1.25))))
#
pval.thresh=0.05/14162
#
P10.Pvalb_Sst.conserved=FindConservedMarkers(ccaP10,"Pvalb","Sst",grouping.var = "DevStage",thresh.use=log(1.25))
P10.Pvalb_Sst.conserved$sign=apply(SubsetColumn(P10.Pvalb_Sst.conserved,"avg_diff"),1,prod)
P10.Pvalb_Sst.conserved.use=subset(P10.Pvalb_Sst.conserved, (max_pval<pval.thresh & sign>0))
P10.Pvalb_Sst = rownames(subset(P10.Pvalb_Sst.conserved.use, P10_avg_diff>log(1.25)))
P10.Sst_Pvalb = rownames(subset(P10.Pvalb_Sst.conserved.use, P10_avg_diff<(-1*log(1.25))))
#
P10.Vip_Id2.conserved=FindConservedMarkers(ccaP10,"Vip","Id2",grouping.var = "DevStage",thresh.use=log(1.25))
P10.Vip_Id2.conserved$sign=apply(SubsetColumn(P10.Vip_Id2.conserved,"avg_diff"),1,prod)
P10.Vip_Id2.conserved.use=subset(P10.Vip_Id2.conserved, (max_pval<pval.thresh & sign>0))
P10.Vip_Id2 = rownames(subset(P10.Vip_Id2.conserved.use, P10_avg_diff>log(1.25)))
P10.Id2_Vip = rownames(subset(P10.Vip_Id2.conserved.use, P10_avg_diff<(-1*log(1.25))))
#
P10.MGE_CGE.conserved=FindConservedMarkers(ccaP10,c("Pvalb","Sst"),c("Vip","Id2"),grouping.var = "DevStage",thresh.use=log(1.25))
P10.MGE_CGE.conserved$sign=apply(SubsetColumn(P10.MGE_CGE.conserved,"avg_diff"),1,prod)
P10.MGE_CGE.conserved.use=subset(P10.MGE_CGE.conserved, (max_pval<pval.thresh & sign>0))
P10.MGE_CGE = rownames(subset(P10.MGE_CGE.conserved.use, P10_avg_diff>log(1.25)))
P10.CGE_MGE = rownames(subset(P10.MGE_CGE.conserved.use, P10_avg_diff<(-1*log(1.25))))

# identification of conserved marker genes along a developmental time course
# PV
g.1a = c(P10.Pvalb_Sst) #P10
g.2a = c(E18.Pvalb_Sst) #E18
g.3a = c(E13.Pvalb_Sst) #E13
Pvalb_all = Reduce(intersect, list(g.1a,g.2a,g.3a))
Pvalb_13 = intersect(g.3a,Pvalb_all)
Pvalb_18 = setdiff(intersect(g.2a,g.1a),Pvalb_13)
Pvalb_10 = setdiff(g.1a,union(Pvalb_18,Pvalb_13))
# SST
g.1a = c(P10.Sst_Pvalb) #P10
g.2a = c(E18.Sst_Pvalb) #E18
g.3a = c(E13.Sst_Pvalb) #E13
Sst_all = Reduce(intersect, list(g.1a,g.2a,g.3a))
Sst_13 = intersect(g.3a,Sst_all)
Sst_18 = setdiff(intersect(g.2a,g.1a),Sst_13)
Sst_10 = setdiff(g.1a,union(Sst_18,Sst_13))
# VIP
g.1a = c(P10.Vip_Id2) #P10
g.2a = c(E18.Vip_Id2) #E18
g.3a = c(E13.Vip_Id2) #E13
Vip_all = Reduce(intersect, list(g.1a,g.2a,g.3a))
Vip_13 = intersect(g.3a,Vip_all)
Vip_18 = setdiff(intersect(g.2a,g.1a),Vip_13)
Vip_10 = setdiff(g.1a,union(Vip_18,Vip_13))
# Id2
g.1a = c(P10.Id2_Vip) #P10
g.2a = c(E18.Id2_Vip) #E18
g.3a = c(E13.Id2_Vip) #E13
Id2_all = Reduce(intersect, list(g.1a,g.2a,g.3a))
Id2_13 = intersect(g.3a,Id2_all)
Id2_18 = setdiff(intersect(g.2a,g.1a),Id2_13)
Id2_10 = setdiff(g.1a,union(Id2_18,Id2_13))
# MGE
g.1a = c(P10.MGE_CGE) #P10
g.2a = c(E18.MGE_CGE) #E18
g.3a = c(E13.MGE_CGE) #E13
MGE_all = Reduce(intersect, list(g.1a,g.2a,g.3a))
MGE_13 = intersect(g.3a,MGE_all)
MGE_18 = setdiff(intersect(g.2a,g.1a),MGE_13)
MGE_10 = setdiff(g.1a,union(MGE_18,MGE_13))
# SST
g.1a = c(P10.CGE_MGE) #P10
g.2a = c(E18.CGE_MGE) #E18
g.3a = c(E13.CGE_MGE) #E13
CGE_all = Reduce(intersect, list(g.1a,g.2a,g.3a))
CGE_13 = intersect(g.3a,CGE_all)
CGE_18 = setdiff(intersect(g.2a,g.1a),CGE_13)
CGE_10 = setdiff(g.1a,union(CGE_18,CGE_13))

# averaging datasets
ccaP10=SetIdent(ccaP10, ident.use=ccaP10@meta.data$MapStage)
ccaE18=SetIdent(ccaE18, ident.use=ccaE18@meta.data$MapStage)
ccaE13=SetIdent(ccaE13, ident.use=ccaE13@meta.data$MapStage)
ccaP10.avg=AverageExpression(ccaP10,return.seurat = T)
ccaP10.avg=ScaleData(NormalizeData(ccaP10.avg))
ccaE18.avg=AverageExpression(ccaE18,return.seurat = T)
ccaE18.avg=ScaleData(NormalizeData(ccaE18.avg))
ccaE13.avg=AverageExpression(ccaE13,return.seurat = T)
ccaE13.avg=ScaleData(NormalizeData(ccaE13.avg))
save(ccaP10.avg,ccaE18.avg,ccaE13.avg,file = "~/cca.avg.Rda")

# bar-graph-heat-maps visualization - Extended Data Figure 
#PV SST
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("P10 Pvalb","P10 Sst")),genes.use = c(Pvalb_13,Pvalb_18,Pvalb_10)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_P_10_E.pdf", width=0.5,height=length(c(Pvalb_13,Pvalb_18,Pvalb_10))/20,dpi=600)
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("P10 Pvalb","P10 Sst")),genes.use = rev(c(Sst_13,Sst_18,Sst_10))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_S_10_E.pdf", width=0.5,height=length(c(Sst_13,Sst_18,Sst_10))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("E18 Pvalb","E18 Sst")),genes.use = c(Pvalb_13,Pvalb_18)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_P_18_E.pdf", width=0.5,height=length(c(Pvalb_13,Pvalb_18))/20,dpi=600)
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("E18 Pvalb","E18 Sst")),genes.use = rev(c(Sst_13,Sst_18))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_S_18_E.pdf", width=0.5,height=length(c(Sst_13,Sst_18))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE13.avg,ident.use = c("E13 Pvalb","E13 Sst")),genes.use = c(Pvalb_13)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_P_13_E.pdf", width=0.5,height=length(c(Pvalb_13))/20,dpi=600)
DoHeatmap(SubsetData(ccaE13.avg,ident.use = c("E13 Pvalb","E13 Sst")),genes.use = rev(c(Sst_13))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_S_13_E.pdf", width=0.5,height=length(c(Sst_13))/20,dpi=600)
#VIP IDs
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("P10 Vip","P10 Id2")),genes.use = c(Vip_13,Vip_18,Vip_10)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_V_10_E.pdf", width=0.5,height=length(c(Vip_13,Vip_18,Vip_10))/20,dpi=600)
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("P10 Vip","P10 Id2")),genes.use = rev(c(Id2_18,Id2_10))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_I_10_E.pdf", width=0.5,height=length(c(Id2_13,Id2_18,Id2_10))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("E18 Vip","E18 Id2")),genes.use = c(Vip_13,Vip_18)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_V_18_E.pdf", width=0.5,height=length(c(Vip_13,Vip_18))/20,dpi=600)
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("E18 Vip","E18 Id2")),genes.use = rev(c(Id2_18))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_I_18_E.pdf", width=0.5,height=length(c(Id2_13,Id2_18))/20,dpi=600)
#
#MGE IDs
ccaP10.avg@ident <- plyr::mapvalues(ccaP10.avg@ident, from = c("P10 Pvalb","P10 Sst", "P10 Vip","P10 Id2"), to = c("MGE","MGE","CGE","CGE"))
ccaE18.avg@ident <- plyr::mapvalues(ccaE18.avg@ident, from = c("E18 Pvalb","E18 Sst", "E18 Vip","E18 Id2"), to = c("MGE","MGE","CGE","CGE"))
ccaE13.avg@ident <- plyr::mapvalues(ccaE13.avg@ident, from = c("E13 Pvalb","E13 Sst", "E13 Vip","E13 Id2"), to = c("MGE","MGE","CGE","CGE"))
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("MGE","CGE")),genes.use = c(MGE_13,MGE_18,MGE_10)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_M_10_E.pdf", width=0.5,height=length(c(MGE_13,MGE_18,MGE_10))/20,dpi=600)
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("MGE","CGE")),genes.use = rev(c(CGE_13,CGE_18,CGE_10))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_E_10_E.pdf", width=0.5,height=length(c(CGE_13,CGE_18,CGE_10))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("MGE","CGE")),genes.use = c(MGE_13,MGE_18)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_M_18_E.pdf", width=0.5,height=length(c(MGE_13,MGE_18))/20,dpi=600)
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("MGE","CGE")),genes.use = rev(c(CGE_13,CGE_18))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_E_18_E.pdf", width=0.5,height=length(c(CGE_13,CGE_18))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE13.avg,ident.use = c("MGE","CGE")),genes.use = c(MGE_13)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_M_13_E.pdf", width=0.5,height=length(c(MGE_13))/20,dpi=600)
DoHeatmap(SubsetData(ccaE13.avg,ident.use = c("MGE","CGE")),genes.use = rev(c(CGE_13))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=F)+theme(plot.margin = unit(c(0,-0.2,-1,0.03), "cm"))
ggsave("Heatmap_E_13_E.pdf", width=0.5,height=length(c(CGE_13))/20,dpi=600)
load(file = "~/cca.avg.Rda") 

# bar-graph-heat-maps visualization - Genenames 
#PV SST
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("P10 Pvalb","P10 Sst")),genes.use = c(Pvalb_13,Pvalb_18,Pvalb_10)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_P_10_F.pdf", width=0.5,height=length(c(Pvalb_13,Pvalb_18,Pvalb_10))/20,dpi=600)
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("P10 Pvalb","P10 Sst")),genes.use = rev(c(Sst_13,Sst_18,Sst_10))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_S_10_F.pdf", width=0.5,height=length(c(Sst_13,Sst_18,Sst_10))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("E18 Pvalb","E18 Sst")),genes.use = c(Pvalb_13,Pvalb_18)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_P_18_F.pdf", width=0.5,height=length(c(Pvalb_13,Pvalb_18))/20,dpi=600)
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("E18 Pvalb","E18 Sst")),genes.use = rev(c(Sst_13,Sst_18))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_S_18_F.pdf", width=0.5,height=length(c(Sst_13,Sst_18))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE13.avg,ident.use = c("E13 Pvalb","E13 Sst")),genes.use = c(Pvalb_13)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_P_13_F.pdf", width=0.5,height=length(c(Pvalb_13))/20,dpi=600)
DoHeatmap(SubsetData(ccaE13.avg,ident.use = c("E13 Pvalb","E13 Sst")),genes.use = rev(c(Sst_13))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_S_13_F.pdf", width=0.5,height=length(c(Sst_13))/20,dpi=600)
#VIP IDs
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("P10 Vip","P10 Id2")),genes.use = c(Vip_13,Vip_18,Vip_10)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_V_10_F.pdf", width=0.5,height=length(c(Vip_13,Vip_18,Vip_10))/20,dpi=600)
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("P10 Vip","P10 Id2")),genes.use = rev(c(Id2_18,Id2_10))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_I_10_F.pdf", width=0.5,height=length(c(Id2_13,Id2_18,Id2_10))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("E18 Vip","E18 Id2")),genes.use = c(Vip_13,Vip_18)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_V_18_F.pdf", width=0.5,height=length(c(Vip_13,Vip_18))/20,dpi=600)
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("E18 Vip","E18 Id2")),genes.use = rev(c(Id2_18))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_I_18_F.pdf", width=0.5,height=length(c(Id2_13,Id2_18))/20,dpi=600)
#MGE IDs
load(file = "~/cca.avg.Rda")
ccaP10.avg@ident <- plyr::mapvalues(ccaP10.avg@ident, from = c("P10 Pvalb","P10 Sst", "P10 Vip","P10 Id2"), to = c("MGE","MGE","CGE","CGE"))
ccaE18.avg@ident <- plyr::mapvalues(ccaE18.avg@ident, from = c("E18 Pvalb","E18 Sst", "E18 Vip","E18 Id2"), to = c("MGE","MGE","CGE","CGE"))
ccaE13.avg@ident <- plyr::mapvalues(ccaE13.avg@ident, from = c("E13 Pvalb","E13 Sst", "E13 Vip","E13 Id2"), to = c("MGE","MGE","CGE","CGE"))
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("MGE","CGE")),genes.use = c(MGE_13,MGE_18,MGE_10)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_M_10_F.pdf", width=0.5,height=length(c(MGE_13,MGE_18,MGE_10))/20,dpi=600)
DoHeatmap(SubsetData(ccaP10.avg,ident.use = c("MGE","CGE")),genes.use = rev(c(CGE_13,CGE_18,CGE_10))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_F_10_F.pdf", width=0.5,height=length(c(CGE_13,CGE_18,CGE_10))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("MGE","CGE")),genes.use = c(MGE_13,MGE_18)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_M_18_F.pdf", width=0.5,height=length(c(MGE_13,MGE_18))/20,dpi=600)
DoHeatmap(SubsetData(ccaE18.avg,ident.use = c("MGE","CGE")),genes.use = rev(c(CGE_13,CGE_18))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_F_18_F.pdf", width=0.5,height=length(c(CGE_13,CGE_18))/20,dpi=600)
#
DoHeatmap(SubsetData(ccaE13.avg,ident.use = c("MGE","CGE")),genes.use = c(MGE_13)
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_M_13_F.pdf", width=0.5,height=length(c(MGE_13))/20,dpi=600)
DoHeatmap(SubsetData(ccaE13.avg,ident.use = c("MGE","CGE")),genes.use = rev(c(CGE_13))
          ,slim.col.label = TRUE,group.by="ident",remove.key=T,cex.row=4)+theme(plot.margin = unit(c(0,0.4,-1,0.03), "cm"))
ggsave("Heatmap_F_13_F.pdf", width=0.5,height=length(c(CGE_13))/20,dpi=600)
load(file = "~/cca.avg.Rda")

