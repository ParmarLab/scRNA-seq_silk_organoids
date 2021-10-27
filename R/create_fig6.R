##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Petter Storm
##' @export
create_fig6 <- function() {

  
  #######################################################
  # Figure 6 Single cell transcriptomics mapping DA diversity in VM organoids
  ######################################################
  
  #######################################################
  # Libraries and read raw data
  #######################################################
  
  
  
  data.seurat<-readRDS("raw_data/data.seurat.qc.std.SILK.3groups.RDS")
  Idents(data.seurat)<-data.seurat$orig.ident
  data.seurat$treatment<-factor(data.seurat$treatment,levels = c("std","silk","Lam"))
  data.seurat.x3 <- subset(data.seurat,idents=c("03standardorgday30","05standardorgday30","SilkDay30-09silk"),invert=T)
  
  # Subsample to 5000 cells per treatment
  set.seed(28)
  data.seurat.x3@meta.data$cell_id <- rownames(data.seurat.x3@meta.data)
  new_df <- data.seurat.x3@meta.data %>% group_by(treatment) %>% sample_n(5000) 
  
  data.seurat.qc.std.ss <- subset(data.seurat.x3,cells=new_df$cell_id)
  
  
  data.seurat.qc.std.ss$UMAP1 <- Embeddings(data.seurat.qc.std.ss,reduction = "umap")[,1]
  data.seurat.qc.std.ss$UMAP2 <- Embeddings(data.seurat.qc.std.ss,reduction = "umap")[,2]
  
  #######################################################
  # A umap
  #######################################################

  DimPlot(data.seurat.qc.std.ss[],group.by = "predicted.id",split.by = "treatment",ncol = 3,pt.size = .3)+scale_color_manual(values = palettes$`Tableau 10`$value[c(5,2,4,7,1,6,8)])+NoLegend()
  ggsave("results//6A1.pdf",w=12)
  
  #######################################################
  # B barplots
  #######################################################
  
  
  t1 <-table(data.seurat.qc.std.ss$orig.ident,data.seurat.qc.std.ss$celltypes2)
  
  pt1<- data.frame(prop.table(t1,1),group=factor(rep(c("std","silk","silklam"),each=3),levels = c("std","silk","silklam")))
  pt1$Var2<-factor(pt1$Var2, levels=c("FP3", "FP2" ,"FP1" ,"Dopamine neurons" ,"Astrocytes" ,"VLMC", "Neural stem cells"))
  
  ggplot(pt1,aes(x=Var1,y=100*Freq, fill=Var2))+geom_bar(stat="identity")+
    facet_wrap(~group,scale="free",ncol = 3)+scale_fill_manual(values = palettes$`Tableau 10`$value[c(1,7,4,2,3,8,6)])+ theme(legend.position="bottom")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ theme(legend.title = element_blank())+ylab("Cells in cluster (%)")+theme_cowplot()
  ggsave("results//6B.pdf",w=10,h=7)
  
  #######################################################
  # C Chord
  #######################################################
  
  pdf("results/6C.choordplot.silkd.pdf")
  ClusterIdentityChordPlot(data.seurat.qc.std.ss$treatment,data.seurat$celltypes2)
  dev.off()

  
  #######################################################
  # D
  #######################################################
  
  m1<-reshape2::melt(prop.table(table(data.seurat.x3$predicted.id,factor(data.seurat.x3$orig.ident)),2))
  
  m1$g<-factor(c(rep("std",21),rep("silk",21), rep("silkLam",21)),levels=c("std","silk","silkLam"))
  m1$value<-m1$value*100
  
  ggplot(m1[m1$Var1=="Dopamine neurons",],aes(x=g,y=value,fill=g))+geom_violin()+ylab("Percent cells in cluster")+geom_point()+
    scale_fill_manual(values = palettes$`Tableau 10`$value[c(3,2,4,7,1,6,5,8)])+theme_cowplot()+xlab("")+NoLegend()+facet_wrap(~Var1,scales="free")
  ggsave("results//6D.violin.pdf",w=5,h=5)
  

  #######################################################
  # E
  ######################################################
  
  genes <- c("DCX","STMN2","PBX1","TH")
  
  
  md1 <- data.frame(data.seurat.x3@meta.data,
                    t(data.frame(GetAssay(data.seurat.x3,assay = "RNA")[c(intersect(rownames(data.seurat.x3),genes)),])))
 
  
  md1m <- reshape2::melt(md1[,c("treatment",genes)])
  means<-cbind(aggregate(value~variable+treatment,md1m,mean),sem=aggregate(value~variable+treatment,md1m,se)[,"value"])
  
  ggplot(means,aes(x=treatment,y=value,fill=treatment))+geom_bar(stat = "identity")+facet_wrap(~variable,scales="free")+xlab("")+
    geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.2,
                  position=position_dodge(.9))+NoLegend()+
    ylab("Average expression (all cells)\n Scaled counts")+theme_cowplot()
    
  ggsave("results//6E.pdf")
  
  
  
  #######################################################
  # SILK DAY 60/120
  ######################################################
  
  silk.d60.seurat<-readRDS("raw_data/silk.d60.seurat.rds")
  
  
  Idents(silk.d60.seurat)<-silk.d60.seurat$orig.ident
  
  silk.d60.seurat<-RenameIdents(silk.d60.seurat, "17silkorgday60hg38"  = "silk_17_d60",
                                "18silkOrgday60"     = "silk_18_d60",
                                "19silkOrgday60"      = "silk_19_d60",
                                "20silkorglam111day60hg38" = "silk_20_d60",
                                "2silknewprotocolday60"  = "silk_2_d120",
                                "3silknewprotocolday60"  = "silk_3_d120",
                                "4silknewprotocolday60"    = "silk_4_d120",
                                "6silknewprotocolday60"   = "silk_6_d120",
                                "7silknewprotocolday60"    = "silk_7_d120",
                                "8silknewprotocolday60"  = "silk_8_d120"
  ) 
  
  
  
  #######################################################
  # Compare silk and std DA to fetal
  #######################################################
  
  cult.reclust_dopa<-readRDS("../data.seurat.fetal.rds")
  
  fetal.ae<-AverageExpression(cult.reclust_dopa)$RNA
  colnames(fetal.ae)<-paste0("fetal_",colnames(fetal.ae))
  std.ae<-AverageExpression(data.seurat.qc.std.dop.high_5)$RNA
  colnames(std.ae)<-paste0("std_",colnames(std.ae))
  silk.da <- subset(silk.d60.seurat,subset =  celltypes2=="Dopamine neurons")
  
  # label transfer to silk
  anchors <- FindTransferAnchors(reference = data.seurat.qc.std.dop.high_5, query = silk.da,dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = data.seurat.qc.std.dop.high_5$DA_type5, dims = 1:30)
  silk.da <- AddMetaData(silk.da, metadata = predictions)
  
  Idents(silk.da)<-silk.da$predicted.id
  silk.ae<-AverageExpression(silk.da)$RNA
  colnames(silk.ae)<-paste0("silk_",colnames(std.ae))
  
  is.genes <- intersect(rownames(fetal.ae),intersect(rownames(std.ae),rownames(silk.ae)))
  
  ae.mtx <- cbind(fetal.ae[is.genes,],silk.ae[is.genes,],std.ae[is.genes,])
  
  plotMDS(ae.mtx)
  
  
  #######################################################
  # UMAPS
  #######################################################
  
  
  if(file.exists("raw_data/std_silk.merge.rds")){
    std_silk.merge<-readRDS("raw_data/std_silk.merge.rds")
  }else{
    # PREPARE STD
    data.seurat.qc.std <- readRDS("raw_data/data.seurat.qc.std.20201007.rds")
    
    data.seurat.qc.std$CellID <- colnames(data.seurat.qc.std)
    subsample <- data.seurat.qc.std@meta.data[grepl("60|120",data.seurat.qc.std$groups),]
    subsample <- subsample %>% group_by(orig.ident) %>% sample_n(size = 1800)
    data.seurat.qc.std.d60_120.ss <-data.seurat.qc.std[,subsample$CellID]
    
    # PREPARE SILK
    Idents(silk.d60.seurat)<-silk.d60.seurat$orig.ident
    silk.d60.seurat<-RenameIdents(silk.d60.seurat, "17silkorgday60hg38"  = "silk_17_d60",
                                  "18silkOrgday60"     = "silk_18_d60",
                                  "19silkOrgday60"      = "silk_19_d60",
                                  "20silkorglam111day60hg38" = "silk_20_d60",
                                  "2silknewprotocolday60"  = "silk_2_d120",
                                  "3silknewprotocolday60"  = "silk_3_d120",
                                  "4silknewprotocolday60"    = "silk_4_d120",
                                  "6silknewprotocolday60"   = "silk_6_d120",
                                  "7silknewprotocolday60"    = "silk_7_d120",
                                  "8silknewprotocolday60"  = "silk_8_d120"
                                  ) 
    silk.d60.seurat$orig.ident<-Idents(silk.d60.seurat)
    silk.d60.seurat$CellID <- colnames(silk.d60.seurat)
    subsample.silk <- silk.d60.seurat@meta.data %>% group_by(orig.ident) %>% sample_n(size = 1800)
    silk.d60.seurat.ss <- silk.d60.seurat[,subsample.silk$CellID]
    silk.d60.seurat.ss$groups <- "day060"
    silk.d60.seurat.ss@meta.data[grep("d120",Idents(silk.d60.seurat.ss)),"groups"]<-"day120"
    
    silk.d60.seurat.ss$orig.ident<-Idents(silk.d60.seurat.ss)
    # MERGE
    std_silk.merge <- merge(silk.d60.seurat.ss,data.seurat.qc.std.d60_120.ss)
    
    std_silk.merge <- Seurat::NormalizeData(std_silk.merge,verbose = FALSE)
    std_silk.merge <- FindVariableFeatures(std_silk.merge,selection.method = "vst", nfeatures = 2000)
    std_silk.merge <-  ScaleData(std_silk.merge,verbose = FALSE)
    std_silk.merge <-  RunPCA(std_silk.merge,pc.genes = data.seurat@var.genes, npcs = 60, verbose = FALSE)
    
    std_silk.merge <- std_silk.merge %>%
      RunHarmony("orig.ident", plot_convergence = TRUE)
    
    #
    #
    std_silk.merge <- std_silk.merge %>%
      RunUMAP(reduction = "harmony", dims = 1:20) %>%
      FindNeighbors(reduction = "harmony", dims = 1:20) %>%
      FindClusters(resolution = 0.1) %>%
      identity()
    
    
    std_silk.merge$treat <- "std"
    std_silk.merge@meta.data[grep("silk",std_silk.merge$orig.ident),"treat"]<-"silk"
    
    std_silk.merge$treat_group <- paste(std_silk.merge$treat,std_silk.merge$groups)
  }
  
  
  DimPlot(std_silk.merge[,],ncol = 1,group.by = "celltypes2")+
    scale_color_manual(values = palettes$`Tableau 10`$value[c(3,2,4,7,1,6,5,8)])+ggtitle("")
  ggsave("results/silkd60_120.B.UMAP.allcells.pdf",w=8,h=6)
  
  
  
  
  ggplot(pt1,aes(x=tg,y=100*Freq, fill=Var2))+geom_violin()+ggbeeswarm::geom_quasirandom()+
    facet_wrap(~Var2,scale="free")+scale_fill_manual(values = palettes$`Tableau 10`$value[c(4,1,7,2,3,8,6,5)])+ theme(legend.position="bottom")+
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+ theme(legend.title = element_blank())+ylab("Cells in cluster (%)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave("results/silkd60_120.C1.vln.pdf")
  
  ggplot(pt1,aes(x=tg,y=100*Freq, fill=Var2))+geom_boxplot()+
    facet_wrap(~Var2,scale="free")+scale_fill_manual(values = palettes$`Tableau 10`$value[c(4,1,7,2,3,8,6,5)])+ theme(legend.position="bottom")+
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+ theme(legend.title = element_blank())+ylab("Cells in cluster (%)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave("results/silkd60_120.C2.BP.pdf")
  
  ggplot(pt1[grep("120",pt1$Var1),],aes(x=tg,y=100*Freq, fill=Var2))+geom_boxplot()+
    facet_wrap(~Var2,scale="free")+scale_fill_manual(values = palettes$`Tableau 10`$value[c(4,1,7,2,3,8,6,5)])+ theme(legend.position="bottom")+
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+ theme(legend.title = element_blank())+ylab("Cells in cluster (%)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave("results/silkd60_120.C2.BP_only120.pdf")
  
  
  pt1.120 <- pt1[grep("120",pt1$Var1),]
  pt1.120$group<-factor(pt1.120$group,levels=c("std","silk"))
  
  ggplot(pt1.120[pt1.120$Var2=="Dopamine neurons",],aes(x=group,y=Freq*100,fill=group))+geom_violin(trim = F)+ylab("Percent cells in cluster")+ggbeeswarm::geom_beeswarm()+
    scale_fill_manual(values = palettes$`Tableau 10`$value[c(3,2,4,7,1,6,5,8)])+theme_cowplot()+xlab("")+NoLegend()+facet_wrap(~Var2,scales="free")+ylim(c(0,30))
  ggsave("results//6X.danueroncluster_d120_violin.v2.pdf",w=5,h=5)
  
  
  means<-cbind(aggregate(Freq~Var2+tg,pt1,mean),sem=aggregate(Freq~Var2+tg,pt1,sd)[,"Freq"])
  means <- means[means$Var2=="Dopamine neurons",]
  means <- means[grep("120",means$tg),]
  means$Freq<-means$Freq*100
  means$sem<-means$sem*100
  #p1<-DimPlot(data.seurat.qc.std.dop, reduction = "umap", group.by = "DA_type", pt.size = .1, split.by = 'groups')
  ggplot(means,aes(x=tg,y=Freq,fill=tg))+geom_bar(stat = "identity")+facet_wrap(~Var2,scales="free")+xlab("")+
    geom_errorbar(aes(ymax=Freq+sem,ymin=Freq), width=.2,
                  position=position_dodge(.9))+NoLegend()+
    ylab("Percent cells in cluster")+theme_cowplot()+
    scale_fill_manual(values = palettes$`Tableau 10`$value[c(3,2,4,7,1,6,5,8)])+NoLegend()
  
  ggsave("results//6X.bars.pdf",w=5,h=5)
  
  
  
  #d Intraclass correlation
  #silk
  std_silk.merge$orig.ident<-factor(std_silk.merge$orig.ident)
  t2 <-table(factor(std_silk.merge$orig.ident),std_silk.merge$predicted.id)[-c(1:14,16),]
  pheatmap(cor(t(t2)),color = redblue(200))
  icc(t(t2),model="twoway", type="agreement",unit="single")
  
  #std
  data.seurat.qc.std.day120 <- subset(data.seurat.qc.std,subset = groups=="day120")
  t2 <-table(factor(data.seurat.qc.std.day120$orig.ident),data.seurat.qc.std.day120$celltypes2)
  icc(t(t2),model="twoway", type="agreement",unit="single")
  
  
  
  
  genelist <-read.csv("raw_data/genelist_paper.2021505.silk.csv",header=T)
  Idents(std_silk.merge)<-std_silk.merge$celltypes2
  dpd<-DotPlot(std_silk.merge,features = unique(genelist[,"Gene"]))
  dpd.d<-merge(dpd$data,genelist,by.x="features.plot",by.y="Gene")
  
  dpd.d$plot_labels <- factor(substr(dpd.d$id,1,5),
                              levels = c("Neura","FP1","FP2","FP3","Dopam","Astro","Oligo","VLMC"))
  # Finally plot
  
  pal<-tableau_color_pal(palette = "Tableau 10")(10)
  
  dpd.d$features.plot2 <- factor(dpd.d$features.plot,levels=rev(levels(dpd.d$features.plot)))
  
  ggplot(dpd.d,aes(x=plot_labels,y=features.plot2,size=pct.exp,col=Cluster,alpha=avg.exp.scaled))+
    geom_point()+ scale_alpha(range = c(-1, 2))+
    xlab("") + labs(size = "Percent expressed",col="Cell type")+theme_cowplot()+NoLegend()+
    scale_color_manual(values = pal[c(3,2,1,4,5,8)])+ylab("")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave("results/silkd60_120.D.dotplot.pdf",w=4,h=7)
  
  # DA
  
  std_silk.merge.da <- subset(std_silk.merge,subset = celltypes2=="Dopamine neurons")
  table(std_silk.merge.da$orig.ident,std_silk.merge.da$celltypes2)
  std_silk.merge.da$THEXP <- GetAssayData(std_silk.merge.da)["TH",]
  #std_silk.merge.da <- subset(std_silk.merge.da,subset = THEXP>0)
  
  Idents(std_silk.merge.da)<-std_silk.merge.da$treat_group
  std_silk.merge.da <- std_silk.merge.da[,grep("120",std_silk.merge.da$orig.ident)]

  genes <- c("TH","SLC6A3","SLC18A2","ALDH1A1", "PITX3", "CALB1","KCNJ6","KCNC2", "SCN2A", "SCN1A", "KCNQ2", "KCNA5")
  
  
  md1 <- data.frame(std_silk.merge.da@meta.data,
                    t(data.frame(GetAssay(std_silk.merge.da,assay = "RNA")[c(intersect(rownames(std_silk.merge.da),genes)),])))

  # TH high/low
  
  md1m <- reshape2::melt(md1[,c("treat_group",genes)])
  means<-cbind(aggregate(value~variable+treat_group,md1m,mean),sem=aggregate(value~variable+treat_group,md1m,se)[,"value"])
  
  summary(lm(value~variable*treat_group,md1m))
  
  
  
  #p1<-DimPlot(data.seurat.qc.std.dop, reduction = "umap", group.by = "DA_type", pt.size = .1, split.by = 'groups')
  ggplot(means,aes(x=treat_group,y=value,fill=treat_group))+geom_bar(stat = "identity")+facet_wrap(~variable,scales="free")+xlab("")+
    geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.2,
                  position=position_dodge(.9))+NoLegend()+
    ylab("Average expression (all cells)\n Scaled counts")+theme_cowplot()+scale_fill_few()
  ggsave("results/silkd60_120.barplot_genes.addedgenes.pdf")
  
  }
