##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Petter Storm
##' @export
create_fig2 <- function() {

  
  #######################################################
  # Figure 2 Single cell transcriptomics identifying VM organoid cell types
  ######################################################
  
  
  #######################################################
  # Libraries and read raw data
  #######################################################
  
  data.seurat.qc.std <- readRDS("raw_data/data.seurat.qc.std.20201007.rds")
  data.seurat.qc.std$Cycling <- ifelse(data.seurat.qc.std$Phase=="G1","-","+")
  data.seurat.qc.std$CyclingScore <- data.seurat.qc.std$S.Score+data.seurat.qc.std$G2M.Score
  data.seurat.qc.std$Cycling2 <- ifelse(data.seurat.qc.std$S.Score+data.seurat.qc.std$G2M.Score < 0.05,"-","+")
  
  #######################################################
  # B
  #######################################################
  
  DimPlot(data.seurat.qc.std,split.by = "celltypes2")+scale_color_manual(values = palettes$`Tableau 10`$value[c(1,7,4,2,3,8,6,5)])+theme(legend.position = "top")
  ggsave("results/2B.png",w=7,h=7)
 
 
  #######################################################
  # C
  #######################################################
  
  DimPlot(data.seurat.qc.std, group.by = "groups")+scale_color_manual(values = palettes$`Classic Gray 5`$value)+theme(legend.position = "top")
  ggsave("results/2C.pdf",w=7,h=7)
  
  
  
  #######################################################
  # D
  #######################################################
  
  genelist <-read.csv("raw_data/genelist_paper.2020130.2.csv",header=T)
  
  dpd<-DotPlot(data.seurat.qc.std,features = unique(genelist[,"Gene"]))
  dpd.d<-merge(dpd$data,genelist,by.x="features.plot",by.y="Gene")
  
  dpd.d$plot_labels <- factor(substr(dpd.d$id,1,5),
                              levels = c("Neura","FP1","FP2","FP3","Dopam","Astro","Oligo","VLMC"))
  # Finally plot
  
  pal<-tableau_color_pal(palette = "Tableau 10")(10)
  
  dpd.d$features.plot2 <- factor(dpd.d$features.plot,levels=rev(levels(dpd.d$features.plot)))
  
  ggplot(dpd.d,aes(x=plot_labels,y=features.plot2,size=pct.exp,col=Cluster,alpha=avg.exp.scaled))+
    geom_point()+ scale_alpha(range = c(-1, 2))+
    xlab("") + labs(size = "Percent expressed",col="Cell type")+
    scale_color_manual(values = pal[c(3,2,1,4,5,8)])+ylab("")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_cowplot()
  ggsave("results/2D.png",w=4,h=7)
  
  
  
  #######################################################
  # E
  #######################################################
  
  
  DimPlot(data.seurat.qc.std, group.by = "Cycling")+theme(legend.position = "top")+scale_color_manual(values = c("#05afc9","#c93c05"))
  ggsave("results/2E.pdf")
  
    
  
  #######################################################
  # F
  #######################################################
  
  
  FeaturePlot(data.seurat.qc.std,features="TOP2A",order=T,min.cutoff = "q9")+NoLegend()+NoAxes()
  ggsave("results/2E.TOP2A.png")
  
  FeaturePlot(data.seurat.qc.std,features="CCNB2",order=T,min.cutoff = "q9")+NoLegend()+NoAxes()
  ggsave("results/2E.CCNB2.png")
  
  
  #######################################################
  # F+G
  #######################################################
  
  
  data.seurat.qc.std.no90 <- subset(data.seurat.qc.std,subset = groups=="day090",invert=T)
  m1<-reshape2::melt(prop.table(table(data.seurat.qc.std.no90$celltypes2,factor(data.seurat.qc.std.no90$orig.ident)),2))
  m1$g<-stringr::str_extract(m1$Var2, "day\\d{2,3}")
  m1$value<-m1$value*100
  
  
  my_color_palette <-palettes$`Tableau 10`$value[1:10]
  
  means1 <- m1[grep("FP|stem",m1$Var1),]
  means1$g <- as.integer(factor(means1$g,levels = c("day15"  ,"day30"  ,"day60", "day120"  )))
  means2 <- m1[grep("FP|stem",m1$Var1,invert = T),]
  means2$g <- as.integer(factor(means2$g,levels = c("day15"  ,"day30"  ,"day60", "day120"  )))
  
  
  g1<-ggplot(means1,aes(x=g,y=value,col=Var1,group=Var1))+scale_color_manual(values = my_color_palette[c(1,7,4,6)])+ylab("% cells")+xlab("")+ylab("% cells")+
    geom_point(stat="identity",alpha=.5)+
    ylab("")+geom_smooth(se = F)+  scale_x_continuous(labels = c("0"="day15"  ,"1"="day30"  ,"2"="day60", "3"="day120"  ))+theme(legend.position="top")+theme_cowplot()
  
  
  g2<-ggplot(means2,aes(x=g,y=value,col=Var1,group=Var1))+scale_color_manual(values = my_color_palette[c(2,3,8,5)])+xlab("")+ylab("% cells")+geom_point(stat="identity", alpha=.5)+
    ylab("")+geom_smooth(se = F)+  scale_x_continuous(labels = c("0"="day15"  ,"1"="day30"  ,"2"="day60", "3"="day120"  ))+ 
    theme(legend.position="top")+ theme_cowplot()
  
  
  g1+g2
  ggsave("results/2FG.png",w=14,h=7)
 
  
  #######################################################
  # S2A Distribution umi/mt
  #######################################################
  
  
  data.seurat.qc.std[["percent.mt"]] <- PercentageFeatureSet(data.seurat.qc.std, pattern = "^MT-")
  
  Idents(data.seurat.qc.std)<-data.seurat.qc.std$orig.ident
  
  
  StackedVlnPlot(data.seurat.qc.std,features = c("nFeature_RNA","percent.mt"),pt.size = 0)+theme(axis.text.x = element_text(angle = 90))
  
  
  ggsave("~/Dropbox/Parmar/Alessandro/ale_organoids_paper_figs/nature_comms_revision/high_res_scrnaseq/violon_qc_nfeature_percentmt.pdf",h=8)
  
  
  #######################################################
  # S4A
  #######################################################
  
  DimPlot(data.seurat.qc.std,split.by = "groups")+scale_color_manual(values = palettes$`Tableau 10`$value[c(1,7,4,2,3,8,6,5)])+theme(legend.position = "top")
  ggsave("results/S4A.pdf",w=12,h=7)
  
  
  
  #######################################################
  # S4B
  #######################################################
  
  m1<-reshape2::melt(prop.table(table(factor(data.seurat.qc.std$orig.ident),data.seurat.qc.std$celltypes2),1))
  m1$Var2 <- as.character(m1$Var2)
  m1[m1$Var2=="Neural stem cells","Var2"]<-"FP0"
  m1[m1$Var2=="Oligodendrocytes","Var2"]<-"OPC"
  m1$Var2<-factor(m1$Var2,levels=c("FP0", "FP1", "FP2", "FP3","Dopamine neurons","Astrocytes","OPC","VLMC" ))
  m1$day <- factor(stringr::str_match(m1$Var1,"day\\d{1,3}"),levels=c("day15", "day30" ,"day60", "day90","day120"))
  m1$value<-m1$value*100
  
  means<-cbind(aggregate(value~Var2+day,m1,mean),sem=aggregate(value~Var2+day,m1,se)[,"value"])
  
  
  
  ggplot(means,aes(x=day,y=value,fill=Var2))+geom_bar(stat = "identity",color="black")+facet_wrap(~Var2,scales="free",ncol = 4)+xlab("")+
    geom_errorbar(aes(ymax=value+sem,ymin=value), width=.2,
                  position=position_dodge(.9))+NoLegend()+
    ylab("Percent cells in cluster")+theme_cowplot()+
    scale_fill_manual(values = palettes$`Tableau 10`$value[c(6,4,7,1,2,3,5,8)])+NoLegend()+geom_point(data=m1,aes(x=day,y=value))
  ggsave("results/figure.s4b.pdf",w=12)
  
  
  #######################################################
  # S4C
  #######################################################
  
  DimPlot(data.seurat.qc.std,reduction = "pca")+scale_color_manual(values = palettes$`Tableau 10`$value[c(1,7,4,2,3,8,6,5)])+theme(legend.position = "top")
  ggsave("results/S4A.pdf",w=12,h=7)
  
  
  
  return(NULL)

}
