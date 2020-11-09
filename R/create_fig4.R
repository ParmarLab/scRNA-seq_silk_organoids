##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Petter Storm
##' @export
create_fig4 <- function() {

  
  #######################################################
  # Figure 4 Single cell transcriptomics mapping DA diversity in VM organoids
  ######################################################
  
  #######################################################
  # Libraries and read raw data
  #######################################################
  
  data.seurat.qc.std <- readRDS("raw_data/data.seurat.qc.std.20201007.rds")
  
  Idents(data.seurat.qc.std)<-data.seurat.qc.std$orig.ident
  av.exp <- AverageExpression(data.seurat.qc.std,assays = "RNA")$RNA
  
  
  # remove stanrd 5, it is silk
  av.exp<- av.exp[,!(colnames(av.exp) %in% "5standardorgday30")]
  
  
  cor.exp <- as.data.frame(cor(av.exp,method = "spearman"))
  
  
  annotation <- data.frame(group=tolower(gsub("^\\d+", "", colnames(cor.exp))))
  rownames(annotation)<-colnames(av.exp)
  annotation$day <- "d60"
  annotation[grep("30",rownames(annotation)),"day"] <- "d30"
  annotation[grep("day0",rownames(annotation)),"day"] <- "d0"
  annotation[grep("day15",rownames(annotation)),"day"] <- "d15"
  
  
  #######################################################
  # A
  #######################################################
  
  png("results/4A.png")
  pheatmap(cor.exp,annotation_col = data.frame(group=annotation[,-2],row.names = rownames(annotation)),color = redblue(200))
  dev.off()
  
  #######################################################
  # B
  #######################################################
  
  t1 <- table(data.seurat.qc.std$orig.ident,data.seurat.qc.std$celltypes)
  t2 <- table(data.seurat.qc.std$orig.ident,data.seurat.qc.std$celltypes)
  
  
  t1[grep("1standardorgday30",rownames(t1)),]
  t2[grep("1standardorgday30",rownames(t2)),]
  
  tmp <- subset(data.seurat.qc.std,idents=grep("day30",levels(Idents(data.seurat.qc.std)),value=T))
  DimPlot(tmp,group.by = "celltypes",split.by = "orig.ident",ncol = 2,label=T)
  #write.table(x = t1,file="Figure_variability/cells_per_celltyp.xls", quote=FALSE, sep='\t', col.names = NA)
  
  
  t1 <-t1[rowSums(t1)>0,]
  
  pt1<-melt(prop.table(t1,1))
  pt1$group <- gsub("\\d*standardorg","",pt1$Var1)
  pt1 <- pt1[pt1$group!="day15",]
  pt1 <- pt1[pt1$group!="day90",]
  
  # remove stanrd 5, it is silk
  pt1 <- pt1[pt1$Var1!="5standardorgday30",]
  
  pt1$group<-factor(pt1$group,levels = c("day30", "day60","day120"))
  
  pt1$Var2<-factor(pt1$Var2, levels=c("FP3", "FP2" ,"FP1" ,"Dopamine neurons" ,"Astrocytes" ,"VLMC", "Neural stem cells" ,"Oligodendrocytes"))
  
  ggplot(pt1,aes(x=Var1,y=value, fill=Var2))+geom_bar(stat="identity")+
    facet_wrap(~group,scale="free")+scale_fill_manual(values = palettes$`Tableau 10`$value[c(4,1,7,2,3,8,6,5)])+ theme(legend.position="bottom")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ theme(legend.title = element_blank())
  
  ggsave("results/4C.pdf")
  
  #######################################################
  # C
  #######################################################
  
  se <- SummarizedExperiment(as.matrix(log2(av.exp)[,grep("day15|day90", colnames(av.exp),invert = T)]),
                             colData=annotation[grep("day15|day90", annotation$group,invert = T),])
  # the call to DESeqTransform() is needed to
  # trigger our plotPCA method.
  pdf("results//4A.1.pdf")
  plotPCA( DESeqTransform( se ),intgroup="group" )
  dev.off()
  

}
