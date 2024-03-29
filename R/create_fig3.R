##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Petter Storm
##' @export
##' 

create_fig3 <- function() {

  
  #######################################################
  # Figure 3 Single cell transcriptomics mapping DA diversity in VM organoids
  ######################################################
  
  #######################################################
  # Libraries and read raw data
  #######################################################
  
  
  data.seurat.qc.std.dop.high_5<-readRDS("raw_data/data.seurat.qc.std.dop.high_5.rds")
  
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }
  
  system("mkdir -p results")
  
  
  #######################################################
  # Agarwal - read/prepare data
  #######################################################
  
  if(file.exists("raw_data/agarwal.all.rds")){
    Agarwal<-readRDS("raw_data/agarwal.all.rds")
  }else{
    # READ 10X
    dirs <- list.dirs("~/Dropbox/Parmar/GSE140231_RAW/",full.names = T)[-c(1,2)]
    names(dirs)<-gsub("^.*_", "", dirs)
    human<-Read10X(dirs)
    
    # READ METADATA
    metad <- data.frame(readxl::read_excel("~/Dropbox/Parmar/GSE140231_RAW/metadata-all.xlsx"))
    metad$Samplenames<-paste0(metad$Samplenames,"-1")
    setdiff(metad$Samplenames,colnames(human))
    rownames(metad)<-metad$Samplenames
    
    # RUN SEURAT
    Agarwal<-CreateSeuratObject(human,meta.data = metad)
    Agarwal <- NormalizeData(Agarwal)
    Agarwal <- FindVariableFeatures(Agarwal, selection.method = "vst", nfeatures = 6000)
    all.genes <- rownames(Agarwal)
    Agarwal <- ScaleData(Agarwal, features = all.genes)
    Agarwal <- RunPCA(Agarwal, features = VariableFeatures(object = Agarwal))
    Agarwal <- FindNeighbors(Agarwal, dims = 1:10)
    Agarwal <- FindClusters(Agarwal, resolution = 0.5)
    Agarwal <- RunUMAP(Agarwal, dims = 1:10)
    Agarwal$Level_1_cell_type_region <- paste0(Agarwal$Level_1_cell_type," - ",Agarwal$Brain_Region)
    saveRDS(Agarwal, file="raw_data/agarwal.all.rds")
  }
  
  Agarwal$Level_1_cell_type_region <- paste0(Agarwal$Level_1_cell_type," - ",Agarwal$Brain_Region)
  
  #######################################################
  # A Agarwal on organiods
  #######################################################
  
  anchors <- FindTransferAnchors(reference = Agarwal, query = data.seurat.qc.std, 
                                 dims = 1:30)
  
  predictions <- TransferData(anchorset = anchors, refdata = Agarwal$Level_1_cell_type_region, k.weight = 20,
                              dims = 1:30)
  
  
  data.seurat.qc.std <- AddMetaData(data.seurat.qc.std, metadata = predictions)
  
  DimPlot(data.seurat.qc.std,group.by = "predicted.id")
  
  
  
  
  #######################################################
  # B
  #######################################################
  
  
  t1<-table(data.seurat.qc.std$predicted.id , data.seurat.qc.std$celltypes2)
  pt1m<-melt(prop.table(t1,2))
  
  colnames(pt1m)<-c("linnar","org","value")
  ggplot(pt1m,aes(x=org,y=linnar,size=value,color=org))+geom_point()+scale_radius(range = c(1,40))+theme_cowplot()+
    scale_color_manual(values = palettes$`Tableau 10`$value[c(1,7,4,2,3,8,6,5)])+xlab("")+ylab("")+NoLegend()+ggtitle("")
  
  
  #######################################################
  # D
  #######################################################
  
  
  m1<-reshape2::melt(prop.table(table(data.seurat.qc.std.dop.high_5$EaLa,data.seurat.qc.std.dop.high_5$orig.ident),2))
  m2<-reshape2::melt(prop.table(table(data.seurat.qc.std.dop.high_5$EaLa,data.seurat.qc.std.dop.high_5$groups),2))
  
  m1$Var1<-factor(m1$Var1)
  m1 <- m1[complete.cases(m1),]
  m1$D <- paste0("Day",gsub("\\d*standardorgday","",as.character(m1$Var2),))
  
  m1$D<-factor(m1$D, levels= c("Day15","Day30", "Day60", "Day90","Day120"))
  m1$value<-m1$value*100
  tgc <- summarySE(m1, measurevar="value", groupvars=c("Var1","D"))
  
  tgcm1<-cbind(se=tgc[order(tgc$D,tgc$Var1),][,c("se")],m2)
  tgcm1$value <- tgcm1$value*100
  
  ggplot(tgcm1,aes(x=Var2,y=value,col=Var1,group=Var1))+theme(legend.position = "none")+
    geom_point()+geom_line(size=1.2)+scale_color_manual(values = c("darkgrey","#ce9f08"))+
    ylab("% of cells")+xlab("")+ geom_errorbar(aes(ymax=value+se,ymin=value-se), width=.1)+theme_cowplot()
  
  ggsave("results/3D.png",w=7,h=7)
  
  
  
  #######################################################
  # E
  #######################################################
  
  genes <- c("ID3","NES","SOX2","RFX4","LMX1A","NR4A2","SCN2A","SLC18A2","KCNC2","OTX2","TH","PITX3")
  
  md1 <- data.frame(data.seurat.qc.std.dop.high_5@meta.data,
                    t(data.frame(GetAssay(data.seurat.qc.std.dop.high_5,slot = "count")[c(intersect(rownames(data.seurat.qc.std.dop.high_5),genes),"TH"),])))
  se <- function(x) sqrt(var(x)/length(x))
  
  
  # TH high/low
  md1m <- reshape2::melt(md1[,c("DA_type","ID3","NES","SOX2","RFX4","LMX1A","NR4A2","SCN2A","SLC18A2","KCNC2","OTX2","TH","PITX3")])
  md1m$DA_type2 <- ifelse(md1m$DA_type=="TH Low","DA early","DA late")
  
  means<-cbind(aggregate(value~variable+DA_type2,md1m,mean),sem=aggregate(value~variable+DA_type2,md1m,se)[,"value"])
  
  ggplot(means,aes(x=DA_type2,y=value,fill=DA_type2))+geom_bar(stat = "identity")+facet_wrap(~variable,scales="free")+
    geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.2,
                  position=position_dodge(.9))+NoLegend()+ylab("")+theme_cowplot()+scale_fill_manual(values = c("darkgrey","#ce9f08"))
  
  ggsave("results/3E.png",w=7,h=7)
  
  #######################################################
  # F
  #######################################################
  
  pal <- few_pal()(5)
  
  DimPlot(data.seurat.qc.std.dop.high_5,group.by = "DA_type5")+scale_color_manual(values = pal)
  ggsave("results/3F.pdf",w=7,h=7)
  
  #######################################################
  # J
  #######################################################
  
  cult.reclust_dopa<-readRDS("../data.seurat.fetal.rds")
  DimPlot(cult.reclust_dopa,,pt.size = 1.02)+scale_color_manual(values = c("#ababab","#595959","#9a532b","#503143","#a2c8ec"))
  DimPlot(cult.reclust_dopa,,pt.size = 1.02)+scale_color_manual(values = c("#ababab","#595959","#f77f00","#d62828","#"))
  ggsave("results/3I.fetalumap.pdf",w=12,h=8)
  

  
  #######################################################
  # H - DopamineHigh - heatmap
  #######################################################
  
  col.pal <- RColorBrewer::brewer.pal(9,"GnBu")
  aveE <- AverageExpression(data.seurat.qc.std.dop.high_5)$RNA
  aveE <- aveE[,order(colnames(aveE),decreasing = T)]
  
  genes<-read.csv("raw_data/genelist_da_subcluster.20200327.csv",header=F)
  genes$General<-toupper(genes$V1)
  genes[genes$V1=="FOXP1","V2"]<-"DA-2"
  genes <- genes[order(genes$V2),]
  genes.toplot <- intersect(genes$General,rownames(aveE))
  setdiff(genes$General,rownames(aveE))
  
  genes<-genes[!duplicated(genes$General),]
  rownames(genes)<-genes$General
  genes$D1 <- 1
  
  col.pal <- RColorBrewer::brewer.pal(9,"GnBu")
  
  colfunc <- colorRampPalette(c("#BE6216","#DCD2D0" ,"#796AA3"))
  
  png("results/3h.png",w=1000,h=400)
  pheatmap(t(aveE[genes.toplot,]),scale="column", color=colfunc(100),
           cluster_rows = F,cluster_cols = F,gaps_row = c(1,2,3),gaps_col = c(12,12+28),
           border_color = "lightgrey")
  dev.off()
  
  
  
  
  return(NULL)

}
