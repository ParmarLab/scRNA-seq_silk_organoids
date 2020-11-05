
#######################################################
# Figure 6 Single cell transcriptomics mapping DA diversity in VM organoids
######################################################

#######################################################
# Libraries and read raw data
#######################################################

library(Seurat)
library(ggplot2)
library(ggthemes)
library(cowplot)


palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]

system("mkdir -p results")

data.seurat<-readRDS("raw_data/data.seurat.qc.std.SILK.3groups.RDS")

#######################################################
# A
#######################################################
Idents(data.seurat)<-data.seurat$orig.ident

data.seurat$treatment<-factor(data.seurat$treatment,levels = c("std","silk","Lam"))

data.seurat.x3 <- subset(data.seurat,idents=c("03standardorgday30","05standardorgday30","SilkDay30-09silk"),invert=T)

DimPlot(data.seurat.x3,group.by = "predicted.id",split.by = "orig.ident",ncol = 3)+scale_color_manual(values = palettes$`Tableau 10`$value[c(3,2,4,7,1,6,5,8)])
ggsave("Figure_silk/umap.by_sample.pdf")


set.seed(28)
data.seurat.x3@meta.data$cell_id <- rownames(data.seurat.x3@meta.data)
new_df <- data.seurat.x3@meta.data %>% group_by(treatment) %>% sample_n(5000) 
table(new_df$celltypes2)

data.seurat.qc.std.ss <- subset(data.seurat.x3,cells=new_df$cell_id)


DimPlot(data.seurat.qc.std.ss[],group.by = "predicted.id",split.by = "treatment",ncol = 2,pt.size = .3)+scale_color_manual(values = palettes$`Tableau 10`$value[c(5,2,4,7,1,6,8)])+NoLegend()

#######################################################
# B
#######################################################

m1<-reshape2::melt(prop.table(table(data.seurat.x3$predicted.id,factor(data.seurat.x3$orig.ident)),2))

m1$g<-c(rep("std",21),rep("silk",21), rep("silkLam",21))
m1$value<-m1$value*100



ggplot(m1,aes(x=g,y=value,fill=Var1))+geom_boxplot()+facet_wrap(~Var1,scales="free")+ylab("Percent cells in cluster")+
  scale_fill_manual(values = palettes$`Tableau 10`$value[c(3,2,4,7,1,6,5,8)])


#######################################################
# C
#######################################################

genes <- c("DCX","STMN2","PBX1","TH")

setdiff(genes,rownames(data.seurat))

Idents(data.seurat)<-data.seurat$celltypes2


for(g in genes){

    FeaturePlot(data.seurat,g,order=T,min.cutoff = "q9",label=F,split.by = "treatment")
}


#######################################################
# D
######################################################



md1 <- data.frame(data.seurat.x3@meta.data,
                  t(data.frame(GetAssay(data.seurat.x3,assay = "RNA")[c(intersect(rownames(data.seurat.x3),genes)),])))
se <- function(x) sqrt(var(x)/length(x))

# TH high/low

md1m <- reshape2::melt(md1[,c("treatment",genes)])
means<-cbind(aggregate(value~variable+treatment,md1m,mean),sem=aggregate(value~variable+treatment,md1m,se)[,"value"])

#p1<-DimPlot(data.seurat.qc.std.dop, reduction = "umap", group.by = "DA_type", pt.size = .1, split.by = 'groups')
ggplot(means,aes(x=treatment,y=value,fill=treatment))+geom_bar(stat = "identity")+facet_wrap(~variable,scales="free")+xlab("")+
  geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.2,
                position=position_dodge(.9))+NoLegend()+
  ylab("Average expression (all cells)\n Scaled counts")+scale_fill_fivethirtyeight()
