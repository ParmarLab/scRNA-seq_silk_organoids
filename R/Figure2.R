
#######################################################
# Figure 2 Single cell transcriptomics identifying VM organoid cell types
######################################################

#a, Uniform manifold approximation and projection (UMAP) plot showing clustering of 91,034 analyzed cells from VM organoids at day 15, 30, 60, 90, and 120. Cell type assignments are indicated.
#b, UMAP plot of cells color-coded by organoid of origin. c, Dot plot showing expression levels of indicated genes per cluster. Indicated genes are established markers for neural progenitors, 
#floor plate progenitors, DA neurons, astrocytes, oligodendrocyte and vascular leptomeningeal cells. d, Silhouette plot and scores showing the difference between mean distance to cells within the same and nearest clusters. 
#e, f, Representative curves showing proportion of cells as a percentage for each cluster during VM organoid differentiation (day 15–120). g, RNA velocity plotted on UMAP representing location of estimated future cell state.
#h, Developmental trajectory from pluripotency to terminally differentiated stages reconstructed using SPRING in VM organoid. Pseudocells are color-graded by total count. Bottom right, pseudocells are color-coded by cell type assignments.  i, SPRING plot colored (purple) by marker gene expression of emerging cellular clusters. 



#######################################################
# Libraries and read raw data
#######################################################

library(Seurat)
library(ggplot2)
library(ggthemes)
library(cowplot)

data.seurat.qc.std <- readRDS("raw_data/data.seurat.qc.std.20201007.rds")
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]


data.seurat.qc.std$Cycling <- ifelse(data.seurat.qc.std$Phase=="G1","-","+")

system("mkdir -p results")

#######################################################
# A
#######################################################

DimPlot(data.seurat.qc.std, group.by = "celltypes2")+scale_color_manual(values = palettes$`Tableau 10`$value[c(1,7,4,2,3,8,6,5)])+theme(legend.position = "top")
ggsave("results/2A.png",w=7,h=7)


#######################################################
# B
#######################################################

DimPlot(data.seurat.qc.std, group.by = "groups")+scale_color_manual(values = palettes$`Classic Gray 5`$value)+theme(legend.position = "top")
ggsave("results/2B.png",w=7,h=7)

#######################################################
# C
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
ggsave("results/2C.png",w=4,h=7)


#######################################################
# D
#######################################################


DimPlot(data.seurat.qc.std, group.by = "Cycling")+theme(legend.position = "top")
ggsave("results/2D.png",w=7,h=7)

#######################################################
# E
#######################################################


FeaturePlot(data.seurat.qc.std,features="TOP2A",order=T,min.cutoff = "q9")+NoLegend()+NoAxes()
ggsave("results/2E.TOP2A.png")

FeaturePlot(data.seurat.qc.std,features="CCNB2",order=T,min.cutoff = "q9")+NoLegend()+NoAxes()
ggsave("results/2E.CCNB2.png")


#######################################################
# F
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
