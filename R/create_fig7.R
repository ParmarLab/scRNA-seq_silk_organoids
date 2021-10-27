
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Petter Storm
##' @export
create_fig7 <- function() {
  
  
  #######################################################
  # Figure 7 Stress response
  ######################################################
  

  suppressPackageStartupMessages(library(escape))
  data.seurat<-readRDS("raw_data/data.seurat.qc.std.SILK.3groups.RDS")
  Idents(data.seurat)<-data.seurat$orig.ident

  data.seurat$treatment<-factor(data.seurat$treatment,levels = c("std","silk","Lam"))

  data.seurat.x3 <- subset(data.seurat,idents=c("03standardorgday30","05standardorgday30","SilkDay30-09silk"),invert=T)
  data.seurat.x3$groups<-ifelse(data.seurat.x3$treatment=="std","std","lam")
  
  data.seurat.x3$CellID <- rownames(data.seurat.x3@meta.data)
  df <- data.seurat.x3@meta.data %>% group_by(treatment) %>% sample_n(size = 2500)
  data.seurat.x3.ss <- data.seurat.x3[,df$CellID]
  data.seurat.x3.ss<-subset(data.seurat.x3.ss,subset = groups != "silk")
  data.seurat.x3.ss$groups<-factor(data.seurat.x3.ss$groups)
data.seurat.x3$groups<-data.seurat.x3$treatment

GS.go <- getGeneSets(library = "C5",gene.sets = c("GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"  ))

GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING<-read.delim("http://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING&fileType=txt")



# day 60/120
std_silk.merge<-readRDS("raw_data/std_silk.merge.rds")
df2 <- data.frame(std_silk.merge@meta.data %>% group_by(treat_group) %>% sample_n(size = 5000))
std_silk.merge.ss <- std_silk.merge[,df2[df2$groups=="day120","CellID"]]
std_silk.merge.ss$groups<-factor(std_silk.merge.ss$treat_group)
DimPlot(std_silk.merge.ss,group.by = "groups")

ES.seurat.go <- enrichIt(obj = std_silk.merge.ss, gene.sets = GS.go, groups = 1000, cores = 4)

data.seurat.x3.ss.md <- cbind(data.seurat.x3.ss@meta.data,ES.seurat.go)
bp1<-ggplot(data.seurat.x3.ss.md,aes(x=groups,y=GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING,fill=groups))+
  geom_boxplot()+theme_cowplot()+scale_fill_manual(values = c("#65997C","#9FA1A0"))+ coord_flip()+NoLegend()+xlab("")+ggtitle("All cells")

bp1

#######################################################
# Stress comparison
#######################################################


#ALL CELLS
DimPlot(data.seurat.x3,group.by = "predicted.id",split.by = "orig.ident",ncol = 3)+scale_color_manual(values = palettes$`Tableau 10`$value[c(3,2,4,7,1,6,5,8)])

genelist.stress <- c("YIPF5","ARCN1","GORASP2","PGK1","ALDOA","BNIP3")
data.seurat.x3[["stress"]] <- PercentageFeatureSet(data.seurat.x3, features = genelist.stress)

VlnPlot(data.seurat.x3, features = c("stress"),pt.size = 0)
source("~/Dropbox/scripts/single_cell_functions.R")
dpdplot(data.seurat.x3,genelist = genelist.stress)



#Extended Data Fig. 12: Glycolysis and ER stress across culture systems.
# https://www.nature.com/articles/s41586-020-1962-0/figures/17

genelist.stress2 <- c("PGK1","ARCN1","GORASP2")

VlnPlot(data.seurat.x3,features = genelist.stress2,split.by = "treatment")
FeaturePlot(data.seurat.x3,features = genelist.stress2,split.by = "treatment",min.cutoff = "q9",order=T)


#######################################################
# Stress comparison GSEA
#######################################################


stress.gs <- list(stress=genelist.stress)

DimPlot(data.seurat.x3.ss,group.by = "treatment")


ES.seurat <- enrichIt(obj = data.seurat.x3.ss, gene.sets = GS.c2, groups = 1000, cores = 4)
data.seurat.x3.ss <- Seurat::AddMetaData(data.seurat.x3, ES.seurat)
multi_dittoPlot(data.seurat.x3.ss, vars = c("BIOCARTA_STRESS_PATHWAY"), 
                group.by = "groups", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

ES.seurat.go <- enrichIt(obj = data.seurat.x3.ss, gene.sets = GS.go, groups = 1000, cores = 4)
data.seurat.x3.ss <- Seurat::AddMetaData(data.seurat.x3, ES.seurat.go)
multi_dittoPlot(data.seurat.x3.ss, vars = c("GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"  , "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS",
                                            "GOBP_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS","GOBP_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS",
                                            "GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS"  ), 
                group.by = "groups", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))




ES.seurat <- enrichIt(obj = std_silk.merge.ss, gene.sets = GS.c2, groups = 1000, cores = 4)
std_silk.merge.ss <- Seurat::AddMetaData(std_silk.merge.ss, ES.seurat)
multi_dittoPlot(std_silk.merge.ss, vars = c("BIOCARTA_STRESS_PATHWAY"), 
                group.by = "groups", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

ES.seurat.go <- enrichIt(obj = std_silk.merge.ss, gene.sets = GS.go, groups = 1000, cores = 4)
std_silk.merge.ss <- Seurat::AddMetaData(std_silk.merge.ss, ES.seurat.go)
multi_dittoPlot(std_silk.merge.ss, vars = c("GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"  , "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS",
                                            "GOBP_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS","GOBP_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS",
                                            "GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS"  ), 
                group.by = "groups", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))


std_silk.merge.ss$treat_group_relevel<-factor(std_silk.merge.ss$treat_group,levels = c("std day120","std day060", "silk day060", "silk day120"))
summary(lm(GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING~treat_group_relevel,std_silk.merge.ss@meta.data))

std_silk.merge.ss.no_day60 <- std_silk.merge.ss[,grep("120",std_silk.merge.ss$treat_group)]

multi_dittoPlot(std_silk.merge.ss.no_day60, vars = c("GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"  ), 
                group.by = "groups", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", ncol = 1,nrow = 1,
                theme = theme_cowplot() + theme(plot.title = element_text(size = 14)))
ggsave("results/violin.stress_response.pdf",w=8,h=10)


ES2 <- data.frame(std_silk.merge.ss.no_day60[[]], Idents(std_silk.merge.ss.no_day60))
ES2$treat_group_cellt <- paste(ES2$treat_group,ES2$celltypes2)

ridgeEnrichment(ES2, gene.set = "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING", group = "celltypes2", add.rug = TRUE,facet = "treat_group")




ES.seurat <- enrichIt(obj = std_silk.merge, gene.sets = GS.go, groups = 1000, cores = 4)
std_silk.merge <- Seurat::AddMetaData(std_silk.merge, ES.seurat)
std_silk.merge.no_day60 <- std_silk.merge[,grep("120",std_silk.merge$treat_group)]
ES2 <- data.frame(std_silk.merge.no_day60[[]], Idents(std_silk.merge.no_day60))


ES3 <- ES2[grep("120",ES2$treat_group),]
ridgeEnrichment(ES3[ES3$celltypes2!="Neural stem cells",], gene.set = "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING", group = "celltypes2", add.rug = TRUE,facet = "treat_group")
summary(lm(GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING~treat*celltypes2,ES3[ES3$celltypes2!="Neural stem cells",]))


ES3$treat_group2<-ifelse(ES3$treat_group=="std day120","VM Org","Silk-lam VM Org")
r1<-ridgeEnrichment(ES3[ES3$celltypes2!="Neural stem cells",], gene.set = "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING", group = "treat_group2", add.rug = TRUE)+theme_cowplot()+ylab("")+scale_fill_manual(values = c("#65997C","#9FA1A0"))+ggtitle("All cells")
ggsave("results/ridge.integrated_stress.pdf",w=10,h=6)


r2<-ridgeEnrichment(ES3[ES3$celltypes2=="Dopamine neurons",], gene.set = "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING", group = "treat_group2", add.rug = TRUE)+theme_cowplot()+ylab("")+scale_fill_manual(values = c("#65997C","#9FA1A0"))+ggtitle("Dopamine neurons")
r1+r2

bp1<-ggplot(ES3,aes(x=treat_group2,y=GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING,fill=treat_group2))+geom_boxplot()+theme_cowplot()+scale_fill_manual(values = c("#65997C","#9FA1A0"))+ coord_flip()+NoLegend()+xlab("")+ggtitle("All cells")
bp2<-ggplot(ES3[ES3$celltypes2=="Dopamine neurons",],aes(x=treat_group2,y=GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING,fill=treat_group2))+geom_boxplot()+theme_cowplot()+scale_fill_manual(values = c("#65997C","#9FA1A0"))+ coord_flip()+NoLegend()+xlab("")+ggtitle("")

plot_grid(r1,r2,bp1,bp2,nrow = 2,labels = "AUTO")
ggsave("results/INTEGRATED_STRESS_RESPONSE_x4.pdf",w=16,h=10)



#INDIVID GENES
ale.final.old <- c("GSTP1","SOD1","ENO1","PGK1","PDIA3","PARK7")
ale.final<-c("PGK1", "PDIA3","ENO1","HERPUD1","PPP1R15A","DDIT3")
std_silk.merge.da <- std_silk.merge.no_day60[,std_silk.merge.no_day60$celltypes2=="Dopamine neurons"]
Idents(std_silk.merge.da)<-std_silk.merge.da$treat_group
vp1<-StackedVlnPlot(std_silk.merge.da,features =  ale.final[1:3],pt.size = 0,colors = c("#9FA1A0","#65997C"))
ggsave("results.vp.pdf",h=10)
vp2<-StackedVlnPlot(std_silk.merge.da,features =  ale.final[4:6],pt.size = 0,colors = c("#9FA1A0","#65997C"))
ggsave("results.vp2.pdf",h=10)

plot_grid(bp2,vp1,vp2,ncol = 3)
ggsave("results/boxplot.violinplot.stress_genes.da_neurons.pdf",w=15,h=6)

de1 <- FindMarkers(std_silk.merge.da,ident.1 = "std day120",ident.2 = "silk day120",features = intersect(GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING$GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING,rownames(std_silk.merge.ss)),logfc.threshold = 0)
de1$genes <- rownames(de1)

openxlsx::write.xlsx(de1,file="results/de_genes.da_neurons.integrated_stress.xlsx",overwrite = T)
plot_grid(StackedVlnPlot(std_silk.merge.da,features = rownames(de1)[1:10]),StackedVlnPlot(std_silk.merge.da,features = rownames(de1)[11:19]))
ggsave("results/de_genes.da_neurons.vln.pdf",w=8,h=12)


# LOOK AT INDIVIDUAL GENES

multi_dittoPlot(std_silk.merge.ss,vars = genelist.stress2,group.by = "treat_group",ncol = 3)+theme_cowplot()
ggsave("results/violin.stress_genes.pdf",w=8,h=6)



multi_dittoPlot(std_silk.merge.ss,var = intersect(GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING$GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING,rownames(std_silk.merge.ss)),group.by = "treat_group")
StackedVlnPlot(std_silk.merge.ss,features = intersect(GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING$GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING,rownames(std_silk.merge.ss)),group.by = "treat_group")


Idents(std_silk.merge.ss)<-std_silk.merge.ss$treat_group
Idents(std_silk.merge)<-std_silk.merge$treat_group

de1 <- FindMarkers(std_silk.merge.ss,ident.1 = "std day120",ident.2 = "silk day120",features = intersect(GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING$GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING,rownames(std_silk.merge.ss)))

DotPlot(std_silk.merge.ss,genelist = rownames(de1))
StackedVlnPlot(std_silk.merge,features =  rownames(de1))

GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS<-read.delim("https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS&fileType=txt")
de2 <- FindMarkers(std_silk.merge.ss,ident.1 = "std day120",ident.2 = "silk day120",features = intersect(GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS$GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS,rownames(std_silk.merge.ss)))
StackedVlnPlot(std_silk.merge.ss.no_day60,features =  rownames(de2)[1:10])
ggsave("top_hits_GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS.pdf",w=5,h=12)
de2$gene <- rownames(de2)
openxlsx::write.xlsx(de2,file="top_hits_GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS.xlsx",asTable = T,overwrite = T)


GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS<-read.delim("https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS&fileType=txt")
de3 <- FindMarkers(std_silk.merge.ss,ident.1 = "std day120",ident.2 = "silk day120",features = intersect(GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS$GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS,rownames(std_silk.merge.ss)))
StackedVlnPlot(std_silk.merge.ss.no_day60,features =  rownames(de3)[1:10])
ggsave("top_hits_GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS.pdf",w=5,h=12)
de3$gene <- rownames(de3)
openxlsx::write.xlsx(de3,file="top_hits_GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS.pdf.xlsx",asTable = T,overwrite = T)


ale.genes <- read.csv("~/Documents/ale_stress_genes.csv")
setdiff(ale.genes$ATF4,rownames(std_silk.merge.ss.no_day60))
plot_grid(StackedVlnPlot(std_silk.merge.ss.no_day60,features =  ale.genes$ATF4,pt.size = 1),StackedVlnPlot(std_silk.merge.ss.no_day60,features =  ale.genes$ATF4,pt.size = 0),ncol=2,labels = list("Dots","No dots"))

ggsave("ale_stress_genes.pdf",w=8,h=14)

ale.final <- c("GSTP1","SOD1","ENO1","PGK1","PDIA3","PARK7")
StackedVlnPlot(std_silk.merge.ss.no_day60,features =  ale.final[1:3],pt.size = 0)
ggsave("results.vp.pdf",h=10)
StackedVlnPlot(std_silk.merge.ss.no_day60,features =  ale.final[4:6],pt.size = 0,)
ggsave("results.vp2.pdf",h=10)




de.sel <- FindMarkers(std_silk.merge.ss,ident.1 = "std day120",ident.2 = "silk day120",features = ale.final)

#######################################################
# slingshot
#######################################################

library(slingshot)
library(RColorBrewer)


rd <- Embeddings(data.seurat.x3.ss,"umap")[,1:2]
cl <- factor(data.seurat.x3.ss$celltypes2)
condition <-factor(data.seurat.x3.ss$celltypes2)

pto <- slingshot(data=rd, clusterLabels = cl,start.clus="Neural stem cells",end.clus=c("Astrocytes","VLMC","Dopamine neurons"))

pt <- slingPseudotime(pto)
nms <- colnames(pt)

ptdf <- data.frame(pt,data.seurat.x3.ss@meta.data)
ptdf.m <- reshape2::melt(ptdf[,c("curve1","curve2","curve3", "curve4", "celltypes2","groups")])

ptdf.m$celltypes2<-factor(ptdf.m$celltypes2, levels=c("FP3", "FP2", "FP1", "Dopamine neurons", "Astrocytes", "VLMC", "Neural stem cells"))
ptdf.m$celltypes3 <- ptdf.m$celltypes2

ptdf.m[grep("FP",ptdf.m$celltypes2),"celltypes3"]<-"FP1"
#ggplot(ptdf,aes(x=curve2,y=curve3,col=seurat_clusters)) +geom_quasirandom()+scale_color_manual(values = thePalette)
ggplot(ptdf.m,aes(x=variable,y=value,col=celltypes2)) +geom_quasirandom(dodge.width = .8)+facet_wrap(~groups,nrow = 2)+
  scale_color_manual(values = palettes$`Tableau 10`$value[c(1,7,4,2,3,8,6,5)])+coord_flip()+theme_cowplot()+ylab("Pseudotime")+xlab("")
ggsave(file="results/beeswarm.pseudotime.std_silk.pdf",w=10,h=6)

ggplot(ptdf.m,aes(x=variable,y=value,col=celltypes3)) +geom_quasirandom(dodge.width = .8)+facet_wrap(~groups,nrow = 2)+
  scale_color_manual(values = palettes$`Tableau 10`$value[c(4,2,3,8,6,5)])+coord_flip()+theme_cowplot()+ylab("Pseudotime")+xlab("")
ggsave(file="results/beeswarm.pseudotime.std_silk.v2.pdf",w=10,h=6)







