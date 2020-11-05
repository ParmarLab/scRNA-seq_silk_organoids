
# load packages
library(Seurat)
library(ggplot2)
library(ggthemes)
library(cowplot)

palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]


# create folder for output
system("mkdir -p results")


source("R/Figure2.R")
source("R/Figure3.R")
source("R/Figure4.R")
source("R/Figure6.R")