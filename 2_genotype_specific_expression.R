
# required packages

load.lib<-c("dplyr", "openxlsx", "ggplot2", "BiocManager", "UpSetR", "openxlsx", "Cairo")

install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

BiocManager::install(c("edgeR", "ComplexHeatmap", "topGO"))

bioc_pckg = c("edgeR", "ComplexHeatmap", "topGO") 

lapply(bioc_pckg, require, character.only = TRUE)

rm(list=ls())
set.seed(123)

combined.logFC = read.csv("inputs/combined.logFC.csv")
row.names(combined.logFC) = combined.logFC$X
combined.logFC = combined.logFC[, -1]


######################################################################################
######################################################################################
            #  heat response of different genotypes on pollen  #
######################################################################################
######################################################################################


pollen = combined.logFC[, c("e3.1", "e5.1", "e5.2", "e9.1")]
colnames(pollen) = c("Ha.3017","Moneymaker", "Red Setter","Ha.3042")

condition <- function(x) {
  ifelse(x == 0, 0, 1)
}
upset.pollen <-
  as.data.frame(lapply(pollen, condition))
pdf("outputs/genotypes_response_in_pollen.pdf",
    width = 6.13,
    height = 5.60)
upset(
  upset.pollen,
  nsets = 4,
  sets.x.label = "significant genes",
  mb.ratio = c(0.7, 0.3),
  point.size = 2, line.size = 0.75,
  matrix.color = "black", main.bar.color = "black",
  text.scale = c(1.3, 1.3, 1, 1, 2, 1.3),
  order.by = "freq"
)
dev.off()

pollen_mat = as.matrix(pollen)

pdf("outputs/genotypes_response_in_pollen_heatmap.pdf",
    width = 6.13,
    height = 5.60)
Heatmap(
  pollen_mat, name = "logFC",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 0
)
dev.off()


######################################################################################
######################################################################################
      #  heat response of different genotypes on leaf  #
######################################################################################
######################################################################################

colnames(combined.logFC)
leaf = combined.logFC[, c("e1.1", "e2.2", "e6.1", "e10.2")]
colnames(leaf) = c("M82","Jinlingmeiyu", "Moneymaker","LA1698")

condition <- function(x) {
  ifelse(x == 0, 0, 1)
}
upset.leaf <-
  as.data.frame(lapply(leaf, condition))

pdf("outputs/genotypes_response_in_leaf.pdf",
    width = 6.13,
    height = 5.60)
upset(
  upset.leaf,
  nsets = 4,
  sets.x.label = "significant genes",
  mb.ratio = c(0.7, 0.3),
  point.size = 2, line.size = 0.75,
  matrix.color = "black", main.bar.color = "black",
  text.scale = c(1.3, 1.3, 1, 1, 2, 1.3),
  order.by = "freq"
)
dev.off()

leaf_mat = as.matrix(leaf)

pdf("outputs/genotypes_response_in_leaf_heatmap.pdf",
    width = 6.13,
    height = 5.60)
Heatmap(
  leaf_mat, name = "logFC",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 0
)
dev.off()

