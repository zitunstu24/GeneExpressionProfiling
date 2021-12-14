
library(dplyr)
library(UpSetR)
library(ComplexHeatmap)
library(openxlsx)

######################################################################################
######################################################################################
         #  heat response of different tissues on moneymaker genotype  #
######################################################################################
######################################################################################

rm(list=ls())
set.seed(123)

combined.logFC = read.csv("inputs/combined.logFC.csv")
row.names(combined.logFC) = combined.logFC$X
combined.logFC = combined.logFC[, -1]

moneymaker = combined.logFC[, c("e5.1", "e6.1", "e6.3", "e8.1")]
colnames(moneymaker) = c("pollen", "leaf", "anther", "stem")

condition <- function(x) {
  ifelse(x == 0, 0, 1)
}
upset.moneymaker <-
  as.data.frame(lapply(moneymaker, condition))

pdf("outputs/tissue_heat_response.pdf",
    width = 6.13,
    height = 5.60)
upset(
  upset.moneymaker,
  nsets = 4,
  sets.x.label = "significant genes",
  mb.ratio = c(0.7, 0.3),
  point.size = 2, line.size = 0.75,
  matrix.color = "black", main.bar.color = "black",
  text.scale = c(1.3, 1.3, 1, 1, 2, 1.3),
  order.by = "freq"
)
dev.off()

moneymaker_mat = as.matrix(moneymaker)

pdf("outputs/tissue_heat_response.pdf",
    width = 6.13,
    height = 5.60)
Heatmap(
  moneymaker_mat, name = "logFC",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 0
)
dev.off()

