
##############################################
        # reqiured packages #
##############################################


load.lib<-c("dplyr", "openxlsx", "ggplot2", "BiocManager", "UpSetR", "openxlsx")


install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

BiocManager::install(c("edgeR", "ComplexHeatmap", "topGO"))

bioc_pckg = c("edgeR", "ComplexHeatmap", "topGO") 

lapply(bioc_pckg, require, character.only = TRUE)
