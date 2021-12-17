
##############################################
        # reqiured packages #
##############################################

install.packages(c("dplyr", "openxlsx", "ggplot2", "BiocManager", "UpSetR", "openxlsx"))

BiocManager::install(c("edgeR", "ComplexHeatmap", "topGO"))

pckg = c("dplyr", "openxlsx", "ggplot2", "BiocManager", "UpSetR", "openxlsx", "edgeR", "ComplexHeatmap", "topGO") 

lapply(pckg, require, character.only = TRUE)
