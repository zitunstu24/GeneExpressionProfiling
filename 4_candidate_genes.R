# required packages

load.lib<-c("dplyr", "openxlsx", "ggplot2", "BiocManager", "UpSetR", "openxlsx")

install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

BiocManager::install(c("edgeR", "ComplexHeatmap", "topGO"))

bioc_pckg = c("edgeR", "ComplexHeatmap", "topGO") 

lapply(bioc_pckg, require, character.only = TRUE)

######################################################################################
######################################################################################
                         # candidate genes  #
######################################################################################
######################################################################################




rm(list=ls())

set.seed(123)

master.table = read.csv("inputs/master-table.csv")

master.table = master.table %>%
  filter(species == "Solanum lycopersicum")


combined.logFC = read.csv("inputs/combined.logFC.csv")
row.names(combined.logFC) = combined.logFC$X
combined.logFC = combined.logFC[, -1]

condition_1 <- function(x) {
  ifelse(x > 1, 1, 0)
}

input_mat = data.frame(lapply(combined.logFC, condition_1))
row.names(input_mat) = row.names(combined.logFC)
index = cbind(input_mat, total = rowSums(input_mat))

index = index %>%
  mutate(gene_id = row.names(index))
candidate = filter(index, total >24)
candidates = candidate$gene_id
candidate_exp = combined.logFC[candidates,]
row.names(candidate_exp) = candidate$gene_id
colnames(candidate_exp)
names(candidate_exp) <- gsub(x = names(candidate_exp), pattern = "\\.", replacement = "-") 
gene_id = read.csv("inputs/candidate_new.csv")
row.names(candidate_exp) = gene_id$id

col = "e1-1"
tissue = unique(master.table[master.table$sample.group == col, ]$tissue)

stress = unique(master.table[master.table$sample.group == col, ]$stress_type)
duration = unique(master.table[master.table$sample.group == col, ]$stress.duration)
temp = unique(master.table[master.table$sample.group == col, ]$temperature)
genotype = unique(master.table[master.table$sample.group == col, ]$genotype.name)

design.logFC = data.frame(tissue, stress, duration, temp, genotype)

col = colnames(candidate_exp)

for (i in 2:length(col)) {
  tissue = unique(master.table[master.table$sample.group == col[i], ]$tissue)
  stress = unique(master.table[master.table$sample.group == col[i], ]$stress_type)
  duration = unique(master.table[master.table$sample.group == col[i], ]$stress.duration)
  temp = unique(master.table[master.table$sample.group == col[i], ]$temperature)
  genotype = unique(master.table[master.table$sample.group == col[i], ]$genotype.name)
  if(is.null(genotype)) {
    genotype = unique(master.table[master.table$sample.group == col[i], ]$genotype.accession)
    design.rest = data.frame(tissue, stress, duration, temp, genotype)
    design.logFC = rbind(design.logFC, design.rest)
  }
  else {
    design.rest = data.frame(tissue, stress, duration, temp, genotype)
    design.logFC = rbind(design.logFC, design.rest)
    print(col[i])
  }
}
genotype = design.logFC$genotype
genotype = as.data.frame(genotype)
stress = design.logFC$stress
stress = as.data.frame(stress)
duration = design.logFC$duration
duration = as.data.frame(duration)
temp = design.logFC$temp
temp = as.data.frame(temp)

ann_gen = HeatmapAnnotation(
  df = cbind(duration, temp),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 7))

tissue = design.logFC$tissue
tissue = as.data.frame(tissue)
ann_tissue = HeatmapAnnotation(df =  cbind(genotype, tissue),
                               annotation_name_gp = gpar(fontsize = 7))



candidate_exp_mat = as.matrix(candidate_exp)

pdf("outputs/candidate_genes_expression.pdf",
    width = 6.13,
    height = 5.60)
Heatmap(+
  candidate_exp_mat,
  name = "logFC",
  top_annotation = ann_gen,
  bottom_annotation = ann_tissue,
  column_title_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 5),
  row_names_gp = gpar(fontsize = 6)
)
dev.off()