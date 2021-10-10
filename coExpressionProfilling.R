library(dplyr)
library(tidyverse)
library(WGCNA)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(genefilter)
master.table = read.csv("inputs/master.table.csv")
master.table = master.table[1:190,]
read.counts = read.delim("inputs/combined_abundance.tsv")
row.names(read.counts) = read.counts$target_id
colnames(read.counts)
read.counts = read.counts[,c(grep("_est_counts", colnames(read.counts)))]
#colnames(read.counts) = gsub("_est_counts","", colnames(read.counts))
colnames(read.counts)

sample.name = "m82-1-0.1"
sra.in.name = master.table[master.table$sample.name == sample.name, ]$sra_run_id
selected.column = as.data.frame(read.counts[,paste0(sra.in.name, "_est_counts")])
colnames(selected.column) = sample.name
row.names(selected.column) = row.names(read.counts)

sample.name = master.table$sample.name
for( i in 2:length(sample.name)) {
  sra.in.name = master.table[master.table$sample.name == sample.name[i], ]$sra_run_id
  columns = as.data.frame(read.counts[,paste0(sra.in.name, "_est_counts")])
  colnames(columns) = sample.name[i]
  row.names(columns) = row.names(read.counts)
  selected.column = cbind(selected.column,columns)
  print(i)
}
input = as.matrix(selected.column)
design.input = data.frame(samples = names(selected.column)) %>%
  mutate(
    group = gsub("[.].*","", samples) 
  )
d = DESeqDataSetFromMatrix(round(input), design.input, design = ~group)
d = DESeq(d)
vsd = getVarianceStabilizedData(d)
rv_vsd = rowVars(vsd)
q75 <- quantile( rowVars(vsd), .75)
q95 <- quantile( rowVars(vsd), .95)
normalised.d = vsd[rv_vsd>q95, ]
pivot.d = data.frame(normalised.d) %>%
  mutate(
    id = row.names(normalised.d)
  ) %>&
  pivot_longer(-id)
pivot.d %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "samples",
    y = "normalized expression"
  )
input_mat = t(normalised.d)
picked_power <- 12
temp_cor <- cor
cor <- WGCNA::cor
netwk <- blockwiseModules(
  input_mat,
  power = picked_power,
  networkType = "signed",
  deepSplit = 2,
  pamRespectsDendro = F,
  minModuleSize = 20,
  maxBlockSize = 5000,
  reassignThreshold = 0,
  mergeCutHeight = 0.4,
  saveTOMs = T,
  saveTOMFileBase = "ER",
  numericLabels = T,
  verbose = 3
)
mergeColors <- labels2colors(netwk$colors)
dendro <-  plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergeColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
save(dendro,"./Results/figures/dendro.png")
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
write.csv(module_df, "./data/Generated/gene_module_colors.csv", row.names = F)
MEs0 <- moduleEigengenes(input_mat, mergeColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order <- names(MEs0) %>%
  gsub("ME","",.)
MEs0$treatment = row.names(MEs0)
mME <- MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME","", name),
    name = factor(name, levels = module_order)
  )
mME %>% ggplot(., aes(x = treatment, y=name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Module-trait Relationship", x= "sample", y ="Module", fill = "corr")
save(mME, "./Results/figures/module_sample_relation.png")
module_of_interest <- c("pink", "green")
submod <- module_df %>%
  subset(colors %in% module_of_interest)
row.names(module_df) <- module_df$gene_id
subexpr <- df[submod$gene_id,]
submod_df <- data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )
submod_df %>% ggplot(., aes(x=name, y=value, group = gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "sample", y = "expression")
genes_of_interest <- module_df %>%
  subset(colors %in% module_of_interest)
expr_of_interest <- df[genes_of_interest$gene_id,]
TOM <- TOMsimilarityFromExpr(t(expr_of_interest),
                             power = picked_power)
row.names(TOM) <- row.names(expr_of_interest)
colnames(TOM) <- row.names(expr_of_interest)
edge_list <- data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )
write.csv(edge_list, "./data/Generated/edge_list.csv", row.names = F)
