master.table = read.csv("./inputs/master.table.csv")
master.table = master.table[1:190,]
combined.abundance = read.csv("./inputs/combined_abundance.tsv", sep = "\t")
row.names(combined.abundance) = combined.abundance$target_id

library(edgeR)
library(dplyr)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
# Example: Calculate DEGs and Fold-Changes for sample group e1-1
sample.group  = "e1-1"
sras.in.group = master.table[master.table$sample.group == sample.group, ]$sra_run_id
control.group = unique(master.table[master.table$sample.group == sample.group, ]$respective_control)
sras.in.control = master.table[master.table$sample.group == control.group, ]$sra_run_id

selected.abundances = combined.abundance[,c(paste0(sras.in.control, "_est_counts"),
                                            paste0(sras.in.group, "_est_counts"))]
edger.group = factor(c(rep(1, length(sras.in.control)),
                       rep(2, length(sras.in.group))))

y <- DGEList(counts = selected.abundances, group = edger.group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
res <- topTags(et, n = Inf)
res <- as.data.frame(res)
res <- res %>%
  mutate(gene_id = row.names(res))
res <- res[order(row.names(res)), ]
res <- filter(res, FDR < 0.01)
dfr <- res[, c(1, 3, 4, 5)]
dfr <- dfr %>%
  filter(logFC < -1 | logFC > 1)
names(dfr)[names(dfr) == "logFC"] <- paste0(sample.group, "_logFC")
names(dfr)[names(dfr) == "PValue"] <- paste0(sample.group, "_PValue")
names(dfr)[names(dfr) == "FDR"] <- paste0(sample.group, "_FDR")
# Example: Loop through all sample groups (then we can do as above):

sample.group = unique(master.table$sample.group)
for (i in 3:length(sample.group)) {
  sras.in.group = master.table[master.table$sample.group == sample.group[i], ]$sra_run_id
  
  control.group = unique(master.table[master.table$sample.group == sample.group[i], ]$respective_control)
  if (control.group == "") {
    next
  }
  else {
    sras.in.control = master.table[master.table$sample.group == control.group, ]$sra_run_id
    col.in.combined.abundance = colnames(combined.abundance)
    col.interest = c(paste0(sras.in.control, "_est_counts"),
                     paste0(sras.in.group, "_est_counts"))
    
    if (length(intersect(col.in.combined.abundance, col.interest)) == 0) {
      next
    }
    else {
      selected.abundances = combined.abundance[c(
        paste0(sras.in.control, "_est_counts"),
        paste0(sras.in.group, "_est_counts"))]
      
      #this is for those have no replication
      if (length(sras.in.group) == 1) {
        edger.group = factor(c(
          rep(1, length(sras.in.control)),
          rep(2, length(sras.in.group))))
        bcv <- 0.1
        y <- DGEList(counts = selected.abundances, group = edger.group)
        et <- exactTest(y, dispersion = bcv ^ 2)
        res <- topTags(et, n = Inf)
        res <- as.data.frame(res)
        res <- res %>%
          mutate(gene_id = row.names(res))
        res <-
          res[order(row.names(res)), ] #ordering the genes in row from 00 chromosome to higher
        res <- filter(res, FDR < 0.01)
        df <-
          res[, c(1, 3, 4, 5)]      # we need logFC, p, fdr value and gene_id column
        df <- df %>%
          filter(logFC < -1 | logFC > 1)
        names(df)[names(df) == "logFC"] <- paste0(sample.group[i], "_logFC")
        names(df)[names(df) == "PValue"] <- paste0(sample.group[i], "_PValue")
        names(df)[names(df) == "FDR"] <- paste0(sample.group[i], "_FDR")
        dfr = full_join(dfr, df, by = "gene_id") #different samples have different DEGs cbind doesnot work when row is different that is why full_join
        print(sample.group[i])
      }
      #this is for those have replication
      else {
        edger.group = factor(c(
          rep(1, length(sras.in.control)),
          rep(2, length(sras.in.group))))
        
        y <- DGEList(counts = selected.abundances, group = edger.group)
        y <- calcNormFactors(y)
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        et <- exactTest(y)
        
        res <- topTags(et, n = Inf)
        res <- as.data.frame(res)
        res <- res %>%
          mutate(gene_id = row.names(res))
        res <- res[order(row.names(res)), ] #ordering the genes in row from 00 chromosome to higher
        res <- filter(res, FDR < 0.01)
        df <-res[, c(1, 3, 4, 5)]      # we need logFC, p, fdr value and gene_id column
        df <- df %>%
          filter(logFC < -1 | logFC > 1)
        names(df)[names(df) == "logFC"] <- paste0(sample.group[i], "_logFC")
        names(df)[names(df) == "PValue"] <- paste0(sample.group[i], "_PValue")
        names(df)[names(df) == "FDR"] <- paste0(sample.group[i], "_FDR")
        dfr = full_join(dfr, df, by = "gene_id") #different samples have different DEGs cbind doesnot work when row is different that is why full_join
        print(sample.group[i])
      }
    }
  }
}
#heatmap for DEGs across all samples for dynamics
dfr[is.na(dfr)] <- 0
row.names(dfr) = dfr$gene_id
dfr = dplyr::select(dfr,-c('gene_id'))
combined.logFC = dfr[, c(grep("_logFC", colnames(dfr)))]
colnames(combined.logFC) = gsub("_logFC", "", colnames(combined.logFC))
colnames(combined.logFC)
col = "e1-1"
tissue = unique(master.table[master.table$sample.group == col, ]$tissue)
temp = unique(master.table[master.table$sample.group == col, ]$temperature)
genotype = unique(master.table[master.table$sample.group == col, ]$genotype.name)
design.logFC = data.frame(tissue, temp, genotype)

col = colnames(combined.logFC)

for (i in 2:length(col)) {
  tissue = unique(master.table[master.table$sample.group == col[i], ]$tissue)
  temp = unique(master.table[master.table$sample.group == col[i], ]$temperature)
  genotype = unique(master.table[master.table$sample.group == col[i], ]$genotype.name)
  if(genotype == "") {
    genotype = unique(master.table[master.table$sample.group == col[i], ]$genotype.accession)
    design.rest = data.frame(tissue, temp, genotype)
    design.logFC = rbind(design.logFC, design.rest)
  }
  else {
    design.rest = data.frame(tissue, temp, genotype)
    design.logFC = rbind(design.logFC, design.rest)
    print(col[i])
  }
}
genotype = design.logFC$genotype
genotype = as.data.frame(genotype)
temp = design.logFC$temp
temp = as.data.frame(temp)
ann_temp = HeatmapAnnotation(
  df = temp,
  show_annotation_name = FALSE,
  show_legend = TRUE,
  col = list(type = c("black", "red", "white"))
)
tissue = design.logFC$tissue
tissue = as.data.frame(tissue)
ann_tissue = HeatmapAnnotation(df = cbind(genotype, tissue),
                               annotation_name_gp = gpar(fontsize = 7))

combined.logFC.heatmap = as.matrix(combined.logFC)
ht = Heatmap(
  combined.logFC.heatmap,
  name = "logFC",
  top_annotation = ann_temp,
  bottom_annotation = ann_tissue,
  show_row_names = FALSE,
  column_title_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 5)
)


#heatmap for selected genes
library(readxl)
master.gene.table <- read_excel("master.table.2.xlsx", sheet = 2)
interested.genes = master.gene.table$genes
gene.logFC = combined.logFC
row.names(gene.logFC) =  gsub("[.].*", "",  row.names(gene.logFC)) #selected genes were not extended with ".1.2" then it should be cut off
selected.gene.logFC = gene.logFC[interested.genes, ]
selected.gene.logFC = na.omit(selected.gene.logFC)

master.gene.table = master.gene.table[-c(25:27, 46), ]
gene.type = master.gene.table$function.type
gene.type = as.data.frame(gene.type)
ann_gene.type = rowAnnotation(df = gene.type)

selected.gene.heatmap = as.matrix(selected.gene.logFC)
h = Heatmap(
  selected.gene.heatmap,
  name = "logFC",
  col = colorRamp2(c(-10, 0, 10), c("blue", "azure", "red")),
  top_annotation = ann_stress,
  bottom_annotation = ann_tissue,
  left_annotation = ann_gene.type,
  row_names_gp = gpar(fontsize = 8),
  row_title_gp = gpar(fontsize = 7),
  column_title_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 5)
)
pdf(h, "outputs/figures/heatmap.heat.pdf")
#upset plot for DEGs
condition_1 <- function(x) {
  ifelse(x > 1, 1, 0)
}
upset.upregulated <-
  as.data.frame(lapply(combined.logFC, condition_1))
u = upset(
  upset.upregulated,
  nsets = 97,
  mb.ratio = c(0.55, 0.45),
  sets.x.label = "Upregulated Genes",
  order.by = "freq"
)
pdf(u, "outputs/figures/upregulated.heat.pdf")
condition_2 <- function(x) {
  ifelse(x < -1, 1, 0)
}
upset.downregulated <-
  as.data.frame(lapply(combined.logFC, condition_2))
d = upset(
  upset.downregulated,
  nsets = 97,
  sets.x.label = "Downregulated Genes",
  mb.ratio = c(0.55, 0.45),
  order.by = "freq"
)
pdf(d, "outputs/figures/downregulated.heat.pdf")

