library(edgeR)
library(dplyr)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(proxy)
library(openxlsx)

rm(list=ls())
set.seed(123)

master.table = read.csv("inputs/master-table.csv")
master.table = master.table %>%
  filter(species == "Solanum lycopersicum")

combined.abundance = read.csv("inputs/combined_abundance.tsv", sep = "\t")
row.names(combined.abundance) = combined.abundance$target_id

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
    
    if (length(intersect(col.in.combined.abundance, col.interest)) == 0) {  #if some sra id are missing in combined table
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

dfr[is.na(dfr)] <- 0
row.names(dfr) = dfr$gene_id
dfr = dplyr::select(dfr,-c('gene_id'))
combined.logFC = dfr[, c(grep("_logFC", colnames(dfr)))]  ##fold change values were extracted
colnames(combined.logFC) = gsub("_logFC", "", colnames(combined.logFC)) 
colnames(combined.logFC)
write.csv(combined.logFC, "inputs/combined.logFC.csv", row.names = T)

#### design of different factors e.g genotypes, tissues, duration and heat levels for annotation
col = "e1-1"
tissue = unique(master.table[master.table$sample.group == col, ]$tissue)
stress = unique(master.table[master.table$sample.group == col, ]$stress_type)
duration = unique(master.table[master.table$sample.group == col, ]$stress.duration)
temp = unique(master.table[master.table$sample.group == col, ]$temperature)
genotype = unique(master.table[master.table$sample.group == col, ]$genotype.name)

design.logFC = data.frame(tissue, stress, duration, temp, genotype)

col = colnames(combined.logFC)

for (i in 2:length(col)) {
  tissue = unique(master.table[master.table$sample.group == col[i], ]$tissue)
  stress = unique(master.table[master.table$sample.group == col[i], ]$stress_type)
  duration = unique(master.table[master.table$sample.group == col[i], ]$stress.duration)
  temp = unique(master.table[master.table$sample.group == col[i], ]$temperature)
  genotype = unique(master.table[master.table$sample.group == col[i], ]$genotype.name)
  if(genotype == "") {
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
tissue = design.logFC$tissue
tissue = as.data.frame(tissue)
duration = design.logFC$duration
duration = as.data.frame(duration)
temp = design.logFC$temp
temp = as.data.frame(temp)

ann_gen = HeatmapAnnotation(
  df = cbind(duration, temp),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 7))


ann_tissue = HeatmapAnnotation(df =  cbind(genotype, tissue),
                               annotation_name_gp = gpar(fontsize = 7))

DEGs_exp = as.matrix(combined.logFC)
pdf("outputs/DEGs.pdf",
    width = 6.13,
    height = 5.60)
Heatmap(
  DEGs_exp,
  name = "logFC",
  top_annotation = ann_gen,
  show_row_names = FALSE,
  bottom_annotation = ann_tissue,
  column_title_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 5)
)
dev.off()

