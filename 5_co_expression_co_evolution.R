
# required packages

load.lib<-c("dplyr", "openxlsx", "ggplot2", "BiocManager", "UpSetR", "openxlsx", "corrplot")

install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

BiocManager::install(c("edgeR", "ComplexHeatmap", "topGO"))

bioc_pckg = c("edgeR", "ComplexHeatmap", "topGO") 

lapply(bioc_pckg, require, character.only = TRUE)

######################################################################################
######################################################################################
                 # co expression profile of candidate genes  #
                 # co evolution profile of candidate genes   #
                 #         mining two aspects                #
######################################################################################
######################################################################################

rm(list=ls())

set.seed(123)

combined.logFC = read.csv("inputs/combined.logFC.csv")
row.names(combined.logFC) = combined.logFC$X
combined.logFC = combined.logFC[, -1]

candidates_id = read.csv("inputs/candidate_new.csv")


candidates = candidates_id$gene_id
cand_exp = combined.logFC[candidates,]
cand_exp_t = t(cand_exp)
colnames(cand_exp_t) = candidates_id$id
co_cand = cor(cand_exp_t)

df = read.csv("inputs/new_profile.csv", na.strings=c("","NA"))
summary(df)
colnames(df)
df = df[-c(2,18:21), -c(1,3,71)]
row.names(df) = df$pid
df = df[, -1]
dim(df)

A = df
table(is.na(A))
A[is.na(A)] = 0

summary(A)
A = data.frame(A)
A = as.numeric(A)
A[] <- lapply(A, function(x) as.numeric(as.character(x)))
lapply(A, class)


sim.jac <- matrix(0, nrow=nrow(A), ncol=nrow(A))
rownames(sim.jac) <- rownames(A)
colnames(sim.jac) <- rownames(A)
pairs <- t(combn(1:nrow(A), 2))
for (i in 1:nrow(pairs)){
  num <- sum(sapply(1:ncol(A), function(x)(min(A[pairs[i,1],x],A[pairs[i,2],x]))))
  den <- sum(sapply(1:ncol(A), function(x)(max(A[pairs[i,1],x],A[pairs[i,2],x]))))
  sim.jac[pairs[i,1],pairs[i,2]] <- num/den
  sim.jac[pairs[i,2],pairs[i,1]] <- num/den 
}
sim.jac[which(is.na(sim.jac))] <- 0
diag(sim.jac) <- 1

pdf("outputs/co-evolution_profile.pdf",
    width = 6.13,
    height = 5.60)
pheatmap::pheatmap(sim.jac)
dev.off()


####Mining two aspects

tmp = intersect(colnames(co_cand), colnames(sim.jac))
co_cand_tmp = co_cand[tmp,tmp]
sim_jac_tmp = sim.jac[tmp,tmp]


tmp_hclust = hclust(as.dist(1-co_cand_tmp))

pdf("outputs/co-expression_profile.pdf",
    width = 6.13,
    height = 5.60)
Heatmap(co_cand_tmp, name= "Corr", 
        cluster_columns = tmp_hclust,
        cluster_rows = tmp_hclust,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))

dev.off()


pdf("outputs/co-evoltion_profile_with_row_order.pdf",
    width = 6.13,
    height = 5.60)
Heatmap(sim_jac_tmp, name = "Sim",
        cluster_columns = tmp_hclust,
        cluster_rows = tmp_hclust,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8)
)
dev.off()
