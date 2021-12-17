

######################################################################
            # gene ontology enrichment analysis #
######################################################################

# required packages

load.lib<-c("dplyr", "openxlsx", "ggplot2", "BiocManager", "UpSetR", "openxlsx")

install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

BiocManager::install(c("edgeR", "ComplexHeatmap", "topGO"))

bioc_pckg = c("edgeR", "ComplexHeatmap", "topGO") 

lapply(bioc_pckg, require, character.only = TRUE)



gene.ID.GO <- readMappings(file = "inputs/final_go_terms.txt")
gene.universe <- names(gene.ID.GO)
gene.universe
genes.of.interest <-  read.csv("inputs/candidate_new.csv")
genes.of.interest$gene_id = gsub("[.].*", "", genes.of.interest$gene_id)
genes.of.interest <- as.character(genes.of.interest$gene_id)
genes.of.interest
gene.list <-
  factor(as.integer(gene.universe %in% genes.of.interest))

gene.list
names(gene.list) <- gene.universe
my.GO.data <- new(
  "topGOdata",
  description = "Tomato",
  ontology = "BP",
  allGenes = gene.list,
  nodeSize = 100,
  annot = annFUN.gene2GO,
  gene2GO = gene.ID.GO
)
allGO = genesInTerm(my.GO.data)
my.GO.data
sig.genes <- sigGenes(my.GO.data)
str(sig.genes)
result.fisher <-
  runTest(my.GO.data, algorithm = "classic", statistic = "fisher")

all.res <- GenTable(
  my.GO.data,
  classicFisher = result.fisher,
  orderBy = "Fisher",
  ranksOf = "classicFisher",
  topNodes = 20,
  numChar = 50
)

p.adj=round(p.adjust(all.res$classicFisher,method="BH"), digits = 40)
all.res = data.frame(all.res, p.adj)
write.csv(all.res, "outputs/GO_candidates.csv", row.names = F)


options(scipen = 999)
options("scipen"=1000, "digits"=8)
ntop <- 15
ggdata <- all.res[1:ntop,]
names(ggdata)[names(ggdata) == "classicFisher"] <- "p_value"
names(ggdata)[names(ggdata) == "Significant"] <- "Counts"

pdf("outputs/GO_enrichment.pdf",
    width = 6.13,
    height = 5.60)

ggplot(ggdata,
           aes(x = Counts, y = Term, size=Counts, color = -log10(p.adj))) +
    geom_point()+
    xlab('Counts') + ylab('') +
    labs(
      title = 'Candidate genes GO biological process')+
  theme(axis.text.y = element_text(angle = 0, size = 10, face = 'bold', vjust = 0.5))
dev.off()
  
