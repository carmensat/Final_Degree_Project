library(topGO)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(readr)
library(ggplot2)

danr_genes <- read_lines("danr_top50_targets.txt")

# Convert gene symbols to Entrez IDs for topGO
gene_map <- AnnotationDbi::select(org.Dm.eg.db, 
                                  keys = danr_genes,
                                  columns = c("ENTREZID", "SYMBOL"),
                                  keytype = "SYMBOL")

all_genes <- keys(org.Dm.eg.db, keytype = "ENTREZID")
gene_list <- factor(as.integer(all_genes %in% gene_map$ENTREZID))
names(gene_list) <- all_genes

# Create topGO data object
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = gene_list,
              geneSelectionFun = function(x) x == 1,
              annot = annFUN.org,
              mapping = "org.Dm.eg.db",
              ID = "entrez")

# Enrichment test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Extract top results
go_table <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)

# View results  
print(go_table)

go_table$classicFisher <- as.numeric(go_table$classicFisher)

go_plot <- ggplot(go_table, aes(x = reorder(Term, -classicFisher), y = -log10(classicFisher))) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "GO Enrichment Analysis", x = "GO Term", y = "-log10(p-value)") +
  theme_minimal(base_size = 25) +
  theme(axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

go_plot

