# Author: Carmen Samedi

# ------------------------------------
# This script performs GO Enrichment Analysis
# of danr Target Genes in Drosophila

# --- Load Required Libraries ---
library(topGO)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(readr)
library(ggplot2)

# --- Read Input Gene List ---
# danr_top50_targets.txt: one gene symbol per line (Drosophila symbols)
target_genes <- read_lines("danr_top50_targets.txt")

# --- Map Gene Symbols to Entrez IDs ---
gene_map <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys     = target_genes,
  columns  = c("ENTREZID", "SYMBOL"),
  keytype  = "SYMBOL"
)

# Remove any genes without Entrez ID mapping
entrez_targets <- gene_map$ENTREZID[!is.na(gene_map$ENTREZID)]

# --- Define Gene Universe and Prepare Binary Factor ---
all_entrez <- keys(org.Dm.eg.db, keytype = "ENTREZID")
gene_list <- factor(as.integer(all_entrez %in% entrez_targets))
names(gene_list) <- all_entrez

# --- Initialize topGOdata Object ---
GOdata <- new(
  "topGOdata",
  ontology           = "BP",
  allGenes           = gene_list,
  geneSelectionFun   = function(x) x == 1,
  annot              = annFUN.org,
  mapping            = "org.Dm.eg.db",
  ID                 = "entrez"
)

# --- Run Fisher's Exact Test for GO Enrichment ---
result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# --- Extract and Tidy Top Results ---
top_n <- 10
go_table <- GenTable(
  GOdata,
  classicFisher = result_fisher,
  orderBy       = "classicFisher",
  topNodes      = top_n
)

# Convert p-values to numeric, suppress warnings from "less than" notation
go_table$classicFisher <- as.numeric(sub("<", "", go_table$classicFisher))

# --- Plot Top GO Terms ---
plot <- ggplot(go_table, aes(x = reorder(Term, -classicFisher), y = -log10(classicFisher))) +
  geom_col(fill = "#4477AA", width = 0.7) +
  coord_flip() +
  labs(
    title    = "Top GO Enrichment: danr Target Genes",
    x        = "GO Biological Process",
    y        = expression(-log[10](italic("p-value")))
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.title  = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", hjust = 0.5)
  )

print(plot)
