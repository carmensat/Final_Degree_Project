# ------------------------------------
# Lamina Single-Cell RNA-seq
# Author: Carmen Samedi

# ------------------------------------
# This script performs Seurat analysis, DE test and
# GO enrichment plotting
# ------------------------------------

# --- Load required libraries ---
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(arrow)
library(dplyr)
library(plotly)
library(AnnotationHub)
library(GenomicRanges)
library(pheatmap)
library(readxl)
library(cowplot)
library(topGO)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(presto)
library(rtracklayer)

# --- Define file paths and output directory ---
output_dir <- "pdf_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)

# --- Load datasets ---
cistopic_obj <- readRDS("GSE163697_cisTopicObject_OL.Rds")
rna_data    <- readRDS("rna_data.RDS")
lamina      <- readRDS("lamina_subset.RDS")
p_50        <- readRDS("p_50.rds")
lamina_2    <- readRDS("lamina_p50.rds")
lamina_full <- readRDS("lamina_tot.rds")

# --- Set Seurat identities ---
Idents(rna_data) <- "Neuropil"
Idents(lamina)   <- "Neuropil"

# --- Basic Visualization (UMAPs, FeaturePlots) ---
pdf(file.path(output_dir, "rna_dimplot.pdf"), width = 10, height = 8)
DimPlot(rna_data, reduction = "umap", group.by = "Neuropil", label = TRUE) + NoLegend()
dev.off()

genes_of_interest <- c("dan", "danr", "bsh", "erm", "svp", "bab2", "zfh1", "ap", "pdm3")

# Feature plots for candidate markers
pdf(file.path(output_dir, "rna_featureplots.pdf"), width=12, height=8)
FeaturePlot(rna_data, features = genes_of_interest, label = TRUE)
dev.off()

# --- Lamina UMAP and Heatmap Visualizations ---
pdf(file.path(output_dir, "lamina_umap.pdf"), width=10, height=8)
DimPlot(lamina, reduction = "umap", group.by = "Neuropil", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(output_dir, "lamina_heatmap.pdf"), width=10, height=8)
DoHeatmap(lamina, group.by='seurat_clusters', features=genes_of_interest)
dev.off()

# --- Differential Gene Expression & GO Analysis ---
deg <- read.csv("differentially_expressed_genes.csv")
go_results <- read.csv("go_enrichment_results.csv")

# Plot top 10 enriched GO terms
top_n <- 10
go_plot <- ggplot(head(go_results, top_n), 
                  aes(x = reorder(Term, -classicFisher), y = -log10(classicFisher))) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "Top GO Enrichment Analysis", x = "GO Term", y = "-log10(p-value)") +
  theme_minimal()
ggsave(filename = file.path(output_dir, "go_enrichment_barplot.pdf"), plot = go_plot, width = 10, height = 8)

# --- Subset and Reclustering Lamina Cells ---
lamina_2 <- lamina_2 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = genes_of_interest) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:30, reduction = "pca")

# Cluster annotation mapping
cluster_to_lamina <- c("0" = "L1", "1" = "L2", "2" = "L3", "3" = "L4-L5",
                       "4" = "L3", "5" = "L1", "6" = "L4-L5", "7" = "L2")
lamina_2 <- RenameIdents(lamina_2, cluster_to_lamina)

# Save UMAP and marker plots
pdf(file.path(output_dir, "lamina50_umap.pdf"), width=10, height=8)
DimPlot(lamina_2, label = TRUE) + NoLegend()
dev.off()

pdf(file.path(output_dir, "lamina50_heatmap.pdf"), width=10, height=8)
DoHeatmap(lamina_2, features=genes_of_interest)
dev.off()

# --- Differential Expression: danr+ vs danr-, dan+ vs dan- ---
for (gene in c("danr", "dan")) {
  group_col <- paste0(gene, "_group")
  lamina_2[[group_col]] <- ifelse(FetchData(lamina_2, vars=gene)[,1] > 0, paste0(gene, "_high"), paste0(gene, "_low"))
  Idents(lamina_2) <- group_col
  markers <- FindMarkers(lamina_2, ident.1 = paste0(gene, "_high"), ident.2 = paste0(gene, "_low"))
  saveRDS(markers, file = paste0("differentially_expressed_gene_lamina_", gene, ".rds"))
}
