############################################################
# Single-cell RNA-seq workflow
#
# Author: Andrew Winters
#
# Methods
#   - Cell Ranger import
#   - SoupX ambient RNA correction
#   - Seurat normalization, clustering and UMAP
#   - SingleR cell type annotation
#   - Differential expression analysis (neurons, Control vs Treatment)
#
# Inputs
#   10x Genomics Cell Ranger outputs (filtered and raw feature-barcode
#   matrices) for a Control and a Treatment sample, plus a reference
#   sample (.h5) used only for the initial clustering.
#
# Outputs
#   UMAPs, cell type distribution, DEG tables, violin/feature/box plots,
#   neuronal DEG dot plot, App feature plots and a four-panel figure
#
# Requirements
#   R packages listed below; reading the .h5 reference also needs hdf5r.
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(SingleR)
  library(celldex)
  library(SummarizedExperiment)
  library(Matrix)
  library(future)
  library(tidyverse)   # dplyr, tidyr, tibble, ggplot2
  library(ggbreak)
  library(patchwork)
  library(cowplot)
  library(RColorBrewer)
})


############################################################
# 1. Paths and settings
############################################################

base_dir   <- "path/to/cellranger_output"
output_dir <- "path/to/output_files"

counts_dir    <- file.path(base_dir, "01_analysis", "cellranger_count")
control_dir   <- file.path(counts_dir, "sample_control")
treatment_dir <- file.path(counts_dir, "sample_treatment")

# Reference dataset: included for the initial clustering that feeds SoupX,
# excluded from all downstream analyses
reference_h5 <- file.path(
  counts_dir, "ref",
  "SC3_v3_NextGem_DI_CellPlex_Neurons_30K_Brain_1_count_sample_feature_bc_matrix.h5"
)

group_colors <- c(Control = "blue", Treatment = "red")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# Reproducibility and parallel processing
set.seed(123)
future::plan(strategy = "multicore", workers = 8)
options(future.globals.maxSize = 30 * 1024^3)


############################################################
# 2. Load Cell Ranger data
############################################################

# toc = table of counts (filtered cells); tod = table of droplets (raw)
toc1 <- Read10X(file.path(control_dir, "filtered_feature_bc_matrix"))
tod1 <- Read10X(file.path(control_dir, "raw_feature_bc_matrix"))

toc2 <- Read10X(file.path(treatment_dir, "filtered_feature_bc_matrix"))
tod2 <- Read10X(file.path(treatment_dir, "raw_feature_bc_matrix"))

toc3 <- Read10X_h5(reference_h5)
if (is.list(toc3)) toc3 <- toc3[[1]]   # keep the first (Gene Expression) matrix

# Make barcodes unique across samples
colnames(toc1) <- paste0(colnames(toc1), "_Control")
colnames(tod1) <- paste0(colnames(tod1), "_Control")
colnames(toc2) <- paste0(colnames(toc2), "_Treatment")
colnames(tod2) <- paste0(colnames(tod2), "_Treatment")
colnames(toc3) <- paste0(colnames(toc3), "_Reference")


############################################################
# 3. Initial clustering (all three samples, used by SoupX)
############################################################

initial_metadata <- data.frame(
  sample = c(rep("Control",   ncol(toc1)),
             rep("Treatment", ncol(toc2)),
             rep("Reference", ncol(toc3))),
  row.names = c(colnames(toc1), colnames(toc2), colnames(toc3))
)

combined_seurat <- CreateSeuratObject(
  counts       = cbind(toc1, toc2, toc3),
  min.cells    = 3,
  min.features = 200
)
combined_seurat <- AddMetaData(combined_seurat, metadata = initial_metadata)

combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:10)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:10)


############################################################
# 4. Ambient RNA correction (SoupX)
############################################################

sc1 <- SoupChannel(tod = tod1, toc = toc1, calcSoupProfile = TRUE)
sc2 <- SoupChannel(tod = tod2, toc = toc2, calcSoupProfile = TRUE)

# Look up each cell's cluster by barcode. Cells removed by min.features = 200
# in the initial object get their own "filtered" group, because SoupX needs
# a cluster for every cell in toc.
clusters_for <- function(toc) {
  cl <- as.character(combined_seurat$seurat_clusters[colnames(toc)])
  cl[is.na(cl)] <- "filtered"
  setNames(cl, colnames(toc))
}

sc1 <- setClusters(sc1, clusters_for(toc1))
sc2 <- setClusters(sc2, clusters_for(toc2))

# Automatic contamination estimate for the Control sample (reported for
# reference only)
sc1 <- autoEstCont(
  sc1,
  tfidfMin           = 0.5,
  soupQuantile       = 0.95,
  maxMarkers         = 200,
  contaminationRange = c(0.01, 0.8),
  rhoMaxFDR          = 0.2,
  priorRho           = 0.02,
  priorRhoStdDev     = 0.1,
  doPlot             = TRUE,
  forceAccept        = FALSE,
  verbose            = TRUE
)

# A fixed contamination fraction of 0.2 is applied to both samples
sc1 <- setContaminationFraction(sc1, 0.2)
sc2 <- setContaminationFraction(sc2, 0.2)

corrected_counts1 <- adjustCounts(sc1)
corrected_counts2 <- adjustCounts(sc2)


############################################################
# 5. Build the corrected Control + Treatment object
############################################################

combined_corrected_counts <- cbind(corrected_counts1, corrected_counts2)

corrected_metadata <- data.frame(
  sample = c(rep("Control",   ncol(corrected_counts1)),
             rep("Treatment", ncol(corrected_counts2))),
  row.names = colnames(combined_corrected_counts)
)

seurat_obj <- CreateSeuratObject(
  counts    = combined_corrected_counts,
  meta.data = corrected_metadata
)

# Remove cells with zero total counts after correction
zero_count_cells <- Matrix::colSums(GetAssayData(seurat_obj, layer = "counts")) == 0
cat("Number of cells with zero counts:", sum(zero_count_cells), "\n")
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!zero_count_cells])

cat("Samples in the object:", unique(seurat_obj$sample), "\n")

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)


############################################################
# 6. Cell type annotation (SingleR)
############################################################

mouse_ref        <- MouseRNAseqData(ensembl = FALSE, cell.ont = "all")
mouse_ref_labels <- colData(mouse_ref)$label.main

test_data <- GetAssayData(seurat_obj, layer = "counts")
cat("Test cells:", ncol(test_data), "| Test genes:", nrow(test_data),
    "| Reference labels:", length(mouse_ref_labels), "\n")

singleR_results <- SingleR(
  test   = test_data,
  ref    = mouse_ref,
  labels = mouse_ref_labels
)

stopifnot(identical(rownames(singleR_results), colnames(seurat_obj)))

singleR_results$Sample_Type <- ifelse(
  grepl("Control", rownames(singleR_results)), "Control", "Treatment"
)

seurat_obj$Cell_Type <- singleR_results$pruned.labels
seurat_obj$Cell_Type[is.na(seurat_obj$Cell_Type)] <- "Unassigned"


############################################################
# 7. Cell type distribution
############################################################

label_summary <- as.data.frame(
  table(singleR_results$pruned.labels, singleR_results$Sample_Type)
)
colnames(label_summary) <- c("Cell_Type", "Sample_Type", "Count")
label_summary <- label_summary[label_summary$Count > 0, ]
print(label_summary)

cell_type_plot <- ggplot(
  label_summary,
  aes(x = reorder(Cell_Type, -Count), y = Count, fill = Sample_Type)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_break(c(1100, 2500)) +
  scale_y_break(c(650, 1000)) +
  scale_y_break(c(350, 550)) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  labs(title = "Cell Type Distribution from SingleR",
       x = "Cell Type", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cell_type_distribution.png", plot = cell_type_plot, width = 8, height = 6)


############################################################
# 8. UMAPs by cell type and by sample
############################################################

# Distinct colors for each cell type, gray for "Unassigned"
assigned_types <- setdiff(sort(unique(seurat_obj$Cell_Type)), "Unassigned")
palette <- c(brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"))
if (length(assigned_types) > length(palette)) {
  palette <- colorRampPalette(palette)(length(assigned_types))
}
cell_type_colors <- c(
  setNames(palette[seq_along(assigned_types)], assigned_types),
  Unassigned = "gray"
)

umap_plot_cell_type <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by  = "Cell_Type",
  pt.size   = 0.5
) +
  scale_color_manual(values = cell_type_colors) +
  ggtitle("UMAP by Cell Type") +
  theme_minimal()

umap_plot_sample <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by  = "sample",
  pt.size   = 0.5,
  alpha     = 0.3
) +
  scale_color_manual(values = group_colors) +
  ggtitle("UMAP by Sample") +
  theme_minimal()

ggsave("combined_umap_plot.png",
       plot = umap_plot_cell_type / umap_plot_sample,
       width = 10, height = 12, dpi = 300)

print(head(seurat_obj@meta.data))
print(table(seurat_obj$Cell_Type, seurat_obj$sample))


############################################################
# 9. Neurons: subset, UMAP and differential expression
############################################################

Idents(seurat_obj) <- "Cell_Type"

neuron_subset <- subset(seurat_obj, idents = "Neurons")
Idents(neuron_subset) <- "sample"
print(neuron_subset)

neuron_subset <- NormalizeData(neuron_subset)
neuron_subset <- FindVariableFeatures(neuron_subset)
neuron_subset <- ScaleData(neuron_subset)
neuron_subset <- RunPCA(neuron_subset)
neuron_subset <- RunUMAP(neuron_subset, dims = 1:20)

neuron_umap_plot <- DimPlot(neuron_subset, reduction = "umap", label = TRUE,
                            cols = group_colors) +
  ggtitle("UMAP of Neurons Only")

ggsave("neuron_umap_plot.png", plot = neuron_umap_plot, width = 8, height = 6, dpi = 300)

# Positive avg_log2FC = higher in Control
deg_neurons <- FindMarkers(
  object   = neuron_subset,
  ident.1  = "Control",
  ident.2  = "Treatment",
  min.pct  = 0.1,
  test.use = "bimod"
)

filtered_deg_neurons <- deg_neurons %>%
  mutate(abs_logFC = abs(avg_log2FC)) %>%
  filter(p_val_adj < 0.05 & abs_logFC > 1.5)

write.csv(filtered_deg_neurons, file = "filtered_deg_neurons.csv", row.names = TRUE)

filtered_genes <- rownames(filtered_deg_neurons)

top_gene_names <- filtered_deg_neurons %>%
  arrange(p_val_adj, desc(abs_logFC)) %>%
  head(5) %>%
  rownames()


############################################################
# 10. Plots of the significant neuronal DEGs
############################################################

if (length(filtered_genes) == 0) {

  message("No neuronal DEGs pass p_val_adj < 0.05 and |log2FC| > 1.5; skipping DEG plots.")

} else {

  deg_title <- "Top Neuron Markers (adj. p < 0.05 & |log2FC| > 1.5)"

  vln_plot <- VlnPlot(seurat_obj, features = top_gene_names, group.by = "sample",
                      cols = group_colors) &
    theme_minimal()
  vln_plot <- vln_plot + plot_annotation(title = paste("Violin Plots:", deg_title))
  ggsave("violin_plot_top_neuron_markers.png", plot = vln_plot, width = 8, height = 6)

  feature_plot <- FeaturePlot(seurat_obj, features = top_gene_names,
                              cols = c("lightgrey", "blue")) &
    theme_minimal()
  feature_plot <- feature_plot + plot_annotation(title = paste("Feature Plots:", deg_title))
  ggsave("feature_plot_top_neuron_markers.png", plot = feature_plot, width = 8, height = 6)

  # Box plots of all significant neuronal DEGs (log-normalized expression)
  expression_long <- FetchData(seurat_obj, vars = c("sample", filtered_genes)) %>%
    rownames_to_column(var = "cell_id") %>%
    pivot_longer(cols = all_of(filtered_genes),
                 names_to = "gene", values_to = "expression")

  box_plot <- ggplot(expression_long, aes(x = sample, y = expression, fill = sample)) +
    geom_boxplot() +
    facet_wrap(~ gene, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(title = "Box Plots of Significant Neuron DEGs (adj. p < 0.05 & |log2FC| > 1.5)",
         x = "Sample", y = "Expression Level") +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave("box_plot_filtered_neuron_markers.png", plot = box_plot,
         width = 10, height = 8, dpi = 300)
}


############################################################
# 11. Multipanel figure
#   A. UMAP by cell type
#   B. Dot plot of the top 10 neuronal DEGs (Control vs Treatment)
#   C. UMAP by sample
#   D. App expression in neurons, Control (left) and Treatment (right)
############################################################

# Non-CNS / unassigned labels are left out of the figure only
figure_drop_types <- c("B cells", "NK cells", "T cells", "Granulocytes", "Unassigned")
figure_obj <- seurat_obj[, !seurat_obj$Cell_Type %in% figure_drop_types]

# Neurons on the full-dataset UMAP (neuron_subset has its own re-run UMAP)
neurons_figure <- seurat_obj[, seurat_obj$Cell_Type == "Neurons"]

# Panel A: UMAP by cell type
panel_a <- DimPlot(figure_obj, reduction = "umap", group.by = "Cell_Type", pt.size = 0.5) +
  ggtitle(NULL) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 4)))

# Panel C: UMAP by sample
panel_c <- DimPlot(figure_obj, reduction = "umap", group.by = "sample",
                   pt.size = 0.5, alpha = 0.3, cols = group_colors) +
  ggtitle(NULL) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

# Panel B: dot plot of the top 10 neuronal DEGs (unscaled average expression)
top_10_genes <- head(filtered_genes, 10)

neurons_figure$Cell_Sample <- factor(
  paste("Neurons", neurons_figure$sample, sep = " - "),
  levels = c("Neurons - Treatment", "Neurons - Control")   # Control drawn on top
)

if (length(top_10_genes) > 0) {
  panel_b <- DotPlot(neurons_figure, features = top_10_genes,
                     group.by = "Cell_Sample", scale = FALSE) +
    xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          legend.position = "bottom")

  ggsave("dot_plot_neurons_top10_by_sample.png", plot = panel_b,
         width = 12, height = 4.25, dpi = 300)
} else {
  panel_b <- plot_spacer()
  message("No neuronal DEGs pass the cutoffs; dot plot panel left empty.")
}

# Panel D: App feature plots, one per sample, on a shared color scale
if ("App" %in% rownames(neurons_figure)) {
  app_max <- max(FetchData(neurons_figure, vars = "App"))

  app_plots <- lapply(names(group_colors), function(s) {
    suppressMessages(
      FeaturePlot(neurons_figure[, neurons_figure$sample == s],
                  features = "App", reduction = "umap", pt.size = 0.5) +
        scale_colour_gradient(low = "lightgrey", high = "blue",
                              limits = c(0, app_max), name = "Average Expression") +
        ggtitle(NULL)
    )
  })

  panel_d <- wrap_plots(app_plots, ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  ggsave("feature_plot_App_by_sample.png", plot = panel_d,
         width = 12, height = 6, dpi = 300)
} else {
  panel_d <- plot_spacer()
  message("App not found in the data; App panel left empty.")
}

multipanel_figure <- (panel_a | panel_b) / (panel_c | panel_d)

ggsave("multipanel_figure.png", plot = multipanel_figure,
       width = 16, height = 13.5, dpi = 300)


############################################################
# 12. Cell type markers (all cell types)
############################################################

# Idents are still Cell_Type here
markers <- FindAllMarkers(seurat_obj, logfc.threshold = 0.75)
print(head(markers))

top_markers <- head(markers[order(markers$p_val_adj), ], n = 8)
print(top_markers)

top_marker_genes <- unique(top_markers$gene)

violin_plot <- VlnPlot(seurat_obj, features = top_marker_genes,
                       group.by = "sample", pt.size = 0.1, cols = group_colors) &
  theme_minimal()
violin_plot <- violin_plot +
  plot_annotation(title = "Violin Plot of Top Differentially Expressed Genes")

ggsave("violin_plot_top_markers.png", plot = violin_plot, width = 10, height = 8, dpi = 300)

feature_plots <- lapply(top_marker_genes, function(gene) {
  FeaturePlot(seurat_obj, features = gene, cols = c("lightgrey", "blue")) +
    ggtitle(gene) +
    theme_minimal()
})
names(feature_plots) <- top_marker_genes

for (gene in top_marker_genes) {
  ggsave(paste0("feature_plot_", gene, ".png"), plot = feature_plots[[gene]],
         width = 10, height = 8, dpi = 300)
}

combined_feature_plot <- plot_grid(plotlist = feature_plots, ncol = 2)
ggsave("combined_feature_plot.png", plot = combined_feature_plot,
       width = 15, height = 12, dpi = 300)


############################################################
# Session info
############################################################

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
