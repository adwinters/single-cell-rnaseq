
############################################################
# Single-cell RNA-seq workflow
#
# Author: Andrew Winters
#
# Methods
# - Cell Ranger import
# - SoupX contamination correction
# - Seurat clustering
# - SingleR annotation
# - Differential expression analysis
#
# Inputs:
# 10x Genomics Cell Ranger outputs
#
# Outputs:
# UMAPs, cell type annotations, DEG tables
############################################################


library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(future)
library(readr)
library(DoubletFinder)
library(SoupX)
library(diem)
library(tidyverse)
library(DropletUtils)
library(SingleR)
library(celldex)
library(ggplot2)
library(ggbreak)
library(ggrepel)
library(patchwork)
library(cowplot)
library(SummarizedExperiment)
library(scran)
library(RColorBrewer)


setwd("sample_dir/output_files/")

sample_dir <- "path/to/cellranger_output"

sample1_dir <- file.path(
  sample_dir,
  "01_analysis",
  "cellranger_count",
  "sample_control"
)

sample2_dir <- file.path(
  sample_dir,
  "01_analysis",
  "cellranger_count",
  "sample_treatment"
)

sample3_dir <- file.path(
  sample_dir,
  "01_analysis",
  "cellranger_count",
  "ref"
)

# Set seed for reproducibility in parallel processing
set.seed(123)  # Choose a number to set the random seed

# Set up multicore processing and memory limit
future::plan(strategy = 'multicore', workers = 8)
options(future.globals.maxSize = 30 * 1024 ^ 3)

# Set file paths for Sample 1 (Control)
filtered_dir1 <- file.path(sample1_dir, "filtered_feature_bc_matrix")
raw_dir1 <- file.path(sample1_dir, "raw_feature_bc_matrix")

# Set file paths for Sample 2 (Treatment)
filtered_dir2 <- file.path(sample2_dir, "filtered_feature_bc_matrix")
raw_dir2 <- file.path(sample2_dir, "raw_feature_bc_matrix")

# Set file paths for Sample 3 (Reference)
# Reference dataset included for initial clustering
# and annotation support but excluded from
# downstream differential expression analyses.
filtered_dir3 <- file.path(sample3_dir, "SC3_v3_NextGem_DI_CellPlex_Neurons_30K_Brain_1_count_sample_feature_bc_matrix.h5")

# Load the toc and tod data for Sample 1 (Control)
toc1 <- Read10X(filtered_dir1)  # Filtered counts for Sample 1
tod1 <- Read10X(raw_dir1)       # Raw counts for Sample 1

# Load the toc and tod data for Sample 2 (Treatment)
toc2 <- Read10X(filtered_dir2)  # Filtered counts for Sample 2
tod2 <- Read10X(raw_dir2)       # Raw counts for Sample 2

# Load the toc data for Sample 3 (Reference)
toc3 <- Read10X_h5(filtered_dir3)  # Load the filtered counts for Sample 3

# Assuming toc3 is a list and you want to use the first matrix
# Extract the matrix you need from toc3
toc3_matrix <- toc3[[1]]  # Change the index if you want a different matrix

# Rename cell barcodes to make them unique for each sample
colnames(toc1) <- paste0(colnames(toc1), "_Control")
colnames(tod1) <- paste0(colnames(tod1), "_Control")
colnames(toc2) <- paste0(colnames(toc2), "_Treatment")
colnames(tod2) <- paste0(colnames(tod2), "_Treatment")
colnames(toc3_matrix) <- paste0(colnames(toc3_matrix), "_Reference")

# Create a Seurat object using the combined data (Control, Treatment, and Reference samples)
combined_seurat <- CreateSeuratObject(counts = cbind(toc1, toc2, toc3_matrix), min.cells = 3, min.features=200)

# Add metadata for all three samples (Control, Treatment, Reference)
metadata1 <- data.frame(sample = rep("Control", ncol(toc1)))   # Metadata for Sample 1
metadata2 <- data.frame(sample = rep("Treatment", ncol(toc2)))   # Metadata for Sample 2
metadata3 <- data.frame(sample = rep("Reference", ncol(toc3_matrix)))  # Metadata for Sample 3

# Set row names of metadata to match cell barcodes
rownames(metadata1) <- colnames(toc1)
rownames(metadata2) <- colnames(toc2)
rownames(metadata3) <- colnames(toc3_matrix)

# Combine metadata for all samples
combined_metadata <- rbind(metadata1, metadata2, metadata3)

# Add the combined metadata to the Seurat object
combined_seurat <- AddMetaData(combined_seurat, metadata = combined_metadata)

# Perform standard preprocessing in Seurat (normalization, variable features, scaling, PCA)
combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat)

# Perform clustering using Seurat's clustering functions
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:10)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)  # Adjust resolution as needed

# Run UMAP for dimensionality reduction
combined_seurat <- RunUMAP(combined_seurat, dims = 1:10)


# Create SoupChannel objects for both samples
#sc1 <- SoupChannel(tod = tod1, toc = toc1, calcSoupProfile = TRUE)
#sc2 <- SoupChannel(tod = tod2, toc = toc2, calcSoupProfile = TRUE)


# Create SoupChannel objects for both samples with adjusted parameters
sc1 <- SoupChannel(tod = tod1, toc = toc1, calcSoupProfile = TRUE, tfidf.cutoff = 0.2, soupQuantile = 0.75)
sc2 <- SoupChannel(tod = tod2, toc = toc2, calcSoupProfile = TRUE, tfidf.cutoff = 0.5, soupQuantile = 0.95)
#sc3 <- SoupChannel(tod = toc3_matrix, toc = toc3_matrix,
                   #calcSoupProfile = TRUE,
                   #tfidf.cutoff = 0.001,
                   #soupQuantile = 0.75)

# Assign clusters to SoupChannel objects
sc1 <- setClusters(sc1, setNames(combined_seurat$seurat_clusters[1:ncol(toc1)], colnames(toc1)))
sc2 <- setClusters(sc2, setNames(combined_seurat$seurat_clusters[(ncol(toc1)+1):(ncol(toc1)+ncol(toc2))], colnames(toc2)))
#sc3 <- setClusters(sc3, setNames(combined_seurat$seurat_clusters[(ncol(toc1) + ncol(toc2) + 1):ncol(combined_seurat)], colnames(toc3_matrix)))



# Estimate contamination fractions using the clustering information
sc1 <- autoEstCont(
  sc1,
  tfidfMin = 0.5,                 # Minimum tf-idf value for marker genes
  soupQuantile = 0.95,          # Quantile to filter for high expression in soup
  maxMarkers = 200,             # Maximum number of marker genes to consider
  contaminationRange = c(0.01, 0.8),  # Range for contamination fractions
  rhoMaxFDR = 0.2,              # False discovery rate for rho estimation
  priorRho = 0.02,              # Prior mode for contamination fraction
  priorRhoStdDev = 0.1,         # Prior standard deviation
  doPlot = TRUE,                # Generate a plot of density estimates
  forceAccept = FALSE,          # Whether to allow high contamination fractions
  verbose = TRUE                 # Enable verbose output
)


# Perform soup correction for the combined data after estimating contamination fractions
sc1 <- setContaminationFraction(sc1, 0.2)
corrected_counts1 <- adjustCounts(sc1)

sc2 <- setContaminationFraction(sc2, 0.2)
corrected_counts2 <- adjustCounts(sc2)


# Combine corrected counts for Control and Treatment only
combined_corrected_counts <- cbind(corrected_counts1, corrected_counts2)

#save.image(file = "R_workspace_after_combined_corrected_counts.RData")
#load("R_workspace_after_combined_corrected_counts.RData")

# Create metadata for the combined samples (Control and Treatment)
# Ensure row names of combined_metadata match the column names of combined_corrected_counts
combined_metadata <- data.frame(
    sample = c(rep("Control", ncol(corrected_counts1)), 
               rep("Treatment", ncol(corrected_counts2))),
    row.names = c(colnames(corrected_counts1), colnames(corrected_counts2)))

# Create the Seurat object with the combined corrected counts and metadata
combined_seurat <- CreateSeuratObject(counts = combined_corrected_counts, meta.data = combined_metadata)

# Identify cells with zero counts across all features
zero_counts <- Matrix::rowSums(GetAssayData(combined_seurat, slot = "counts")) == 0

# Get the cell names with zero counts
cells_to_remove <- colnames(combined_seurat)[zero_counts]

# Print the number of cells to be removed
cat("Number of cells with zero counts:", length(cells_to_remove), "\n")

# Subset the Seurat object to keep only cells with non-zero counts
combined_seurat <- subset(combined_seurat, cells = colnames(combined_seurat)[!zero_counts])


# Remove the Reference sample if present
combined_seurat_Control_Treatment <- subset(
  combined_seurat,
  subset = sample != "Reference"
)

# Verify the remaining samples
cat("Remaining samples in the Seurat object:", unique(combined_seurat_Control_Treatment$sample), "\n")


# Perform standard preprocessing in Seurat (normalization, variable features, scaling, PCA)
combined_seurat_Control_Treatment <- NormalizeData(combined_seurat_Control_Treatment)
combined_seurat_Control_Treatment <- FindVariableFeatures(combined_seurat_Control_Treatment)
combined_seurat_Control_Treatment <- ScaleData(combined_seurat_Control_Treatment)
combined_seurat_Control_Treatment <- RunPCA(combined_seurat_Control_Treatment)

combined_seurat_Control_Treatment <- RunUMAP(combined_seurat_Control_Treatment, dims = 1:10)



# Cell classification with SingleR
# Use the mouse reference dataset (you can change to another dataset if needed)
# Load a mouse reference dataset from SingleR
mouse_ref <- MouseRNAseqData(ensembl = FALSE, cell.ont = "all")



#save.image(file = "R_workspace_cell_classification.RData")
#load("R_workspace_cell_classification.RData")

# Extract the test data
test_data <- GetAssayData(combined_seurat_Control_Treatment, slot = "counts")

# Check the dimensions of the test data
cat("Number of test cells:", ncol(test_data), "\n")
cat("Number of test genes:", nrow(test_data), "\n")

# Extract labels from mouse_ref
mouse_ref_labels <- colData(mouse_ref)$label.main  # Use the appropriate label column

# Check the length of the labels
cat("Number of reference labels:", length(mouse_ref_labels), "\n")



# Run SingleR classification
singleR_results <- SingleR(test = test_data,
                            ref = mouse_ref,
                            labels = mouse_ref_labels)

# Check results
head(singleR_results)

# Create a new column in singleR_results for Sample Type
# Extracting sample type from the row names of singleR_results
singleR_results$Sample_Type <- ifelse(grepl("Control", rownames(singleR_results)), "Control", "Treatment")

# Ensure the row names of singleR_results match the Seurat object cell names
if (!identical(rownames(singleR_results), colnames(combined_seurat_Control_Treatment))) {
  stop("Row names of singleR_results do not match cell names in combined_seurat_Control_Treatment.")
}

# Add the Cell_Type predictions to the Seurat object
singleR_df <- as.data.frame(singleR_results)
combined_seurat_Control_Treatment$Cell_Type <- singleR_df$pruned.labels

# Replace NA values in Cell_Type with "Unassigned"
combined_seurat_Control_Treatment$Cell_Type[is.na(combined_seurat_Control_Treatment$Cell_Type)] <- "Unassigned"

# Create a summary table of predicted labels, grouped by Cell_Type and Sample_Type
label_summary <- as.data.frame(table(singleR_results$pruned.labels, singleR_results$Sample_Type))
colnames(label_summary) <- c("Cell_Type", "Sample_Type", "Count")

# Print the label_summary to check the new structure
print(label_summary)

# Remove zero counts
label_summary <- label_summary[label_summary$Count > 0, ]

# Create the cell type distribution plot
cell_type_plot <- ggplot(label_summary, aes(x = reorder(Cell_Type, -Count), y = Count, fill = Sample_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_break(c(1100, 2500)) +
  scale_y_break(c(650, 1000)) +
  scale_y_break(c(350, 550)) +
  scale_fill_manual(values = c("Control" = "blue", "Treatment" = "red")) +
  theme_minimal() +
  labs(title = "Cell Type Distribution from SingleR", x = "Cell Type", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to a file
ggsave("cell_type_distribution.png", plot = cell_type_plot, width = 8, height = 6)


# Ensure "Unassigned" is in the cell types
unique_cell_types <- unique(combined_seurat_Control_Treatment$Cell_Type)

# Number of distinct cell types excluding "Unassigned"
n_colors_needed <- length(unique_cell_types) - 1

# Use a combination of "Dark2" and "Set1" to get 13 distinct vibrant colors
distinct_colors <- c(brewer.pal(8, "Dark2"), brewer.pal(9, "Set1")[1:5])

# Combine with gray for "Unassigned"
cell_type_colors <- ifelse(unique_cell_types == "Unassigned", 
                           "gray", 
                           distinct_colors[1:n_colors_needed])



# Create the UMAP plot with conditional coloring and different shapes for cell types
umap_plot_cell_type <- DimPlot(
  combined_seurat_Control_Treatment,
  reduction = "umap",
  group.by = "Cell_Type",
  label = FALSE,
  pt.size = 0.5
) +
  scale_color_manual(values = setNames(cell_type_colors, unique(combined_seurat_Control_Treatment$Cell_Type))) +
  ggtitle("UMAP by Cell Type") +
  theme_minimal()

# Create the UMAP plot with conditional coloring for samples
umap_plot_sample <- DimPlot(
  combined_seurat_Control_Treatment,
  reduction = "umap",
  group.by = "sample",  # Group by sample
  label = FALSE,
  pt.size = 0.5,
  alpha=0.3,
) +
  scale_color_manual(values = c("Control" = "blue", "Treatment" = "red")) +
  ggtitle("UMAP by Sample") +
  theme_minimal()

# Combine the plots into one column with two rows
combined_plot <- umap_plot_cell_type / umap_plot_sample

# Save the combined plot to a file
ggsave("combined_umap_plot.png", plot = combined_plot, width = 10, height = 12, dpi = 300)  # High-resolution plot


#save.image(file = "R_workspace_after_combined_umap.RData")
#load("R_workspace_after_combined_umap.RData")


metadata_df <- as.data.frame(combined_seurat_Control_Treatment@meta.data)
print(metadata_df)



# Check the unique identities in the object
print(unique(combined_seurat_Control_Treatment$sample))           # Check for sample identities
print(unique(combined_seurat_Control_Treatment$Cell_Type))        # Check for cell type identities


# Set active identities to "Cell_Type"
combined_seurat_Control_Treatment

combined_seurat_Control_Treatment <- SetIdent(combined_seurat_Control_Treatment, value = combined_seurat_Control_Treatment$Cell_Type)




# Check the unique identities to confirm they are correctly set
print(unique(Idents(combined_seurat_Control_Treatment)))

# Subset the Seurat object for neurons only
neuron_subset <- subset(combined_seurat_Control_Treatment, idents = "Neurons")



# Check if the subset was successful
print(neuron_subset)

# Normalize, find variable features, and scale the data
neuron_subset <- NormalizeData(neuron_subset)
neuron_subset <- FindVariableFeatures(neuron_subset)
neuron_subset <- ScaleData(neuron_subset)
neuron_subset <- RunPCA(neuron_subset)

# Run UMAP with specified dimensions (choose either 10 or 20 as appropriate)
neuron_subset <- RunUMAP(neuron_subset, dims = 1:20)  # Adjust the number of dimensions as needed

# Create the UMAP DimPlot
neuron_umap_plot <- DimPlot(neuron_subset, reduction = "umap", label = TRUE) + ggtitle("UMAP of Neurons Only")

# Save the DimPlot to a file
ggsave("neuron_umap_plot.png", plot = neuron_umap_plot, width = 8, height = 6, dpi = 300)




# Check unique identities in the original Seurat object for DEGs
unique_idents <- unique(Idents(combined_seurat_Control_Treatment))
print(unique_idents)

# DEGs
# Ensure the identities are correctly set before finding DEGs
deg_neurons <- FindMarkers(object = neuron_subset, 
                           ident.1 = "Control", 
                           ident.2 = "Treatment", 
                           min.pct = 0.1,  
                           test.use = "bimod")


# Convert deg_neurons to a data frame and add a column for absolute fold change
deg_neurons_df <- as.data.frame(deg_neurons) %>%
  mutate(abs_logFC = abs(avg_log2FC))

# Filter for corrected p-values below 0.05 and abs(fold change) > 1.5
filtered_deg_neurons <- deg_neurons_df %>%
  filter(p_val_adj < 0.05 & abs_logFC > 1.5)

# Write the filtered data to a CSV file
write.csv(filtered_deg_neurons, file = "filtered_deg_neurons.csv", row.names = TRUE)

# Select the top markers based on adjusted p-value and absolute log fold change
top_markers <- filtered_deg_neurons %>%
  arrange(p_val_adj, desc(abs_logFC)) %>%
  head(5)  # Adjust this number if you want more markers

# Extract the gene names for plotting
top_gene_names <- rownames(top_markers)


# Generate violin plots for the top markers
vln_plot <- VlnPlot(combined_seurat_Control_Treatment, features = top_gene_names, group.by = "sample") +
  ggtitle("Violin Plots for Top Neuron Markers (p-value < 0.05 & |Fold Change| > 1.5)") +
  theme_minimal()

# Save the violin plot to a file
ggsave("violin_plot_top_neuron_markers.png", plot = vln_plot, width = 8, height = 6)

# Generate feature plots for the top markers
feature_plot <- FeaturePlot(combined_seurat_Control_Treatment, features = top_gene_names, cols = c("lightgrey", "blue")) +
  ggtitle("Feature Plots for Top Neuron Markers (p-value < 0.05 & |Fold Change| > 1.5)") +
  theme_minimal()

# Save the feature plot to a file
ggsave("feature_plot_top_neuron_markers.png", plot = feature_plot, width = 8, height = 6)




#Boxplot
# Generate expression data for the filtered markers
# Ensure the filtered markers exist in the expression data
filtered_genes <- rownames(filtered_deg_neurons)
expression_data <- as.data.frame(GetAssayData(combined_seurat_Control_Treatment, slot = "data")[filtered_genes, ])

# Add sample information correctly
expression_data$cell_id <- rownames(expression_data)  # Create a column for cell identifiers
expression_data <- expression_data %>%
  rownames_to_column(var = "cell_id") %>%
  left_join(data.frame(cell_id = Cells(combined_seurat_Control_Treatment), 
                        sample = combined_seurat_Control_Treatment$sample), 
            by = "cell_id")

# Reshape the data into long format for ggplot2
expression_long <- expression_data %>%
  pivot_longer(cols = -c(cell_id, sample), names_to = "gene", values_to = "expression")

# Generate box plots for the filtered markers
box_plot <- ggplot(expression_long, aes(x = sample, y = expression, fill = sample)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = "free_y") +  # Create separate panels for each gene
  labs(title = "Box Plots for Filtered Neuron Markers (p-value < 0.01 & |Fold Change| > 1.5)",
       x = "Sample Type",
       y = "Expression Level") +
  theme_minimal() +
  theme(legend.position = "none")

# Save the box plot to a file
ggsave(
  "box_plot_filtered_neuron_markers.png",
  plot = box_plot,
  width = 10,
  height = 8,
  dpi = 300
)









markers <- FindAllMarkers(combined_seurat_Control_Treatment, logfc.threshold = 0.75)
head(markers)

# View the markers found
head(markers)

# View the top markers
top_markers <- head(markers[order(markers$p_val_adj), ], n = 8)  # Top 8 markers based on adjusted p-value
print(top_markers)



# Create a violin plot for the top markers
violin_plot <- VlnPlot(combined_seurat_Control_Treatment, 
                        features = rownames(top_markers), 
                        group.by = "sample", 
                        pt.size = 0.1) + 
  ggtitle("Violin Plot of Top Differentially Expressed Genes") +
  theme_minimal()

# Save the violin plot
ggsave("violin_plot_top_markers.png", plot = violin_plot, width = 10, height = 8, dpi = 300)




# Create feature plots for the top markers
feature_plots <- lapply(rownames(top_markers), function(gene) {
  FeaturePlot(combined_seurat_Control_Treatment, features = gene, cols = c("lightgrey", "blue")) +
    ggtitle(gene) +
    theme_minimal()
})

# Save each feature plot to a file
for (i in seq_along(feature_plots)) {
  ggsave(paste0("feature_plot_", rownames(top_markers)[i], ".png"), 
         plot = feature_plots[[i]], 
         width = 10, height = 8, dpi = 300)
}



# Combine the feature plots into one plot
combined_feature_plot <- plot_grid(plotlist = feature_plots, ncol = 2)

# Save the combined feature plot
ggsave("combined_feature_plot.png", plot = combined_feature_plot, width = 15, height = 12, dpi = 300)








