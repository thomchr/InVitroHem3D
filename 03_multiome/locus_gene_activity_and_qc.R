# Load required libraries
source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
library(Seurat)
library(Signac)
library(dplyr)
obj=readRDS(file.path(DATA_DIR, "snATACseq/HPC_HE_seurat_ATACOnly.GRanges_clean.rds"))
obj$cell.type = factor(obj$cell.type,levels =  c("1_Stromal","2_Endo", "3_HE" ,"4_HPC","5_LMPP_Myeloid","6_LMPP_Lymphoid","7_EMP","8_EMP_Myeloid"  ,"9_EMP_Mk" ,"10_EMP_Ery" ))
# -----------------------------------------------------------------
# 1. Calculate Average Global Accessibility per Cell by Cluster
# -----------------------------------------------------------------
# By default, Seurat/Signac stores the total number of fragments in peaks
# as 'nCount_ATAC' (or 'nCount_peaks') and the number of accessible peaks 
# as 'nFeature_ATAC' in the metadata. 
# Let's extract the metadata and summarize it by cluster
accessibility_metrics <- obj@meta.data %>%
  group_by(cell.type) %>%
  summarize(
    Avg_Fragments_Per_Cell = mean(nCount_ATAC), # Adjust if your column is named differently
    Avg_Peaks_Per_Cell = mean(nFeature_ATAC),
    Cell_Count = n()
  )

print("Average Accessibility per Cell:")
print(accessibility_metrics)

library(ggplot2)

# Generate Violin Plots for Fragments and Peaks per cell
# We set pt.size = 0 to avoid overplotting dots for the large clusters
vln_plot <- VlnPlot(
  object = obj, 
  features = c("nCount_ATAC", "nFeature_ATAC"), 
  group.by = "cell.type", 
  pt.size = 0, 
  ncol = 1
)

# Optional: Add a boxplot inside the violin for better distribution summary
vln_plot <- vln_plot + 
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(vln_plot)

# Save the plot for your rebuttal
ggsave("HE_cluster_accessibility_metrics.pdf", plot = vln_plot, width = 8, height = 10)

# -----------------------------------------------------------------
# 2. Calculate Average Gene Activity Score per Cell by Cluster
# -----------------------------------------------------------------
gene.activities <- GeneActivity(obj)
obj[['Activity']] <- CreateAssayObject(counts = gene.activities)



# -----------------------------------------------------------------
# 1. Calculate Average Gene Activity per Cell
# -----------------------------------------------------------------
# Define the name of your gene activity assay (often "RNA" or "Activity" in Signac)
assay_name <- 'Activity'
DefaultAssay(obj) <- assay_name

# Extract the normalized data matrix
# Note: If you are using Seurat v5, use layer = "data". For Seurat v4, use slot = "data".
activity_matrix <- GetAssayData(obj, assay = assay_name, layer = "data")

# Calculate the average score across all genes for each individual cell
# colMeans calculates the mean of each column (cell)
obj$avg_gene_activity_per_cell <- colMeans(activity_matrix)

# -----------------------------------------------------------------
# 2. Summarize the Metrics by Cluster
# -----------------------------------------------------------------
activity_summary <- obj@meta.data %>%
  group_by(cell.type) %>%
  summarize(
    Mean_Gene_Activity = mean(avg_gene_activity_per_cell),
    Median_Gene_Activity = median(avg_gene_activity_per_cell),
    Cell_Count = n()
  )

print("Cluster-level Summary of Gene Activity:")
print(activity_summary)

# -----------------------------------------------------------------
# 3. Visualize for the Reviewer
# -----------------------------------------------------------------
# Create a Violin Plot to show the distribution of gene activity per cell
vln_activity <- VlnPlot(
  object = obj, 
  features = "avg_gene_activity_per_cell", 
  group.by = "cell.type", 
  pt.size = 0 # Hides individual dots to prevent overplotting on large clusters
) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Average Gene Activity Score per Cell") +
  ylab("Average Gene Activity")

# Display the plot
print(vln_activity)

# Save the plot (optional)
ggsave("HE_cluster_gene_activity.pdf", plot = vln_activity, width = 8, height = 6)

# Plot A: Fragments per cell
p_fragments <- VlnPlot(
  object = obj, 
  features = "nCount_ATAC", 
  group.by = "cell.type", 
  pt.size = 0
) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  coord_flip() + # This swaps the axes!
  theme(legend.position = "none") + # Remove legend to save space
  ggtitle("Fragments per Cell") +
  xlab("Cell Type") + # This will appear on the left (y-axis) due to coord_flip
  ylab("Count") # This will appear on the bottom (x-axis)

# Plot B: Peaks per cell
p_peaks <- VlnPlot(
  object = obj, 
  features = "nFeature_ATAC", 
  group.by = "cell.type", 
  pt.size = 0
) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  coord_flip() + 
  theme(legend.position = "none", axis.title.y = element_blank()) + # Hide y-axis title since it's redundant
  ggtitle("Peaks per Cell") +
  ylab("Count")

# Plot C: Gene Activity per cell
# (Assuming you successfully added the avg_gene_activity_per_cell to your metadata)
p_activity <- VlnPlot(
  object = obj, 
  features = "avg_gene_activity_per_cell", 
  group.by = "cell.type", 
  pt.size = 0
) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  coord_flip() + 
  theme(legend.position = "none", axis.title.y = element_blank()) + 
  ggtitle("Avg Gene Activity per Cell") +
  ylab("Score")

# -----------------------------------------------------------------
# 2. Combine and Display
# -----------------------------------------------------------------

# Use patchwork to place them side-by-side using the '|' operator
combined_plot <- p_fragments | p_peaks | p_activity

# Display the plot
print(combined_plot)

VlnPlot(
  object = obj, 
  features = "avg_gene_activity_per_cell", 
  group.by = "cell.type", 
  pt.size = 0
) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  coord_flip() + 
  theme(legend.position = "none", axis.title.y = element_blank()) + 
  ggtitle("Avg Gene Activity per Cell") +
  ylab("Score")


obj$group = "HE"
obj$group[which(obj$orig.ident == "HPC_seurat_project")] = "HPC"

Idents(obj) <- "group"

obj <- NormalizeData(
  object = obj,
  assay = 'Activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(obj$nCount_ATAC)
)

# 4. Run FindMarkers on this aggregated assay for your target gene
DefaultAssay(obj) <- 'Activity'
a <- FindMarkers(
  object = obj,
  ident.1 = "HE",
  ident.2 = "HPC",
  features = c("VPS37A", "SLC7A2", "PDGFRL", "MTUS1" , "FGL1" )
)
write.table(a,file=file.path(DATA_DIR, "snATACseq/MTUS1Locus_sinal.txt"),col.names=T,row.names=T,sep="\t",quote=F)

a=FindMarkers(
  object = obj,
  ident.1 = "HE",
  ident.2 = "HPC",
  features = c("HAPLN2", "CRABP2", "HDGF" , "NTRK1" )
)
write.table(a,file=file.path(DATA_DIR, "snATACseq/PRCCLocus_sinal.txt"),col.names=T,row.names=T,sep="\t",quote=F)

a=FindMarkers(
  object = obj,
  ident.1 = "HE",
  ident.2 = "HPC",
  features = c("NES", "SH2D2A", "PEAR1" )
)
write.table(a,file=file.path(DATA_DIR, "snATACseq/Rongreuested_sinal.txt"),col.names=T,row.names=T,sep="\t",quote=F)



