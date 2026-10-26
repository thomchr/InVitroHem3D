# Load necessary libraries
source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(GenomicRanges)
library(Matrix)
library(data.table)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

metadata=read.delim(file.path(DATA_DIR, "snATACseq/HPC_HE_seurat.metadata.txt"),header=T,sep="\t")
he_frag_file <- file.path(DATA_DIR, "cellranger_arc/s/ChristopherThom-09272023-HE_Multiome/outs/atac_fragments.tsv.gz")
hpc_frag_file <- file.path(DATA_DIR, "cellranger_arc/s/ChristopherThom-09272023-HPC-pro_Multiome/outs/atac_fragments.tsv.gz")
he_frag <- CreateFragmentObject(he_frag_file)
hpc_frag <- CreateFragmentObject(hpc_frag_file)
peaks <- read.table(file.path(DATA_DIR, "snATACseq/macs2_out/consensus_regions.bed"), header=F)
colnames(peaks)[1:3] <- c("seqnames", "start", "end")
peaks <- GRanges(peaks[,1:3])

# Split metadata and get cell lists for each sample
hpc_metadata <- metadata %>% dplyr::filter(grepl("^HPC_", rownames(.)))
he_metadata <- metadata %>% dplyr::filter(grepl("^HE_", rownames(.)))

# Get prefixed cell barcodes for each sample
prefixed_hpc_cells <- rownames(hpc_metadata)
prefixed_he_cells <- rownames(he_metadata)

# Get unprefixed cell barcodes for FeatureMatrix (assuming fragment files use unprefixed names)
unprefixed_hpc_cells <- gsub("^HPC_", "", prefixed_hpc_cells)
unprefixed_he_cells <- gsub("^HE_", "", prefixed_he_cells)


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'
genome(annotation) <- "hg38"


# --- Process HPC Sample ---
message("Processing HPC sample...")

# Create FeatureMatrix for HPC using only HPC fragment file and HPC cells
# Note: We pass only the HPC fragment file here
hpc_count_matrix <- FeatureMatrix(
  fragments = hpc_frag,
  features = peaks,
  cells = unprefixed_hpc_cells # Use only unprefixed HPC barcodes
)

# Do NOT manually rename the columns here yet! Keep them unprefixed.
# Create Chromatin Assay for HPC
hpc_chrom_assay <- CreateChromatinAssay(
  counts = hpc_count_matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  annotation = annotation,
  fragments = hpc_frag # <-- CRUCIAL ADDITION: This links the file to the assay
)

# 3. Temporarily remove prefixes from metadata so it matches the unprefixed cells
hpc_metadata_unprefixed <- hpc_metadata
rownames(hpc_metadata_unprefixed) <- unprefixed_hpc_cells

# Create Seurat object for HPC
hpc_seurat <- CreateSeuratObject(
  counts = hpc_chrom_assay,
  assay = "ATAC",
  meta.data = hpc_metadata_unprefixed
)

# 5. Safely add the prefixes. This updates the Seurat object AND Signac's fragment mapping!
hpc_seurat <- RenameCells(hpc_seurat, new.names = prefixed_hpc_cells)
message("HPC sample processed.")


# --- Process HE Sample ---
message("Processing HE sample...")

# Create FeatureMatrix for HE using only HE fragment file and HE cells
# Note: We pass only the HE fragment file here
he_count_matrix <- FeatureMatrix(
  fragments = he_frag,
  features = peaks,
  cells = unprefixed_he_cells
)

# 2. Create Chromatin Assay WITH the fragments object
he_chrom_assay <- CreateChromatinAssay(
  counts = he_count_matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  annotation = annotation,
  fragments = he_frag # <-- CRUCIAL ADDITION
)

# 3. Temporarily remove prefixes from metadata 
he_metadata_unprefixed <- he_metadata
rownames(he_metadata_unprefixed) <- unprefixed_he_cells

# 4. Create Seurat object
he_seurat <- CreateSeuratObject(
  counts = he_chrom_assay,
  assay = "ATAC",
  meta.data = he_metadata_unprefixed
)

# 5. Safely add the prefixes back
he_seurat <- RenameCells(he_seurat, new.names = prefixed_he_cells)
message("HE sample processed.")


# --- Merge Seurat Objects ---
message("Merging Seurat objects...")

# Merge the two Seurat objects
# Because your cell barcodes are already prefixed, the merge function
# will correctly keep them distinct (e.g., HPC_CELL001 is different from HE_CELL001)
merged_seurat <- merge(x = hpc_seurat, y = he_seurat, project = "CombinedATAC")

message("Merging complete.")

# Set genome assembly for the merged object

# --- Display Summary ---
message("\nMerged Seurat object created:")
print(merged_seurat)
message("\nMetadata preview (merged):")
print(head(merged_seurat@meta.data))
message("\nAssay features (peaks) preview:")
print(head(rownames(merged_seurat)))

merged_seurat <- RunTFIDF(merged_seurat)
merged_seurat <- FindTopFeatures(merged_seurat, min.cutoff = 'q0')
merged_seurat <- RunSVD(merged_seurat)
merged_seurat <- RunUMAP(object = merged_seurat, reduction = 'lsi', dims = 2:30)

gene.activities <- GeneActivity(merged_seurat, biotypes = NULL,assay = "ATAC")
merged_seurat[['Activity']] <- CreateAssayObject(counts = gene.activities)
activity_matrix <- GetAssayData(merged_seurat, assay = 'Activity', layer = "data")
merged_seurat$avg_gene_activity_per_cell <- colMeans(activity_matrix)

saveRDS(merged_seurat,file=file.path(DATA_DIR, "snATACseq/HPC_HE_seurat_ATACOnly.GRanges_clean.rds"))

library(dplyr)

# Calculate the total nCount_ATAC per cell type
cluster_depths <- obj@meta.data %>%
  group_by(cell.type) %>%
  summarise(Total_Reads = sum(nCount_ATAC)) %>%
  ungroup()

# Calculate the scaling factor for each
# (Targeting 10 million reads as the baseline for visualization)
target_depth <- 10000000 
cluster_depths <- cluster_depths %>%
  mutate(Scale_Factor = target_depth / Total_Reads)

# Print the table to copy into your bash script
print(cluster_depths[,c(1,3)])
