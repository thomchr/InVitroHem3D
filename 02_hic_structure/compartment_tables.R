#!/usr/bin/env Rscript
# Per-gene and per-OCR log2 fold change (HPC / adherent) by compartment class, including the
# loop-engaged ("_looped") subsets from compartment_classes.sh.
#   group averages of log-normalised expression (RNA) or TF-IDF accessibility (ATAC) over all cells
#   in the group, log2FC = log2((HPC + 1e-9) / (HE + 1e-9))
# Writes hiC/ABcompartment/gene_compartment_df.rds and region_compartment_df.rds (inputs of
# fig3e_plot.R and fig3e_stats.R).
source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
suppressMessages({ library(Seurat); library(SeuratDisk); library(Signac); library(Matrix); library(data.table) })

cls <- file.path(DATA_DIR, "hiC/class_peaks")
out <- file.path(DATA_DIR, "hiC/ABcompartment")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

he_celltypes  <- c("1_Stromal", "2_Endo", "3_HE")
hpc_celltypes <- c("4_HPC", "5_LMPP_Myeloid", "6_LMPP_Lymphoid", "7_EMP", "8_EMP_Myeloid", "9_EMP_Mk", "10_EMP_Ery")
levels_cc     <- c("AtoB", "AtoB_looped", "BtoA", "BtoA_looped", "Unchange")

group_means <- function(mat, celltype) {
  grp <- ifelse(celltype %in% he_celltypes, "HE", ifelse(celltype %in% hpc_celltypes, "HPC", NA))
  keep <- !is.na(grp)
  mat <- mat[, keep]; grp <- grp[keep]
  data.table(feature = rownames(mat),
             HE  = Matrix::rowSums(mat[, grp == "HE"])  / sum(grp == "HE"),
             HPC = Matrix::rowSums(mat[, grp == "HPC"]) / sum(grp == "HPC"))
}
log2fc <- function(he, hpc) log2((hpc + 1e-9) / (he + 1e-9))

# ---- genes ------------------------------------------------------------------------------------
rna  <- LoadH5Seurat(file.path(DATA_DIR, "multiome/HPC_HE_seurat_RNAOnly.h5seurat"))
meta <- read.table(file.path(DATA_DIR, "multiome/HPC_HE_seurat.metadata.txt"), header = TRUE)
rna$cell.type <- meta$cell.type
expr <- group_means(GetAssayData(rna, assay = "RNA", slot = "data"), rna$cell.type)
rm(rna); gc()

promoter_genes <- function(f) unique(read.table(f, header = TRUE)$GENENAME)
looped_genes   <- function(c) unique(fread(file.path(cls, paste0(c, ".gene2OCR.txt")), header = FALSE)$V3)  # gene_name

gene_classes <- rbind(
  data.table(Gene_Name = promoter_genes(file.path(DATA_DIR, "hiC/AtoB.promoter.anno.txt")), Compartment_Change = "AtoB"),
  data.table(Gene_Name = looped_genes("AtoB"),                                             Compartment_Change = "AtoB_looped"),
  data.table(Gene_Name = promoter_genes(file.path(DATA_DIR, "hiC/BtoA.promoter.anno.txt")), Compartment_Change = "BtoA"),
  data.table(Gene_Name = looped_genes("BtoA"),                                             Compartment_Change = "BtoA_looped"))
gene_classes <- rbind(gene_classes,
  data.table(Gene_Name = setdiff(expr$feature, gene_classes$Gene_Name), Compartment_Change = "Unchange"))

gene <- merge(gene_classes, expr, by.x = "Gene_Name", by.y = "feature", all.x = TRUE, sort = FALSE)
setnames(gene, c("HE", "HPC"), c("HE_Expr", "HPC_Expr"))
gene[, Log2FC := log2fc(HE_Expr, HPC_Expr)]
gene[, Compartment_Change := factor(Compartment_Change, levels = levels_cc)]
setorder(gene, Compartment_Change)
saveRDS(as.data.frame(gene), file.path(out, "gene_compartment_df.rds"))

# ---- OCRs -------------------------------------------------------------------------------------
atac <- readRDS(file.path(DATA_DIR, "snATACseq/HPC_HE_seurat_ATACOnly.GRanges_clean.rds"))
acc  <- group_means(GetAssayData(atac, assay = "ATAC", layer = "data"), atac$cell.type)
rm(atac); gc()

read_peaks <- function(c) {
  b <- fread(file.path(cls, paste0(c, ".peaks.bed")), header = FALSE)
  data.table(RegionName = paste(b$V1, b$V2, b$V3, sep = "-"), RegionType = c)
}
region_classes <- rbindlist(lapply(levels_cc, read_peaks))
region <- merge(region_classes, acc, by.x = "RegionName", by.y = "feature", all.x = TRUE, sort = FALSE)
region[, Log2FC := log2fc(HE, HPC)]
region[, RegionType := factor(RegionType, levels = levels_cc)]
setorder(region, RegionType)
saveRDS(as.data.frame(region), file.path(out, "region_compartment_df.rds"))
