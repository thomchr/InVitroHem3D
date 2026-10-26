#!/usr/bin/env Rscript
# Loop-constrained vs unconstrained TAD/subTAD/compartment effect sizes.
# For Reviewer #1 point 7 and Reviewer #2 point 3.
#
# Every class is reported RELATIVE to its own conserved/unchanged reference,
# because both modalities carry a large global offset (see the "reference"
# rows) that is a normalisation property, not biology.
#
# Read-only on <DATA_DIR>; all output goes to OUT_DIR.

source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
base <- file.path(DATA_DIR, "hiC")
OUT  <- dirname(dirname(normalizePath(sub("--file=", "", grep("--file=",
          commandArgs(trailingOnly = FALSE), value = TRUE)[1]))))

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

# Hodges-Lehmann shift of `x` relative to `ref`, with a bootstrap CI.
shift <- function(x, ref, B = 2000) {
  x <- x[is.finite(x)]; ref <- ref[is.finite(ref)]
  hl <- median(x, na.rm = TRUE) - median(ref, na.rm = TRUE)
  bs <- replicate(B, median(sample(x, length(x), TRUE)) -
                     median(sample(ref, length(ref), TRUE)))
  c(delta = hl, lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975)))
}

summarise <- function(d, cls, val, ref_level, label) {
  lv <- setdiff(unique(d[[cls]]), ref_level)
  refv <- d[[val]][d[[cls]] == ref_level]
  res <- do.call(rbind, lapply(lv, function(L) {
    v <- d[[val]][d[[cls]] == L]
    s <- shift(v, refv)
    p <- tryCatch(wilcox.test(v, refv)$p.value, error = function(e) NA_real_)
    data.frame(analysis = label, class = L, n = length(v),
               median = median(v, na.rm = TRUE),
               delta_vs_ref = s["delta"], lo = s["lo"], hi = s["hi"],
               wilcox_p = p, row.names = NULL)
  }))
  rbind(data.frame(analysis = label, class = paste0(ref_level, " (reference)"),
                   n = length(refv), median = median(refv, na.rm = TRUE),
                   delta_vs_ref = 0, lo = NA, hi = NA, wilcox_p = NA),
        res[order(-res$delta_vs_ref), ])
}

set.seed(1)
all <- list()

# ---------------------------------------------------------- A/B compartments
g <- as.data.frame(readRDS(file.path(base, "ABcompartment/gene_compartment_df.rds")))
r <- as.data.frame(readRDS(file.path(base, "ABcompartment/region_compartment_df.rds")))

cat("## Are the *_looped classes disjoint from the unconstrained ones?\n")
cat("genes  : nrow", nrow(g), " unique", length(unique(g$Gene_Name)),
    " -> duplicated gene symbols:", sum(duplicated(g$Gene_Name)), "\n")
cat("regions: nrow", nrow(r), " unique", length(unique(r$RegionName)),
    " -> duplicated regions:", sum(duplicated(r$RegionName)), "\n")
dupg <- g[g$Gene_Name %in% g$Gene_Name[duplicated(g$Gene_Name)], ]
cat("duplicated genes with CONFLICTING compartment labels:",
    sum(tapply(dupg$Compartment_Change, dupg$Gene_Name,
               function(z) length(unique(z)) > 1)), "\n\n")

cat("## Total A->B / B->A once looped and unlooped are pooled\n")
pool <- function(d, cc) {
  tab <- table(d[[cc]])
  c(AtoB_all = unname(tab["AtoB"] + tab["AtoB_looped"]),
    AtoB_looped = unname(tab["AtoB_looped"]),
    BtoA_all = unname(tab["BtoA"] + tab["BtoA_looped"]),
    BtoA_looped = unname(tab["BtoA_looped"]))
}
print(rbind(genes = pool(g, "Compartment_Change"), OCRs = pool(r, "RegionType")))
cat("\n")

all$gene_comp <- summarise(g, "Compartment_Change", "Log2FC", "Unchange",
                           "compartment / genes")
all$ocr_comp  <- summarise(r, "RegionType", "Log2FC", "Unchange",
                           "compartment / OCRs")

# ------------------------------------------------------------ TAD and subTAD
tad <- list(
  c("compareTAD/ATAC_TAD_changes.rds",       "type", "Log2FC", "conservedTAD",    "TAD / OCRs (all)"),
  c("compareTAD/ATACloop_TAD_changes.rds",   "type", "Log2FC", "conservedTAD",    "TAD / cREs (looped)"),
  c("compareTAD/RNA_TAD_changes.rds",        "type", "Log2FC", "conservedTAD",    "TAD / genes (all)"),
  c("compareTAD/RNAloop_TAD_changes.rds",    "type", "Log2FC", "conservedTAD",    "TAD / genes (looped)"),
  c("compareTAD/ATAC_subTAD_changes.rds",    "type", "Log2FC", "conservedsubTAD", "subTAD / OCRs (all)"),
  c("compareTAD/ATACloop_subTAD_changes.rds","type", "Log2FC", "conservedsubTAD", "subTAD / cREs (looped)"),
  c("compareTAD/RNA_subTAD_changes.rds",     "type", "Log2FC", "conservedsubTAD", "subTAD / genes (all)"),
  c("compareTAD/RNAloop_subTAD_changes.rds", "type", "Log2FC", "conservedsubTAD", "subTAD / genes (looped)")
)
for (t in tad) {
  d <- as.data.frame(readRDS(file.path(base, t[1])))
  all[[t[5]]] <- summarise(d, t[2], t[3], t[4], t[5])
}

res <- do.call(rbind, all)
res$wilcox_p <- signif(res$wilcox_p, 3)
for (cn in c("median", "delta_vs_ref", "lo", "hi")) res[[cn]] <- round(res[[cn]], 3)
write.csv(res, file.path(OUT, "260827.constrained_vs_unconstrained_effectsizes.csv"),
          row.names = FALSE)

for (a in unique(res$analysis)) {
  cat("\n=====", a, "=====\n")
  print(res[res$analysis == a, -1], row.names = FALSE)
}
cat("\nwrote", file.path(OUT, "260827.constrained_vs_unconstrained_effectsizes.csv"), "\n")
