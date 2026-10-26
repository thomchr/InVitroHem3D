#!/usr/bin/env Rscript
# Fig. 3e — log2 fold change (HPC / adherent) of gene expression and OCR accessibility by
# compartment class, redrawn from the plotted data. Statistics: fig3e_stats.R.
source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
suppressMessages({ library(ggplot2); library(ggsignif); library(patchwork) })

d   <- file.path(DATA_DIR, "hiC/ABcompartment")
out <- OUT_DIR
dir.create(out, showWarnings = FALSE, recursive = TRUE)

gene <- readRDS(file.path(d, "gene_compartment_df.rds"))
reg  <- readRDS(file.path(d, "region_compartment_df.rds"))

lv   <- c("AtoB", "AtoB_looped", "BtoA", "BtoA_looped", "Unchange")
grp_labels <- parse(text = c("AtoB", "AtoB~looped", "BtoA", "BtoA~looped", "No~Delta"))   # plotmath: pdf() cannot encode Δ
cols <- c(AtoB = "#F8766D", AtoB_looped = "#E88A82", BtoA = "#00BA38", BtoA_looped = "#66C07F", Unchange = "grey75")

# brackets as drawn in the figure, bottom to top; two-sided Wilcoxon rank-sum, unadjusted
pairs <- list(c("AtoB", "AtoB_looped"), c("AtoB", "BtoA"), c("AtoB", "Unchange"),
              c("AtoB_looped", "BtoA_looped"), c("AtoB_looped", "Unchange"),
              c("BtoA", "BtoA_looped"), c("BtoA", "Unchange"), c("BtoA_looped", "Unchange"))

panel <- function(df, ylab) {
  df <- df[is.finite(df$y), ]                    # 815 genes without expression values are excluded
  df$g <- factor(df$g, levels = lv)
  ggplot(df, aes(g, y, fill = g)) +
    geom_boxplot(notch = TRUE, outlier.shape = NA, colour = "grey20", linewidth = 0.4) +
    geom_signif(comparisons = pairs, test = "wilcox.test", map_signif_level = TRUE,
                y_position = seq(5.2, 10.1, length.out = length(pairs)),
                tip_length = 0.008, textsize = 2.6, vjust = 0.2) +
    scale_fill_manual(values = cols, labels = grp_labels, name = NULL) +
    scale_y_continuous(breaks = seq(-5, 7.5, 2.5)) +
    coord_cartesian(ylim = c(-5, 11)) +   # display only; statistics use all values
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(), legend.key.size = unit(0.4, "cm"))
}

p_rna  <- panel(data.frame(g = gene$Compartment_Change, y = gene$Log2FC),
                expression(Delta ~ "Gene Expression (" * log[2] * "FC)"))
p_atac <- panel(data.frame(g = reg$RegionType, y = reg$Log2FC),
                expression(Delta ~ "Accessibility (" * log[2] * "FC)"))

p <- (p_rna | p_atac) + plot_layout(guides = "collect")
for (ext in c("pdf", "svg", "png")) {
  f <- file.path(out, paste0("Fig3e_compartment_boxplots.", ext))
  if (ext == "pdf") { pdf(f, width = 6.2, height = 4.2, useDingbats = FALSE); print(p); dev.off() }
  else ggsave(f, p, width = 6.2, height = 4.2, dpi = 300)
}
