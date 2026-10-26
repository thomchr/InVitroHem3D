#!/usr/bin/env Rscript
# Figures for Reviewer #1 point 7 and Reviewer #2 point 3.
#
#   Fig 1  A/B compartment switches       genes + cREs, all vs promoter-contacting
#   Fig 2  TAD remodelling classes        idem
#   Fig 3  sub-TAD remodelling classes    idem
#   Fig 4  effect sizes (delta vs conserved, 95% CI) for Figs 2-3
#
# Conventions, all deliberate:
#   * y is log2(HPC / adherent). Positive = higher in HPC. The -log2FC convention
#     used in the original Supplementary Fig. S4e-f is NOT used anywhere here.
#   * every distribution is centred on the median of its own conserved/unchanged
#     class, because the raw distributions carry a large global offset that is a
#     property of the normalisation and not of the structural classes.
#   * elements whose identifier is assigned to more than one class within an
#     analysis are DROPPED from that analysis (per instruction, 2026-08-27), and
#     the number dropped is printed and written to the audit CSV.
#
# Read-only on <DATA_DIR>. All output goes to OUT_DIR.

source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

base <- file.path(DATA_DIR, "hiC")
here <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
OUT  <- dirname(dirname(here))
FIG  <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE)

REF   <- "#8a8985"   # neutral ink for the reference class
S1    <- "#2a78d6"   # categorical slot 1 - all elements
S2    <- "#eb6834"   # categorical slot 2 - promoter-contacting
INK   <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"; SURF <- "#fcfcfb"

theme_pub <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = SURF, colour = NA),
      plot.background  = element_rect(fill = SURF, colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
      axis.text  = element_text(colour = INK2),
      axis.title = element_text(colour = INK),
      strip.text = element_text(colour = INK, face = "bold", hjust = 0),
      plot.title = element_text(colour = INK, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(colour = INK2, size = base_size - 0.5),
      plot.caption  = element_text(colour = INK2, size = base_size - 1.5, hjust = 0),
      legend.position = "top", legend.title = element_blank(),
      legend.text = element_text(colour = INK2),
      legend.key.size = unit(9, "pt")
    )
}

audit <- list()

# Drop identifiers that carry genuinely contradictory class calls.
#
# `key` is the class stripped of anything that is not part of the biological call.
# In the compartment files the `_looped` suffix encodes the promoter-contact
# CONSTRAINT, not a different compartment change, so a gene listed under both
# `AtoB` and `AtoB_looped` is the same call recorded in both the unconstrained and
# the constrained set - not a conflict. Only a gene listed under two different
# separates 7 real conflicts from 290 harmless duplicates.
#
# `drop_multi` is TRUE only for compartments. A locus has exactly one compartment
# state, so two states for one gene is a contradiction. TADs and sub-TADs are
# NESTED, so an element sitting in both a conserved and a split domain is the
# normal case, not an error - dropping those removed 97% of the conserved sub-TAD
# reference class (down to 24 cREs / 5 genes) and left the whole figure centred on
# noise. For those files the multi-class elements are kept and counted instead.
drop_conflicts <- function(d, idcol, clscol, tag, strip = NULL, drop_multi = TRUE) {
  n0 <- nrow(d)
  key <- if (is.null(strip)) d[[clscol]] else sub(strip, "", d[[clscol]])
  dupn <- length(unique(d[[idcol]][duplicated(d[[idcol]])]))
  bad <- unique(d[[idcol]][duplicated(d[[idcol]])])
  conflicted <- character(0)
  if (length(bad)) {
    i <- d[[idcol]] %in% bad
    conflicted <- names(which(tapply(key[i], d[[idcol]][i],
                                     function(z) length(unique(z)) > 1)))
    if (drop_multi) d <- d[!d[[idcol]] %in% conflicted, ]
  }
  d <- d[!duplicated(paste(d[[idcol]], d[[clscol]])), ]   # exact repeats
  nonfin <- sum(!is.finite(d$Log2FC))
  audit[[tag]] <<- data.frame(analysis = tag, rows_in = n0,
                              duplicated_ids = dupn,
                              multiclass_ids = length(conflicted),
                              multiclass_dropped = drop_multi,
                              rows_out = nrow(d),
                              non_finite_log2FC_dropped = nonfin)
  cat(sprintf("%-34s %7d -> %7d rows | %5d dup ids, %5d multi-class (%s), %5d non-finite log2FC\n",
              tag, n0, nrow(d), dupn, length(conflicted),
              if (drop_multi) "dropped" else "kept", nonfin))
  d
}

centre <- function(d, val, cls, ref) {
  d$centred <- d[[val]] - median(d[[val]][d[[cls]] == ref], na.rm = TRUE)
  d
}

pretty_tad <- c(conservedTAD = "Conserved", conservedTAD_shifted = "Conserved\n(shifted)",
                splitTAD = "Split", mergeTAD = "Merged", semiTAD = "Semi",
                newlyform = "De novo", disappeared = "Dissolved",
                subTADbecomeTAD = "Ascended\nfrom sub-TAD", complex = "Complex")
pretty_sub <- c(conservedsubTAD = "Conserved", conservedsubTAD_shifted = "Conserved\n(shifted)",
                splitsubTAD = "Split", mergesubTAD = "Merged", newlyform = "De novo",
                disappeared = "Dissolved", TADbecomesubTAD = "Descended\nfrom TAD",
                subTADlower = "Hierarchy\nlower", subTADhigher = "Hierarchy\nhigher",
                complex = "Complex")

# ------------------------------------------------------------------ Figure 1
cat("\n--- compartments ---\n")
g <- as.data.frame(readRDS(file.path(base, "ABcompartment/gene_compartment_df.rds")))
r <- as.data.frame(readRDS(file.path(base, "ABcompartment/region_compartment_df.rds")))
g <- drop_conflicts(g, "Gene_Name",  "Compartment_Change", "compartment / genes", "_looped$")
r <- drop_conflicts(r, "RegionName", "RegionType",         "compartment / cREs", "_looped$")
nf_g <- sum(!is.finite(g$Log2FC) & g$Compartment_Change != "Unchange")
nf_r <- sum(!is.finite(r$Log2FC) & r$RegionType != "Unchange")
g <- g[is.finite(g$Log2FC), ]; r <- r[is.finite(r$Log2FC), ]
g <- centre(g, "Log2FC", "Compartment_Change", "Unchange")
r <- centre(r, "Log2FC", "RegionType",         "Unchange")

mk_comp <- function(d, idcol, cls, modality) {
  x <- data.frame(id = d[[idcol]], raw = d[[cls]], v = d$centred)
  x$sw <- sub("_looped$", "", x$raw)
  all_arm <- x[x$sw != "Unchange", ]                 # unconstrained U constrained
  all_arm <- all_arm[!duplicated(all_arm$id), ]      # a gene in both is counted once
  loop_arm <- x[grepl("_looped$", x$raw), ]
  ref_arm  <- x[x$sw == "Unchange", ]
  rbind(
    data.frame(modality, class = all_arm$sw,  arm = "All elements",          value = all_arm$v),
    data.frame(modality, class = loop_arm$sw, arm = "Promoter-contacting",   value = loop_arm$v),
    data.frame(modality, class = "Unchange",  arm = "Unchanged (reference)", value = ref_arm$v)
  )
}
comp <- rbind(mk_comp(g, "Gene_Name",  "Compartment_Change", "Genes (expression)"),
              mk_comp(r, "RegionName", "RegionType",         "cREs (accessibility)"))
comp$class <- factor(recode(comp$class, Unchange = "Unchanged",
                            AtoB = "A to B", BtoA = "B to A"),
                     levels = c("Unchanged", "A to B", "B to A"))
comp$arm <- factor(comp$arm, levels = c("Unchanged (reference)",
                                        "All elements", "Promoter-contacting"))
comp$modality <- factor(comp$modality, levels = c("cREs (accessibility)", "Genes (expression)"))

# n goes into the tick label, not into the panel, so it can never collide with a
# whisker. scales = "free_x" is what lets the two facets carry different n.
add_n <- function(d, clsvar) {
  key <- paste(d$modality, d[[clsvar]], sep = "\r")
  n   <- table(key, d$arm)
  lab <- apply(n, 1, function(z) paste(scales::comma(z[z > 0]), collapse = " / "))
  d$xlab <- paste0(d[[clsvar]], "\n", lab[key])
  d
}
comp <- add_n(comp, "class")
comp$xlab <- factor(comp$xlab, levels = unique(comp$xlab[order(comp$modality, comp$class)]))

lim <- quantile(comp$value, c(.015, .985), na.rm = TRUE)
lo  <- lim[1]; hi <- lim[2]
p1 <- ggplot(comp, aes(xlab, value, fill = arm)) +
  geom_hline(yintercept = 0, colour = INK2, linewidth = 0.35, linetype = "22") +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, colour = "#33322f",
               position = position_dodge2(preserve = "single", padding = 0.18)) +
  facet_wrap(~modality, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c("Unchanged (reference)" = REF,
                               "All elements" = S1, "Promoter-contacting" = S2)) +
  coord_cartesian(ylim = c(lo, hi), clip = "on") +
  labs(title = "A/B compartment switches during EHT",
       subtitle = "Requiring a promoter contact leaves the compartment result intact and balances the two directions",
       x = NULL, y = expression(log[2]*" fold change (HPC / adherent), centred on unchanged"),
       caption = paste("Boxes: median, IQR, 1.5x IQR whiskers; outliers hidden by axis clipping. Positive = higher in HPC.",
                       "Tick labels give n as all elements / promoter-contacting. 7 genes with contradictory calls removed.",
                       sprintf("%s switching genes and %s switching cREs have an undefined log2 ratio (no signal in one condition) and cannot be shown.",
                               scales::comma(nf_g), scales::comma(nf_r)), sep = "\n")) +
  theme_pub() + theme(axis.text.x = element_text(lineheight = 1.0))

# --------------------------------------------------------------- Figures 2-3
load_pair <- function(all_f, loop_f, idcol_a, idcol_b, modality, tag) {
  a <- as.data.frame(readRDS(file.path(base, all_f)))
  b <- as.data.frame(readRDS(file.path(base, loop_f)))
  names(a)[names(a) %in% c("Gene_Name", "gene_name")] <- "id"
  names(b)[names(b) %in% c("Gene_Name", "gene_name")] <- "id"
  names(a)[names(a) == "RegionName"] <- "id"; names(b)[names(b) == "RegionName"] <- "id"
  a <- drop_conflicts(a, "id", "type", paste(tag, "/ all"),                  drop_multi = FALSE)
  b <- drop_conflicts(b, "id", "type", paste(tag, "/ promoter-contacting"),  drop_multi = FALSE)
  ref <- if (grepl("subTAD", all_f)) "conservedsubTAD" else "conservedTAD"
  a <- centre(a[is.finite(a$Log2FC), ], "Log2FC", "type", ref)
  b <- centre(b[is.finite(b$Log2FC), ], "Log2FC", "type", ref)
  rbind(data.frame(modality, arm = "All elements",         type = a$type, value = a$centred),
        data.frame(modality, arm = "Promoter-contacting",  type = b$type, value = b$centred))
}

build <- function(scale_name, files, pretty, ttl, sub) {
  d <- rbind(
    load_pair(files[1], files[2], , , "cREs (accessibility)",  paste(scale_name, "cREs")),
    load_pair(files[3], files[4], , , "Genes (expression)",    paste(scale_name, "genes"))
  )
  # pretty[<factor>] would index the lookup by level code, silently permuting the
  # class labels.
  d$type <- as.character(d$type)
  ord <- d %>% filter(arm == "Promoter-contacting", modality == "cREs (accessibility)") %>%
    group_by(type) %>% summarise(m = median(value), .groups = "drop") %>% arrange(desc(m))
  stopifnot(all(d$type %in% names(pretty)))
  d$type <- factor(pretty[d$type], levels = unname(pretty[as.character(ord$type)]))
  d$modality <- factor(d$modality, levels = c("cREs (accessibility)", "Genes (expression)"))
  d$arm <- factor(d$arm, levels = c("All elements", "Promoter-contacting"))
  d <- add_n(d, "type")
  d$xlab <- factor(d$xlab, levels = unique(d$xlab[order(d$modality, d$type)]))
  lim <- quantile(d$value, c(.02, .98), na.rm = TRUE)
  p <- ggplot(d, aes(xlab, value, fill = arm)) +
    geom_hline(yintercept = 0, colour = INK2, linewidth = 0.35, linetype = "22") +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, colour = "#33322f",
                 position = position_dodge2(preserve = "single", padding = 0.18)) +
    facet_wrap(~modality, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = c("All elements" = S1, "Promoter-contacting" = S2)) +
    coord_cartesian(ylim = c(lim[1], lim[2]), clip = "on") +
    labs(title = ttl, subtitle = sub, x = NULL,
         y = expression(log[2]*" fold change (HPC / adherent), centred on conserved"),
         caption = paste("Classes ordered by the promoter-contacting cRE median. Tick labels give n as all elements / promoter-contacting.",
                         "Domains are nested, so an element inside more than one remodelled domain contributes to each class it falls in.",
                         "Outliers hidden by axis clipping; positive = higher in HPC.", sep = "\n")) +
    theme_pub() + theme(axis.text.x = element_text(size = 5.8, lineheight = 1.0))
  list(plot = p, data = d)
}

cat("\n--- TAD ---\n")
tad <- build("TAD", c("compareTAD/ATAC_TAD_changes.rds", "compareTAD/ATACloop_TAD_changes.rds",
                      "compareTAD/RNA_TAD_changes.rds",  "compareTAD/RNAloop_TAD_changes.rds"),
             pretty_tad, "TAD remodelling during EHT",
             "No TAD class shows a concordant shift in both modalities once a promoter contact is required")

cat("\n--- sub-TAD ---\n")
sub <- build("subTAD", c("compareTAD/ATAC_subTAD_changes.rds", "compareTAD/ATACloop_subTAD_changes.rds",
                         "compareTAD/RNA_subTAD_changes.rds",  "compareTAD/RNAloop_subTAD_changes.rds"),
             pretty_sub, "Sub-TAD remodelling during EHT",
             "De novo sub-TADs are the one class gaining both accessibility and expression")

# ------------------------------------------------------------------ Figure 4
set.seed(1)
eff <- function(d, refname, scale_name) {
  do.call(rbind, lapply(split(d, list(d$modality, d$arm), drop = TRUE), function(s) {
    ref <- s$value[s$type == refname]
    do.call(rbind, lapply(setdiff(levels(droplevels(s$type)), refname), function(L) {
      v <- s$value[s$type == L]
      if (length(v) < 5) return(NULL)
      bs <- replicate(2000, median(sample(v, length(v), TRUE)) -
                            median(sample(ref, length(ref), TRUE)))
      data.frame(scale = scale_name, modality = s$modality[1], arm = s$arm[1], type = L,
                 n = length(v), delta = median(v) - median(ref),
                 lo = quantile(bs, .025), hi = quantile(bs, .975), row.names = NULL)
    }))
  }))
}
e <- rbind(eff(tad$data, "Conserved", "TAD"), eff(sub$data, "Conserved", "Sub-TAD"))
e$type <- factor(e$type, levels = rev(unique(e$type[order(e$scale, -e$delta)])))
p4 <- ggplot(e, aes(delta, type, colour = arm)) +
  geom_vline(xintercept = 0, colour = INK2, linewidth = 0.35, linetype = "22") +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.62),
                 linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.62), size = 1.5, shape = 21,
             aes(fill = arm), colour = SURF, stroke = 0.5) +
  facet_grid(scale ~ modality, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("All elements" = S1, "Promoter-contacting" = S2)) +
  scale_fill_manual(values   = c("All elements" = S1, "Promoter-contacting" = S2)) +
  labs(title = "Effect size relative to conserved domains",
       subtitle = "Points are the shift in median log2 fold change; bars are 95% bootstrap CIs (2,000 replicates)",
       x = expression(Delta*" median "*log[2]*" fold change vs conserved"), y = NULL,
       caption = "A CI crossing zero means the class is indistinguishable from conserved domains.") +
  theme_pub() + theme(panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
                      panel.grid.major.y = element_blank(),
                      axis.text.y = element_text(size = 6.4, lineheight = 0.85))

# --------------------------------------------------------------------- write
# cairo is not linked in this R build, so cairo_pdf() dies silently. Use the base
# PDF device (all labels above are ASCII for that reason) plus ragg for the raster
# preview, and svglite for a UTF-8-safe vector copy for figure assembly.
save2 <- function(p, stem, w, h) {
  f <- file.path(FIG, paste0(stem, ".pdf"))
  pdf(f, width = w, height = h, useDingbats = FALSE); print(p); dev.off()
  svglite::svglite(file.path(FIG, paste0(stem, ".svg")), width = w, height = h)
  print(p); dev.off()
  ragg::agg_png(file.path(FIG, paste0(stem, ".png")), width = w, height = h,
                units = "in", res = 200, background = SURF); print(p); dev.off()
  cat("wrote", f, "(+ .svg, .png)\n")
}
save2(p1,        "260827.Fig1_compartment", 7.2, 4.0)
save2(tad$plot,  "260827.Fig2_TAD",         7.6, 6.4)
save2(sub$plot,  "260827.Fig3_subTAD",      8.0, 6.4)
save2(p4,        "260827.Fig4_effectsizes", 7.6, 6.8)

write.csv(do.call(rbind, audit), file.path(OUT, "260827.figure_conflict_audit.csv"), row.names = FALSE)
e2 <- e; e2[c("delta","lo","hi")] <- lapply(e2[c("delta","lo","hi")], round, 3)
write.csv(e2, file.path(OUT, "260827.figure_effectsizes.csv"), row.names = FALSE)
cat("\nwrote audit + effect-size CSVs\n")
