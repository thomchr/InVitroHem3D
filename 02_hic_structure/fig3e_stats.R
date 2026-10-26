# Fig. 3e — exact statistics for the editor's reporting checklist (KhanhQs.docx, 2026-10-26)
# Inputs are read-only; outputs go to OUT_DIR.
source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
suppressMessages(library(data.table))
d   <- file.path(DATA_DIR, "hiC/ABcompartment")
out <- "<OUT_DIR>/"

gene <- readRDS(file.path(d, "gene_compartment_df.rds"))
reg  <- readRDS(file.path(d, "region_compartment_df.rds"))
sets <- list(RNA  = data.frame(g = gene$Compartment_Change, y = gene$Log2FC),
             ATAC = data.frame(g = reg$RegionType,          y = reg$Log2FC))

# the 8 brackets drawn in the PDFs (read off the figure, bottom to top)
pairs <- list(c("AtoB","AtoB_looped"), c("AtoB","BtoA"), c("AtoB","Unchange"),
              c("AtoB_looped","BtoA_looped"), c("AtoB_looped","Unchange"),
              c("BtoA","BtoA_looped"), c("BtoA","Unchange"), c("BtoA_looped","Unchange"))
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "NS.")))

summ <- list(); pw <- list(); glob <- list()
for (m in names(sets)) {
  x <- sets[[m]]
  cat("\n==", m, ": non-finite y =", sum(!is.finite(x$y)), "\n")
  x <- x[is.finite(x$y), ]
  for (g in levels(x$g)) {
    v <- x$y[x$g == g]; b <- boxplot.stats(v)
    q <- quantile(v, c(.25, .5, .75), type = 7); iqr <- q[3] - q[1]
    summ[[length(summ) + 1]] <- data.frame(modality = m, group = g, n = length(v),
      median = q[2], Q1 = q[1], Q3 = q[3], whisker_lo = b$stats[1], whisker_hi = b$stats[5],
      notch_lo = q[2] - 1.58 * iqr / sqrt(length(v)), notch_hi = q[2] + 1.58 * iqr / sqrt(length(v)),
      n_beyond_whiskers = length(b$out), min = min(v), max = max(v), mean = mean(v))
  }
  a  <- summary(aov(y ~ g, data = x))[[1]]
  kw <- kruskal.test(y ~ g, data = x)
  glob[[m]] <- data.frame(modality = m, N = nrow(x),
    anova_F = a$`F value`[1], df1 = a$Df[1], df2 = a$Df[2], anova_p = a$`Pr(>F)`[1],
    eta2 = a$`Sum Sq`[1] / sum(a$`Sum Sq`),
    kw_H = unname(kw$statistic), kw_df = unname(kw$parameter), kw_p = kw$p.value)
  for (pr in pairs) {
    a1 <- x$y[x$g == pr[1]]; a2 <- x$y[x$g == pr[2]]
    w  <- wilcox.test(a1, a2, conf.int = TRUE)      # two-sided; estimate = Hodges-Lehmann shift
    t  <- t.test(a1, a2)                            # two-sided Welch
    pw[[length(pw) + 1]] <- data.frame(modality = m, group1 = pr[1], group2 = pr[2],
      n1 = length(a1), n2 = length(a2), delta_median = median(a1) - median(a2),
      W = unname(w$statistic), HL_shift = unname(w$estimate),
      HL_lo = w$conf.int[1], HL_hi = w$conf.int[2], wilcox_p = w$p.value,
      welch_t = unname(t$statistic), welch_df = unname(t$parameter), welch_p = t$p.value)
  }
}
summ <- do.call(rbind, summ); glob <- do.call(rbind, glob); pw <- do.call(rbind, pw)
# rank-biserial r = 2*W/(n1*n2) - 1  (positive = group1 tends higher)
pw$r_rankbiserial <- 2 * pw$W / (as.numeric(pw$n1) * pw$n2) - 1
# log10 P from the normal approximation wilcox.test itself uses (no tie correction; for P that underflow to 0)
mu <- as.numeric(pw$n1) * pw$n2 / 2; sg <- sqrt(as.numeric(pw$n1) * pw$n2 * (pw$n1 + pw$n2 + 1) / 12)
pw$wilcox_log10p_approx <- (log(2) + pnorm(-abs((pw$W - mu) / sg), log.p = TRUE)) / log(10)
glob$anova_log10p <- pf(glob$anova_F, glob$df1, glob$df2, lower.tail = FALSE, log.p = TRUE) / log(10)
glob$kw_log10p    <- pchisq(glob$kw_H, glob$kw_df, lower.tail = FALSE, log.p = TRUE) / log(10)
# gene groups overlap (a gene can be both located in and loop-contacted within a switching region)
gl <- split(gene$Gene_Name[is.finite(gene$Log2FC)], gene$Compartment_Change[is.finite(gene$Log2FC)])
rl <- split(reg$RegionName, reg$RegionType)
pw$n_shared <- mapply(function(m, a, b) length(intersect(if (m == "RNA") gl[[a]] else rl[[a]],
                                                         if (m == "RNA") gl[[b]] else rl[[b]])),
                      pw$modality, pw$group1, pw$group2)
pw <- do.call(rbind, lapply(split(pw, pw$modality), function(z) {
  z$wilcox_p_holm <- p.adjust(z$wilcox_p, "holm"); z$wilcox_p_bh <- p.adjust(z$wilcox_p, "BH")
  z$stars_wilcox <- stars(z$wilcox_p); z$stars_welch <- stars(z$welch_p); z }))

options(width = 250, digits = 4)
print(summ, row.names = FALSE); print(glob, row.names = FALSE)
print(pw[, c("modality","group1","group2","n1","n2","delta_median","HL_shift","HL_lo","HL_hi",
             "r_rankbiserial","n_shared","wilcox_p","wilcox_log10p_approx","wilcox_p_holm","stars_wilcox","welch_p","stars_welch")], row.names = FALSE)

fwrite(summ, file.path(out, "261026.Fig3e_box_summary.csv"))
fwrite(glob, file.path(out, "261026.Fig3e_global_tests.csv"))
fwrite(pw,   file.path(out, "261026.Fig3e_pairwise_tests.csv"))
