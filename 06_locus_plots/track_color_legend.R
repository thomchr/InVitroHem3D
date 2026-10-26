#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Cell-type colour legend for pyGenomeTracks panels
#
# Answers Reviewer #1 point 10 ("each ATAC track uses two colors, but does not
# specify what these colors represent") and the open author comment
# "Need to specify colors in legend".
#
# Colours are taken verbatim from
#   06_locus_plots/color_tracks.txt
# and are NOT altered -- they must match the tracks already plotted in
# Figs 3g-h, 4a/b/d/f/g, 6c-f and the corresponding Supplementary panels.
#
# Design constraints, given that the palette is fixed:
#   * every swatch carries its label directly adjacent, so identity never
#     depends on colour alone (several pairs are not CVD-separable -- see
#     the audit note accompanying this script);
#   * swatches get a dark stroke, because #D3D3D3 (Stroma), #EED244 (HPC) and
#     #93E2B8 (LMPP myeloid) have <1.6:1 contrast against white and otherwise
#     disappear on the page;
#   * the adherent / non-adherent split is drawn explicitly, matching how the
#     tracks are stacked in the figures and the dashed divider in Fig. 2a.
#
# NOTE ON NOMENCLATURE: the track file uses "EMP"; the manuscript figure
# (megakaryocyte-erythroid progenitors)". EMP and MEP are different entities in
# the literature. This script keeps the source file's "EMP" and flags the
# discrepancy rather than silently renaming a cell population -- resolve before
#
#   figures/track_color_legend_horizontal.{pdf,svg,png}
#   figures/track_color_legend_vertical.{pdf,svg,png}
# ---------------------------------------------------------------------------

source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
suppressMessages({
  library(ggplot2)
  library(svglite)
})

OUT <- OUT_DIR
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

use_mep <- FALSE   # see NOTE ON NOMENCLATURE above

# --- palette, verbatim from color_tracks.txt -------------------------------
pal <- c(
  "1_Stromal"       = "#D3D3D3",
  "2_Endo"          = "#E48658",
  "3_HE"            = "#B15944",
  "4_HPC"           = "#EED244",
  "5_LMPP_Myeloid"  = "#93E2B8",
  "6_LMPP_Lymphoid" = "#6FC6A1",
  "7_EMP"           = "#F587E2",
  "8_EMP_Myeloid"   = "#B972C4",
  "9_EMP_Mk"        = "#6A3978",
  "10_EMP_Ery"      = "#B73186"
)

stem <- if (use_mep) "MEP" else "EMP"

# display labels: manuscript vocabulary, no numeric prefixes or underscores
lab_long <- c("Stroma", "Endo", "HE", "HPC",
              "LMPP (myeloid)", "LMPP (lymphoid)",
              stem, paste0(stem, " (myeloid)"),
              paste0(stem, " (Mk)"), paste0(stem, " (Ery)"))
lab_short <- c("Stroma", "Endo", "HE", "HPC",
               "LMPP-My", "LMPP-Ly",
               stem, paste0(stem, "-My"), paste0(stem, "-Mk"), paste0(stem, "-Ery"))

d <- data.frame(
  key    = names(pal),
  hex    = unname(pal),
  long   = lab_long,
  short  = lab_short,
  group  = c(rep("Adherent", 3), rep("Non-adherent", 7)),
  stringsAsFactors = FALSE
)

INK    <- "grey15"   # label ink -- never the series colour
MUTED  <- "grey40"   # group headers
STROKE <- "grey25"   # swatch outline
PT     <- 0.35146    # 1 pt in mm, for size = conversions

# --- horizontal strip: one row, divider after HE ---------------------------
dh <- transform(d, x = seq_len(nrow(d)) - 1)

p_h <- ggplot(dh) +
  geom_rect(aes(xmin = x + 0.05, xmax = x + 0.28, ymin = 0.39, ymax = 0.61),
            fill = dh$hex, colour = STROKE, linewidth = 0.3) +
  geom_text(aes(x = x + 0.33, y = 0.50, label = short),
            hjust = 0, vjust = 0.5, size = 6.5 * PT, colour = INK) +
  # adherent / non-adherent divider, mirroring Fig. 2a
  annotate("segment", x = 3.0, xend = 3.0, y = 0.30, yend = 0.70,
           linetype = "dashed", linewidth = 0.3, colour = MUTED) +
  annotate("text", x = 1.5, y = 0.80, label = "Adherent (D9+7)",
           size = 6 * PT, colour = MUTED, fontface = "italic") +
  annotate("text", x = 6.5, y = 0.80, label = "Non-adherent (D9+8)",
           size = 6 * PT, colour = MUTED, fontface = "italic") +
  coord_cartesian(xlim = c(0, 10.15), ylim = c(0.20, 0.92), expand = FALSE) +
  theme_void()

# --- vertical block: grouped, for the side of a multi-panel figure ---------
gap <- 1.5   # clearance so the group header reads as a header, not a row
dv <- d
dv$y <- c(seq_len(3), seq_len(7) + 3 + gap)          # item rows
hdr <- data.frame(y = c(0, 3 + gap), label = c("Adherent (D9+7)", "Non-adherent (D9+8)"))

p_v <- ggplot(dv) +
  geom_rect(aes(xmin = 0.00, xmax = 0.55, ymin = y - 0.32, ymax = y + 0.32),
            fill = dv$hex, colour = STROKE, linewidth = 0.3) +
  geom_text(aes(x = 0.75, y = y, label = long),
            hjust = 0, vjust = 0.5, size = 6.5 * PT, colour = INK) +
  geom_text(data = hdr, aes(x = 0, y = y, label = label),
            hjust = 0, vjust = 0.5, size = 6 * PT, colour = MUTED,
            fontface = "italic") +
  scale_y_reverse() +
  coord_cartesian(xlim = c(-0.1, 5.2), ylim = c(12.0, -0.7), expand = FALSE) +
  theme_void()

# --- write ------------------------------------------------------------------
save_all <- function(p, stem, w, h) {
  # base pdf() device: pure vector, no cairo/X11 dependency. useDingbats=FALSE
  # keeps glyphs as text for Illustrator.
  ggsave(file.path(OUT, paste0(stem, ".pdf")), p, width = w, height = h,
         units = "mm", device = grDevices::pdf, useDingbats = FALSE)
  ggsave(file.path(OUT, paste0(stem, ".svg")), p, width = w, height = h,
         units = "mm", device = svglite::svglite)
  ggsave(file.path(OUT, paste0(stem, ".png")), p, width = w, height = h,
         units = "mm", dpi = 600, bg = "white")
  cat(sprintf("  %-34s %3.0f x %3.0f mm\n", paste0(stem, ".{pdf,svg,png}"), w, h))
}

cat("Writing to ", OUT, "\n", sep = "")
save_all(p_h, "track_color_legend_horizontal", 183, 13)
save_all(p_v, "track_color_legend_vertical",    46, 57)
cat("done\n")
