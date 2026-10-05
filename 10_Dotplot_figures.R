# Run as a standalone script: Rscript 10_Dotplot_figures.R (or source in RStudio)
# Run as a SLURM script:      add your own #SBATCH header and submit: sbatch --wrap="Rscript 10_Dotplot_figures.R"
# Interface area dot plots — tetramer-only and protein-DNA models
# Reads summary_interface_areas.csv from:
#   1. This directory           (tetramer-only AlphaFold models)
#   2. ../protein-dna/model_dir/  (protein+DNA models)
#
# Produces two figures in the style of fig1_interface_dotplot:
#   fig_dotplot_tetramer.pdf/.png    — 6-panel (protein pairs only)
#   fig_dotplot_proteindna.pdf/.png  — 8-panel (protein pairs + AB-DNA, CD-DNA)
#
# Required packages: ggplot2, dplyr, tidyr, scales

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ── paths ─────────────────────────────────────────────────────────────────────
BASE_TET <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) {
    "/path/to/alphafold/tetramer_only"   # set to the directory with the tetramer-only results
  }
)
BASE_DNA <- file.path(
  dirname(BASE_TET),
  "protein-dna/model_dir"   # set to the protein+DNA results directory
)

# ── shared mappings ───────────────────────────────────────────────────────────

pair_labels_pp <- c(
  "A-B" = "PI–AP3",
  "A-C" = "PI–AP1",
  "A-D" = "PI–SEP3",
  "B-C" = "AP3–AP1",
  "B-D" = "AP3–SEP3",
  "C-D" = "AP1–SEP3"
)
pair_labels_dna <- c(
  pair_labels_pp,
  "AB_DNA" = "PI/AP3–DNA",
  "CD_DNA" = "AP1/SEP3–DNA"
)

# Species key → italic scientific label
# Handles both "_1"-suffixed and bare names from the two CSV files
species_labels <- c(
  "arabidopsis" = "A. thaliana",
  "vitis"       = "V. vinifera",
  "rosa"        = "R. chinensis",
  "malus_1"     = "M. domestica",
  "malus"       = "M. domestica",
  "zizyphus"    = "Z. jujuba",
  "zizphus"     = "Z. jujuba",
  "ulmus_1"     = "U. minor",
  "ulmus"       = "U. minor",
  "eleagnus_1"  = "E. mollis",
  "eleagnus"    = "E. mollis",
  "humulus_1"   = "H. lupulus",
  "humulus"     = "H. lupulus",
  "artocarpus"  = "A. heterophyllus",
  "morus_1"     = "M. alba",
  "morus"       = "M. alba",
  "urtica"      = "U. urens"
)

# Canonical order: petalous first (top of plot), then apetalous
sp_order_sci <- c(
  "A. thaliana", "V. vinifera", "R. chinensis", "M. domestica", "Z. jujuba",
  "E. mollis", "U. minor", "H. lupulus", "M. alba", "A. heterophyllus", "U. urens"
)

petalous_sci  <- c("A. thaliana", "V. vinifera", "R. chinensis",
                   "M. domestica", "Z. jujuba")
apetalous_sci <- c("E. mollis", "U. minor", "H. lupulus",
                   "M. alba", "A. heterophyllus", "U. urens")

col_pet  <- "#E05C4B"
col_apet <- "#3A7DBF"

# ── helper: tidy an areas CSV into long form ──────────────────────────────────
tidy_areas <- function(csv_path, pair_map) {
  raw <- read.csv(csv_path, check.names = FALSE)

  # Rename columns using pair_map (keeps only mapped pairs)
  raw_long <- raw %>%
    pivot_longer(-species, names_to = "pair_raw", values_to = "iface_area") %>%
    filter(pair_raw %in% names(pair_map)) %>%
    mutate(
      pair    = pair_map[pair_raw],
      sci     = species_labels[species],
      group   = case_when(
        sci %in% petalous_sci  ~ "Petalous",
        sci %in% apetalous_sci ~ "Apetalous",
        TRUE                   ~ "Unknown"
      ),
      pair    = factor(pair, levels = unname(pair_map)),
      group   = factor(group, levels = c("Petalous", "Apetalous")),
      species = factor(sci, levels = rev(sp_order_sci))
    ) %>%
    filter(!is.na(sci))   # drop rows with no label match

  raw_long
}

# ── helper: build dot plot ────────────────────────────────────────────────────
make_dotplot <- function(df, ncol = 3) {
  means_df <- df %>%
    group_by(pair, group) %>%
    summarise(mean_area = mean(iface_area, na.rm = TRUE), .groups = "drop")

  ggplot(df, aes(x = iface_area, y = species,
                 color = group, shape = group)) +
    geom_point(size = 3.2, alpha = 0.9, stroke = 0.4) +
    geom_vline(
      data = means_df,
      aes(xintercept = mean_area, color = group),
      linetype = "dashed", linewidth = 0.8, alpha = 0.75
    ) +
    scale_color_manual(
      values = c("Petalous" = col_pet, "Apetalous" = col_apet),
      name = NULL
    ) +
    scale_shape_manual(
      values = c("Petalous" = 16, "Apetalous" = 18),
      name = NULL
    ) +
    scale_x_continuous () +
    facet_wrap(~ pair, scales = "free_x", ncol = ncol) +
    labs(
      x = "Interface Area (Å²)",
      y = "Species"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.text         = element_text(face = "bold", size = 15),
      strip.background   = element_rect(fill = "gray92"),
      axis.text.y        = element_text(face = "italic", size = 13, color="black"),
      axis.text.x        = element_text(size = 12, color="black"),
      axis.title         = element_text(size = 20, color="black"),
      legend.position    = "bottom",
      legend.text        = element_text(size = 20, color="black"),
      panel.grid.major.y = element_line(linetype = "dotted", color = "gray80"),
      panel.grid.minor   = element_blank()
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# Figure A — Tetramer-only (6 protein pairs)
# ══════════════════════════════════════════════════════════════════════════════

df_tet <- tidy_areas(
  file.path(BASE_TET, "summary_interface_areas.csv"),
  pair_labels_pp
)

pA <- make_dotplot(
  df_tet,
  ncol  = 3
)
pA
ggsave(file.path(BASE_TET, "fig_dotplot_tetramer.pdf"), pA,
       width = 13, height = 8, device = pdf)
ggsave(file.path(BASE_TET, "fig_dotplot_tetramer.png"), pA,
       width = 13, height = 8, dpi = 200)
cat("Saved: fig_dotplot_tetramer.pdf / .png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Figure B — Protein-DNA (6 protein pairs + AB-DNA + CD-DNA)
# ══════════════════════════════════════════════════════════════════════════════
df_dna <- tidy_areas(
  file.path(BASE_DNA, "summary_interface_areas.csv"),
  pair_labels_dna
)

pB <- make_dotplot(
  df_dna,
  ncol  = 4
)
pB
ggsave(file.path(BASE_TET, "fig_dotplot_proteindna.pdf"), pB,
       width = 17, height = 8, device = pdf)
ggsave(file.path(BASE_TET, "fig_dotplot_proteindna.png"), pB,
       width = 17, height = 8, dpi = 200)
cat("Saved: fig_dotplot_proteindna.pdf / .png\n")

