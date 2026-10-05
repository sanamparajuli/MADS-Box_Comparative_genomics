library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

# ── Read data ─────────────────────────────────────────────────────────────────
df <- read.csv("mads_box_data.csv")

# ── Orders ────────────────────────────────────────────────────────────────────
species_order <- c(
  "Rosa chinensis", "Malus domestica", "Zizyphus jujuba", "Eleagnus mollis",
  "Ulmus minor", "Cannabis sativa", "Humulus lupulus", "Morus alba",
  "Artocarpus heterophyllus", "Urtica urens"
)
family_order <- c("Rosaceae", "Rhamnaceae", "Elaeagnaceae", "Ulmaceae",
                  "Cannabaceae", "Moraceae", "Urticaceae")

# Abbreviated species labels (e.g. "A. heterophyllus")
abbrev <- function(x) sub("^(.).* (.*)$", "\\1. \\2", x)
species_labs <- setNames(abbrev(species_order), species_order)

# ── Long format ───────────────────────────────────────────────────────────────
df_long <- df %>%
  pivot_longer(cols = c(A, B, CD, E),
               names_to = "Class", values_to = "Count") %>%
  mutate(
    Species = factor(Species, levels = species_order),
    Family  = factor(Family,  levels = family_order),
    Class   = factor(Class,   levels = c("A", "B", "CD", "E"))
  )

# ── Palette: soft, harmonious, print-friendly ─────────────────────────────────
class_cols <- c("A"  = "#6FB07F",   # sage green
                "B"  = "#5C8AA6",   # dusty blue
                "CD" = "#C97B7B",   # muted rose
                "E"  = "#E0B16A")   # warm sand
class_labs <- c("A"  = "A (AP1/FUL)",
                "B"  = "B (AP3/PI)",
                "CD" = "C/D (AG/AGL5/AGL11)",
                "E"  = "E (SEP)")

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df_long, aes(x = Species, y = Count, fill = Class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 5, colour = "black") +
  facet_grid(. ~ Family, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_discrete(labels = species_labs) +
  scale_fill_manual(values = class_cols, labels = class_labs,
                    name = "Floral MADS-box class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     breaks = seq(0, 10, 2)) +
  labs(x = NULL, y = "Number of orthologs",
       ) +
  theme_classic(base_size = 20) +
  theme(
    plot.title       = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    axis.text.x      = element_text(angle = 90, hjust = 1, face = "italic", size = 15),
    axis.text.y      = element_text(size = 15),
    axis.title.y     = element_text(face = "bold", size = 16),
    axis.line        = element_line(linewidth = 0.4),
    strip.background = element_rect(fill = "grey", colour = NA),
    strip.text.x     = element_text(face = "bold.italic", size = 10,
                                    colour = "black",
                                    margin = margin(4, 2, 4, 2)),
    strip.placement  = "outside",
    panel.spacing.x  = unit(0.15, "lines"),
    legend.title     = element_text(face = "bold", size = 13),
    legend.text      = element_text(size = 15),
    legend.position  = "right",
    plot.margin      = margin(10, 12, 6, 8)
  )
p
ggsave("mads_box_barplot_pub.png", p, width = 12.5, height = 6.5, dpi = 300)
ggsave("mads_box_barplot_pub.pdf", p, width = 12.5, height = 6.5)
message("Saved.")
