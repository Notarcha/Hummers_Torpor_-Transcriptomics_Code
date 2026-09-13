library(tidyverse)
library(here)

# ==========================================
# 0. FILE PATHS & DIRECTORY CREATION
# ==========================================
output_dir <- here("Outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==========================================
# 1. Flowchart Nodes (Linear 6-Step Pipeline)
# ==========================================
# 6 sequential steps vertically aligned and evenly spaced (Y = 6 to 1)
nodes <- tibble(
  id = 1:6,
  x = rep(2, 6),
  y = 6:1, # Equal 1.0 spacing
  label = c(
    "1. Tissue-Specific WGCNA Analysis\n(Identify Modules per Tissue)",
    "2. Multiple Linear Mixed-Effects Regression (lme4)\nModule Eigengene ~ Metabolic Rate + Surface Body Temp + Time of the Night",
    "3. Retain Trait-Correlated Modules\n Independent effects of MR & Surf Temp \n Exclude modules correlated with Time of the Night",
    "4. KEGG Pathway Enrichment Analysis\nfor each module per tissue type",
    "5. Pool All Enrichment Terms per Tissue Type\nacross correlated modules",
    "6. Pooled KEGG Enrichment Bubble Plot\nAcross Tissues"
  ),
  # Greyscale fill gradient: dark charcoal at top to light grey at bottom
  fill_color = c("#333333", "#555555", "#777777", "#999999", "#BBBBBB", "#DDDDDD"),
  text_color = c("white", "white", "white", "black", "black", "black")
)

# Vertical connecting arrows between sequential steps
main_edges <- tibble(
  x    = rep(2, 5),
  y    = c(5.62, 4.62, 3.62, 2.62, 1.62), # Arrow starts
  xend = rep(2, 5),
  yend = c(5.38, 4.38, 3.38, 2.38, 1.38)  # Arrow ends
)

# ==========================================
# 2. Render Vertical Flowchart
# ==========================================
p_flowchart <- ggplot() +
  # Vertical Connecting Arrows
  geom_segment(
    data = main_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
    linewidth = 0.8, color = "black"
  ) +
  # Flowchart Boxes
  geom_tile(
    data = nodes,
    aes(x = x, y = y, fill = I(fill_color)),
    width = 2.5,
    height = 0.68,
    color = "black", linewidth = 0.6
  ) +
  # Box Labels
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label, color = I(text_color)),
    size = 3.1, fontface = "bold", lineheight = 0.95
  ) +
  coord_cartesian(xlim = c(0.6, 3.4), ylim = c(0.4, 6.6)) +
  theme_void() +
  theme(
    plot.margin = margin(15, 15, 15, 15),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
  ) +
  labs(title = "WGCNA Module Selection & KEGG Enrichment Pipeline")

# Save High-Resolution Image to Outputs folder
ggsave(file.path(output_dir, "wgcna_kegg_pipeline_grayscale.jpeg"), p_flowchart, width = 8.5, height = 9, dpi = 300)
print(p_flowchart)