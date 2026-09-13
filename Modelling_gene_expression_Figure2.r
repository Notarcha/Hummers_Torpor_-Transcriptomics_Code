library(tidyverse)
library(lme4)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(ggpubr)
library(here)

# ==========================================
# 1. FILE PATHS & DATA LOADING
# ==========================================
output_dir <- here("Outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

path_meta_seq  <- here("Results", "RNASeq_metadata.csv")
path_rlog      <- here("Results", "rlog_normalized_counts_STAR.csv")
path_metabolic <- here("WGCNA_data", "metadata_surfTemp_mlO2_min_euthTime.csv")

seq_metadata <- read.csv(path_meta_seq, stringsAsFactors = FALSE)
rlog_counts  <- read.csv(path_rlog, row.names = 1, check.names = FALSE)
metabolic    <- read.csv(path_metabolic, stringsAsFactors = FALSE)

filtered_counts <- rlog_counts[rowSums(rlog_counts, na.rm = TRUE) > 10, ]

# ==========================================
# 2. METADATA ALIGNMENT & DATA TIDYING
# ==========================================
sample_info <- seq_metadata %>%
  left_join(metabolic %>% select(BirdID, Surf_Temp, mean_VO2), by = "BirdID") %>%
  rename(Mean_VO2 = mean_VO2)

common_samples  <- intersect(colnames(filtered_counts), sample_info$X)
filtered_counts <- filtered_counts[, common_samples]
sample_info     <- sample_info[match(common_samples, sample_info$X), ]

sample_info <- sample_info %>%
  arrange(Surf_Temp) %>%
  mutate(Metabolic_State_Label = case_when(
    Metabolic_State == "N" ~ "Normothermy",
    Metabolic_State == "D" ~ "Deep Torpor",
    Metabolic_State == "T" ~ "Transition",
    TRUE ~ Metabolic_State
  )) %>%
  mutate(
    Metabolic_State_Label = factor(Metabolic_State_Label, levels = c("Normothermy", "Transition", "Deep Torpor")),
    Tissue = as.factor(Tissue)
  )

filtered_counts <- filtered_counts[, sample_info$X]

df_long <- filtered_counts %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "X", values_to = "Expression") %>%
  inner_join(sample_info, by = "X")

df_long <- df_long %>%
  group_by(Gene) %>%
  mutate(Z_Score = (Expression - mean(Expression, na.rm = TRUE)) / sd(Expression, na.rm = TRUE)) %>%
  ungroup()

# ==========================================
# 3. FAST SCREENING FOR TEMPERATURE-RESPONSIVE GENES
# ==========================================
ht_opt$message = FALSE

responsive_genes <- df_long %>%
  group_by(Gene) %>%
  summarise(
    pval = tryCatch(
      summary(aov(Expression ~ Metabolic_State_Label, data = cur_data()))[[1]][["Pr(>F)"]][1],
      error = function(e) 1
    ),
    .groups = "drop"
  ) %>%
  filter(pval < 0.01) %>%
  pull(Gene)

df_responsive <- df_long %>% filter(Gene %in% responsive_genes)

# ==========================================
# 4. AIC MODEL EVALUATION & PRECISE POOL CLASSIFICATION
# ==========================================
responsive_gene_list <- unique(df_responsive$Gene)

model_fit_results <- map_dfr(responsive_gene_list, function(g) {
  gene_df <- filter(df_responsive, Gene == g)
  
  fit_linear <- tryCatch(
    suppressWarnings(lmer(Z_Score ~ Surf_Temp + (1 | Tissue), data = gene_df, REML = FALSE)),
    error = function(e) NULL
  )
  
  fit_state <- tryCatch(
    suppressWarnings(lmer(Z_Score ~ Metabolic_State_Label + (1 | Tissue), data = gene_df, REML = FALSE)),
    error = function(e) NULL
  )
  
  if (!is.null(fit_linear) && !is.null(fit_state)) {
    aic_lin <- AIC(fit_linear)
    aic_st  <- AIC(fit_state)
    tibble(Gene = g, AIC_Linear = aic_lin, AIC_State = aic_st, Delta_AIC = aic_lin - aic_st)
  } else {
    NULL
  }
})

gene_state_means <- df_responsive %>%
  group_by(Gene, Metabolic_State_Label) %>%
  summarise(Mean_Z = mean(Z_Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Metabolic_State_Label, values_from = Mean_Z) %>%
  mutate(
    Early_Change = abs(Transition - Normothermy),
    Late_Change  = abs(`Deep Torpor` - Transition),
    State_Variance = apply(select(., Normothermy, Transition, `Deep Torpor`), 1, var),
    Direction = ifelse(`Deep Torpor` > Normothermy, "Direction: Up in Torpor", "Direction: Down in Torpor")
  )

df_model_fits <- model_fit_results %>%
  inner_join(gene_state_means, by = "Gene") %>%
  mutate(Pool_Clean = case_when(
    Delta_AIC < -2 ~ "Continuous Linear",
    Delta_AIC >= -2 & Early_Change >= Late_Change ~ "Early Responders",
    TRUE ~ "Late Responders"
  )) %>%
  mutate(Pool_Clean = factor(Pool_Clean, levels = c(
    "Continuous Linear",
    "Early Responders",
    "Late Responders"
  )))

cat("\n--- Breakdown of Responsive Pools by Direction ---\n")
print(table(df_model_fits$Pool_Clean, df_model_fits$Direction))

df_responsive_plot <- df_responsive %>%
  inner_join(df_model_fits %>% select(Gene, Delta_AIC, Pool_Clean, Direction, State_Variance), by = "Gene")

state_colors <- c(
  "Normothermy" = "#7570B3",
  "Transition"  = "#1B9E77",
  "Deep Torpor" = "#D95F02"
)

# ==========================================
# 5. SAVE SEPARATE CSV FILES
# ==========================================
write.csv(df_model_fits %>% filter(Pool_Clean == "Continuous Linear"), file.path(output_dir, "continuous_linear_responders_genes.csv"), row.names = FALSE)
write.csv(df_model_fits %>% filter(Pool_Clean == "Early Responders"), file.path(output_dir, "early_responders_genes.csv"), row.names = FALSE)
write.csv(df_model_fits %>% filter(Pool_Clean == "Late Responders"), file.path(output_dir, "late_responders_genes.csv"), row.names = FALSE)

# ==========================================
# 6. SELECT TOP 100 UNIQUE GENES PER POOL
# ==========================================
top100_gene_ids <- df_model_fits %>%
  group_by(Pool_Clean) %>%
  slice_max(order_by = State_Variance, n = 100, with_ties = FALSE) %>%
  pull(Gene)

df_top100_plot <- df_responsive_plot %>% filter(Gene %in% top100_gene_ids)

# ==========================================
# 7. FIGURE 4: VIOLINS WITH SPACED STATS
# ==========================================
tukey_results <- df_top100_plot %>%
  group_by(Pool_Clean, Direction) %>%
  tukey_hsd(Z_Score ~ Metabolic_State_Label) %>%
  filter(group1 == "Normothermy" & group2 == "Transition" | 
           group1 == "Transition" & group2 == "Deep Torpor") %>%
  add_significance("p.adj") %>%
  add_y_position(step.increase = 0.12)

p_violins_stats <- ggplot(df_top100_plot, aes(x = Metabolic_State_Label, y = Z_Score)) +
  geom_violin(aes(fill = Metabolic_State_Label), trim = FALSE, alpha = 0.5, color = "black", linewidth = 0.4) +
  geom_boxplot(aes(fill = Metabolic_State_Label), width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.9, color = "black") +
  stat_pvalue_manual(
    tukey_results, 
    label = "p.adj.signif", 
    hide.ns = FALSE,
    tip.length = 0.02,
    bracket.size = 0.5,
    size = 4.5,
    vjust = -0.4
  ) +
  facet_grid(Direction ~ Pool_Clean) +
  scale_fill_manual(values = state_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top 100 Dynamic Genes per Category",
    subtitle = "Z-score distributions with Post-Hoc Tukey HSD markers (* p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001)",
    x = "Metabolic State",
    y = "Per-Tissue Gene Z-Score",
    fill = "Metabolic State"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

ggsave(file.path(output_dir, "top100_genes_refined_violins_stats.png"), p_violins_stats, width = 12, height = 8, dpi = 300)

# ==========================================
# 8. FIGURE 5: COMPLEXHEATMAP (MATCHING TARGET PLOT)
# ==========================================

# 1. Align metadata and expression matrix
heatmap_gene_meta <- df_model_fits %>% 
  filter(Gene %in% top100_gene_ids) %>%
  distinct(Gene, .keep_all = TRUE)

scaled_matrix <- df_top100_plot %>%
  select(Gene, X, Z_Score) %>%
  pivot_wider(names_from = X, values_from = Z_Score, values_fn = mean) %>%
  column_to_rownames("Gene")

scaled_matrix_csv <- as.matrix(scaled_matrix[heatmap_gene_meta$Gene, sample_info$X])

# 2. Expression Heatmap Color Palette (Goldenrod2 -> White -> Purple)
col_expression <- colorRamp2(c(-2, 0, 2), c("goldenrod3", "#FFFFFF", "purple"))

# 3. Top Annotation Color Palettes
vo2_min <- min(sample_info$Mean_VO2, na.rm = TRUE)
vo2_max <- max(sample_info$Mean_VO2, na.rm = TRUE)
col_vo2  <- colorRamp2(c(vo2_min, (vo2_min + vo2_max)/2, vo2_max), c("#EFF3FF", "#6BAED6", "#08519C"))

temp_min <- min(sample_info$Surf_Temp, na.rm = TRUE)
temp_max <- max(sample_info$Surf_Temp, na.rm = TRUE)
col_temp <- colorRamp2(c(temp_min, (temp_min + temp_max)/2, temp_max), c("#2B83BA", "#FFFFBF", "#D7191C"))

col_state <- c(
  "Normothermy" = "#7570B3",
  "Transition"  = "#1B9E77",
  "Deep Torpor" = "#D95F02"
)

# 4. Construct Top Annotation Track
top_annot <- HeatmapAnnotation(
  `Metabolic Rate` = sample_info$Mean_VO2,
  `Surf Temp`      = sample_info$Surf_Temp,
  `Metabolic State` = sample_info$Metabolic_State_Label,
  col = list(
    `Metabolic Rate`  = col_vo2,
    `Surf Temp`       = col_temp,
    `Metabolic State` = col_state
  ),
  annotation_legend_param = list(
    `Metabolic Rate`  = list(title = "Mean VO2"),
    `Surf Temp`       = list(title = "Surf Temp (°C)"),
    `Metabolic State` = list(title = "Metabolic State")
  ),
  annotation_name_side = "right",
  annotation_name_gp   = gpar(fontsize = 12)
)

# 5. Render ComplexHeatmap
ht <- Heatmap(
  scaled_matrix_csv,
  name               = "Per-Tissue Z-Score",
  col                = col_expression,
  top_annotation     = top_annot,
  
  row_split          = factor(heatmap_gene_meta$Pool_Clean, levels = c("Continuous Linear", "Early Responders", "Late Responders")),
  row_title_gp       = gpar(fontface = "bold", fontsize = 12),
  row_title_rot      = 0,
  cluster_row_slices = FALSE,
  
  show_row_names     = TRUE,
  row_names_gp       = gpar(fontsize = 4),
  show_column_names  = FALSE,
  cluster_rows       = TRUE,
  cluster_columns    = FALSE,
  use_raster         = FALSE,
  
  heatmap_legend_param = list(
    title = "Per-Tissue Z-Score",
    at = c(-2, -1, 0, 1, 2)
  )
)

png(file.path(output_dir, "hummingbird_temp_responsive_genes.png"), width = 12, height = 16, units = "in", res = 300)
draw(
  ht, 
  heatmap_legend_side = "right", 
  annotation_legend_side = "right",
  merge_legend = TRUE
)
dev.off()

# ==========================================
# 9. FIGURE 6: UNIFIED DENSITY PLOT
# ==========================================
n_linear <- sum(df_model_fits$Delta_AIC < -2, na.rm = TRUE)
n_equiv  <- sum(df_model_fits$Delta_AIC >= -2 & df_model_fits$Delta_AIC <= 2, na.rm = TRUE)
n_switch <- sum(df_model_fits$Delta_AIC > 2, na.rm = TRUE)
n_total  <- nrow(df_model_fits)

pct_linear <- round((n_linear / n_total) * 100, 1)
pct_equiv  <- round((n_equiv / n_total) * 100, 1)
pct_switch <- round((n_switch / n_total) * 100, 1)

dens <- density(df_model_fits$Delta_AIC, na.rm = TRUE)
max_y  <- max(dens$y)
text_y <- max_y * 1.15

p_aic_density <- ggplot(df_model_fits, aes(x = Delta_AIC)) +
  annotate("rect", xmin = -2, xmax = 2, ymin = 0, ymax = Inf, fill = "grey85", alpha = 0.35) +
  geom_density(fill = "#5ab4ac", alpha = 0.55, color = "#01665e", linewidth = 0.8) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  annotate("text", x = -8, y = text_y, 
           label = paste0("Favors Continuous\nLinear Model\n(ΔAIC < -2)\nN = ", n_linear, " (", pct_linear, "%)"), 
           color = "firebrick", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = 10, y = text_y, 
           label = paste0("Favors Discrete\nSwitch Model\n(ΔAIC > +2)\nN = ", n_switch, " (", pct_switch, "%)"), 
           color = "#01665e", fontface = "bold", size = 3.5, hjust = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Model Selection Dynamics Across Temperature-Responsive Genes",
    subtitle = paste0("Delta AIC distribution comparing Continuous Linear vs. Discrete Switch models (Total N = ", n_total, ")"),
    x = "Delta AIC (AIC_linear - AIC_state)",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "delta_aic_distribution_unified.png"), p_aic_density, width = 10, height = 6, dpi = 300)

df_continuous <- read.csv(file.path(output_dir, "continuous_linear_responders_genes.csv"))
df_early      <- read.csv(file.path(output_dir, "early_responders_genes.csv"))
df_late       <- read.csv(file.path(output_dir, "late_responders_genes.csv"))

cat("--- Gene Counts Per Response Pool ---\n")
cat("Continuous Linear Responders :", nrow(df_continuous), "\n")
cat("Early Responders             :", nrow(df_early), "\n")
cat("Late Responders              :", nrow(df_late), "\n")
cat("Total Responsive Genes       :", nrow(df_continuous) + nrow(df_early) + nrow(df_late), "\n")