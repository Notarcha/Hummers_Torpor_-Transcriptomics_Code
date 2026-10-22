# ==========================================================================
# Figure 2, assembled in R:  A = heatmap (left, full height)   B = within-state p-value histogram vs null (top right)
#                            C = class proportions vs null (bottom right)
# Reads only SAVED outputs of the main script (no model refitting, no new permutations):
#   <predictor>_classified_responsive_genes.csv, <predictor>_permutation_null_summary.csv, <predictor>_perm_partial.rds
# plus the raw expression/metadata (to rebuild the heatmap). Run after the main script has finished for the chosen predictor.
# Needs packages: tidyverse, ComplexHeatmap, circlize, patchwork, here
# ==========================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(ComplexHeatmap); library(circlize); library(patchwork); library(here); library(grid)
})

predictors <- c("Surf_Temp", "log_VO2")   # main figure and the supplementary (metabolic-rate) version; both are built in one run
output_dir <- "/Users/notarcha/Desktop/PNAS_Code/NewOutputs"
out <- function(f) file.path(output_dir, paste0(predictor, "_", f))

classes_shown <- c("Linear", "Early step", "Late step", "Graded / non-monotonic")    # classes drawn in the heatmap
top_n_per_class <- 100
show_gene_names <- TRUE
fig_width  <- 22          # inches ("huge canvas"); reduce if the journal needs a smaller file
fig_height <- 20          # taller because the Graded / non-monotonic block adds ~27 rows

for (predictor in predictors) {      # `out()` and the data prep below read the current `predictor`
  message("=== Building Figure 2 for: ", predictor, " ===")
  
  # ---- load data and rebuild per-tissue Z-scores (same as the main script) ---------------------------
  seq_metadata <- read.csv(here("Results", "RNASeq_metadata.csv"), stringsAsFactors = FALSE)
  rlog_counts  <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"), row.names = 1, check.names = FALSE)
  metabolic    <- read.csv(here("WGCNA_data", "metadata_surfTemp_mlO2_min_euthTime.csv"), stringsAsFactors = FALSE)
  
  sample_info <- seq_metadata %>%
    left_join(metabolic %>% select(BirdID, Surf_Temp, mean_VO2), by = "BirdID") %>%
    rename(Mean_VO2 = mean_VO2)
  common <- intersect(colnames(rlog_counts), sample_info$X)
  sample_info <- sample_info[match(common, sample_info$X), ] %>%
    mutate(
      State  = factor(Metabolic_State, levels = c("N", "T", "D"), labels = c("Normothermy", "Transition", "Deep Torpor")),
      Tissue = forcats::fct_recode(factor(Tissue), `Proximal Gut` = "Gut1", `Medial Gut` = "Gut2",
                                   `Distal Gut` = "Gut3", Pectoral = "Pect"),
      X_raw  = if (predictor == "Surf_Temp") Surf_Temp else log10(Mean_VO2)
    ) %>%
    filter(!is.na(X_raw), !is.na(State)) %>% arrange(X_raw)
  
  counts <- as.matrix(rlog_counts[, sample_info$X])
  zmat <- counts
  for (tt in levels(sample_info$Tissue)) {
    idx <- which(sample_info$Tissue == tt); x <- counts[, idx, drop = FALSE]; mu <- rowMeans(x)
    zmat[, idx] <- (x - mu) / sqrt(rowSums((x - mu)^2) / (length(idx) - 1))
  }
  zmat <- zmat[apply(zmat, 1, function(r) all(is.finite(r))), ]
  
  results <- read.csv(out("classified_responsive_genes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  perm    <- read.csv(out("permutation_null_summary.csv"), stringsAsFactors = FALSE)
  perm_rds <- out("perm_partial.rds")
  if (!file.exists(perm_rds)) stop("Missing ", perm_rds, " (written by the main script during the permutation loop).")
  perm_res <- readRDS(perm_rds)
  stopifnot(length(perm_res) == nrow(perm))     # all permutations finished
  
  # ---- Panel A: heatmap ------------------------------------------------------------------------------
  state_colors <- c("Normothermy" = "#7570B3", "Transition" = "#1B9E77", "Deep Torpor" = "#D95F02")
  top <- results %>% filter(Class_Fine %in% classes_shown) %>%
    group_by(Class_Fine) %>% slice_max(State_Var, n = top_n_per_class, with_ties = FALSE) %>% ungroup()
  si_h <- sample_info[order(sample_info$X_raw, decreasing = TRUE), ]        # warm -> cold
  mat  <- zmat[top$Gene, si_h$X]
  tissue_cols <- setNames(scales::hue_pal()(nlevels(si_h$Tissue)), levels(si_h$Tissue))
  
  ht_opt$message <- FALSE
  top_annot <- HeatmapAnnotation(
    `Metabolic Rate`  = si_h$Mean_VO2, `Surf Temp` = si_h$Surf_Temp,
    `Metabolic State` = si_h$State,    Tissue = si_h$Tissue,
    col = list(`Metabolic State` = state_colors, Tissue = tissue_cols),
    annotation_name_gp = gpar(fontsize = 12))
  ht <- Heatmap(mat, name = "Per-tissue Z",
                col = colorRamp2(c(-2, 0, 2), c("goldenrod3", "white", "purple")),
                top_annotation = top_annot,
                row_split = factor(top$Class_Fine, levels = classes_shown),
                row_title_gp = gpar(fontsize = 14), cluster_row_slices = FALSE,
                cluster_rows = TRUE, cluster_columns = FALSE,
                show_row_names = show_gene_names, row_names_side = "right",
                row_names_gp = gpar(fontsize = 7), show_column_names = FALSE, use_raster = FALSE)
  gA <- grid.grabExpr(draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm")))
  pA <- wrap_elements(full = gA)
  
  # ---- Panel B: within-state p-value histogram vs permutation null -------------------------------------
  bins <- seq(0, 1, 0.05)
  H <- do.call(rbind, lapply(perm_res, function(r) r$hist[[1]]))
  H <- H[complete.cases(H), , drop = FALSE]
  obs_hist <- hist(results$p_within, breaks = bins, plot = FALSE)$counts / nrow(results)
  band <- tibble(mid = bins[-length(bins)] + 0.025, obs = obs_hist,
                 lo = apply(H, 2, quantile, 0.025), mean = colMeans(H), hi = apply(H, 2, quantile, 0.975))
  pB <- ggplot(band, aes(mid)) +
    geom_col(aes(y = obs), width = 0.047, fill = "grey60") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "firebrick", alpha = 0.25) +
    geom_line(aes(y = mean), colour = "firebrick", linewidth = 0.9) +
    labs(x = "p: temperature effect beyond state (within-state slope)", y = "Proportion of responsive genes") +
    theme_minimal(base_size = 15)
  
  # ---- Panel C: class proportions, observed vs permutation null ----------------------------------------
  long_null <- perm %>% select(perm, Linear = frac_linear, `Non-linear` = frac_nonlin, Indistinguishable = frac_indist) %>%
    pivot_longer(-perm, names_to = "Class", values_to = "frac") %>% filter(is.finite(frac))
  long_obs <- tibble(Class = c("Linear", "Non-linear", "Indistinguishable"),
                     frac = c(mean(results$Class_Broad == "Linear"), mean(results$Class_Broad == "Non-linear"),
                              mean(results$Class_Broad == "Indistinguishable")))
  pC <- ggplot(long_null, aes(Class, frac)) +
    geom_violin(fill = "grey85", colour = "grey40") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    geom_point(data = long_obs, colour = "firebrick", size = 5, shape = 18) +
    labs(x = NULL, y = "Proportion of responsive genes") +
    theme_minimal(base_size = 15)
  
  # ---- assemble: A on the left (full height), B over C on the right; bold letters A, B, C ---------------
  fig2 <- wrap_plots(pA, pB, pC, design = "AAB\nAAC") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 28))
  
  ggsave(out("Figure2_assembled.png"),  fig2, width = fig_width, height = fig_height, dpi = 300, bg = "white")
  ggsave(out("Figure2_assembled.pdf"),  fig2, width = fig_width, height = fig_height, bg = "white")
  ggsave(out("Figure2_assembled.tiff"), fig2, width = fig_width, height = fig_height, dpi = 300, bg = "white", compression = "lzw")
  message("Saved Figure 2 (", predictor, ") to: ", output_dir)
}   # end of loop over predictors