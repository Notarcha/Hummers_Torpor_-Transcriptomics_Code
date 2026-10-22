# ==========================================================================
# ONE COMPOSITE FIGURE PER PATHWAY (11 figures), built only from the CSVs the pathway scripts already wrote.
#   A = gene-level heatmap (limma moderated t; per tissue + pooled)
#   B = SELF-CONTAINED set heatmap (Welch / lme; '+' = mroast 'mixed')
#   C = COMPETITIVE set heatmap (camera; pooled column = camera on whitened data)
# Nothing is refitted. Run time: seconds per pathway.
# A failure in one pathway (e.g. a missing file) is reported at the end and does not stop the others.
# Needs: tidyverse, pheatmap, patchwork
# ==========================================================================
suppressPackageStartupMessages({ library(tidyverse); library(pheatmap); library(patchwork); library(grid) })
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate; rename <- dplyr::rename
# NOTE: pheatmap is always called as pheatmap::pheatmap(), because ComplexHeatmap (if loaded in the same R session)
# masks pheatmap() with a version that returns an S4 object and ignores silent = TRUE.

# ---- WHICH PATHWAYS / WHERE ----------------------------------------------------------------------
pathways     <- c("TCA_CYCLE", "OXPHOS", "PPAR", "THERMOGENESIS", "CIRCADIAN", "SPLICEOSOME",
                  "TRANSLATION", "TRANSCRIPTION", "CELL_CYCLE", "PI3K_AKT", "FOXO")
main_fig_dir <- "/Users/notarcha/Desktop/PNAS_Code/NewOutputs"
save_png <- TRUE; save_pdf <- TRUE; save_tiff <- TRUE

# ---- settings (same values as the pathway scripts) -----------------------------------------------
tissue_order <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Lungs", "Heart", "Liver", "Pectoral")
pooled_label_gene <- "All tissues\n(pooled)"
fdr_cut    <- 0.05
set_g_lim  <- 2
gene_t_lim <- 5
my_colors  <- colorRampPalette(c("#b35806", "#f1a340", "#fee0b6", "#f7f7f7", "#d8daeb", "#998ec3", "#542788"))(100)

# layout
fig_width    <- 28           # inches
panel_widths <- c(1, 1.15)   # left (A) : right (B over C)
letter_size  <- 14           # mm (~40 pt)
label_cm     <- 2.4          # height of the letter row above each heatmap
set_row_cm   <- 1.5          # height per set row in B and C
col_label_cm <- 7            # space for the rotated column labels under B and C
# text sizes (points)
fs_gene <- list(number = 10, base = 14, row = 11, col = 14)    # panel A
fs_set  <- list(number = 15, base = 16, row = 17, col = 16)    # panels B and C

# ---- helpers -----------------------------------------------------------------------------------------
find_file <- function(name) {                       # case-insensitive (the folder mixes e.g. "Transcription_" and "TRANSCRIPTION_")
  fl  <- list.files(main_fig_dir)
  hit <- fl[tolower(fl) == tolower(name)]
  if (length(hit)) file.path(main_fig_dir, hit[1]) else NA_character_
}
rd <- function(prefix, suffix, required = TRUE) {
  f <- find_file(paste0(prefix, suffix))
  if (is.na(f)) { if (required) stop("Missing file: ", paste0(prefix, suffix)) else return(NULL) }
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
to_mat <- function(tbl, row_col, col_col, val_col, rows, cols, fill = NA_real_) {
  m <- matrix(fill, length(rows), length(cols), dimnames = list(rows, cols))
  m[cbind(match(tbl[[row_col]], rows), match(tbl[[col_col]], cols))] <- tbl[[val_col]]
  m
}
star_fun  <- function(q) ifelse(is.na(q), "", ifelse(q < 0.01, "**", ifelse(q < fdr_cut, "*", "")))
clean_lab <- function(x) trimws(sub("\\s*\\((non-)?KEGG\\)\\s*$", "", x))      # drop '(KEGG)' / '(non-KEGG)'
is_catchall <- function(x) grepl("^KEGG only", x)

label_row <- function(l) {
  ggplot() + annotate("text", x = 0, y = 0.5, label = l, hjust = 0, vjust = 0.5, size = letter_size, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 25))
}
labelled <- function(l, p) label_row(l) / p + plot_layout(heights = unit(c(label_cm, 1), c("cm", "null")))

# ==========================================================================
make_composite <- function(FILE_PREFIX) {
  gene_tbl        <- rd(FILE_PREFIX, "_gene_level_FDR.csv")
  gene_pooled_tbl <- rd(FILE_PREFIX, "_gene_level_pooled_FDR.csv")
  set_tbl         <- rd(FILE_PREFIX, "_set_level_per_tissue_FDR.csv")
  cam_pooled      <- rd(FILE_PREFIX, "_set_level_pooled_COMPETITIVE_camera.csv")
  pooled_sc       <- rd(FILE_PREFIX, "_set_level_pooled_SELFCONTAINED_lme.csv")
  membership      <- rd(FILE_PREFIX, "_gene_membership_audit.csv", required = FALSE)
  
  # set order = definition order (pooled_sc is written in that order)
  sets_ord <- pooled_sc$Set
  
  # gene -> set (as written by the pathway script) and set level order
  gene_set_map <- gene_tbl %>% distinct(Gene, Set)
  set_genes <- if (!is.null(membership)) unique(membership$Gene[membership$in_data]) else unique(gene_pooled_tbl$Gene)
  gene_set_map <- gene_set_map %>% arrange(match(Gene, set_genes))
  lev <- unique(gene_set_map$Set)
  if (is.null(membership)) lev <- c(intersect(sets_ord, lev), setdiff(lev, sets_ord))
  gene_levels <- c(lev[!is_catchall(lev)], lev[is_catchall(lev)])
  
  # ---------------- PANEL A: gene-level heatmap -------------------------------------------------
  gene_order <- gene_tbl %>% group_by(Gene) %>% summarise(mean_shift = mean(shift_Z), .groups = "drop") %>%
    left_join(distinct(gene_set_map, Gene, Set), by = "Gene") %>%
    mutate(Set = factor(Set, levels = gene_levels)) %>% arrange(Set, mean_shift)
  
  gcols <- c(tissue_order, "All tissues")
  gt  <- to_mat(gene_tbl, "Gene", "Tissue", "t",       gene_order$Gene, gcols)
  gsh <- to_mat(gene_tbl, "Gene", "Tissue", "shift_Z", gene_order$Gene, gcols)
  gst <- to_mat(mutate(gene_tbl, star = ifelse(q < fdr_cut, "*", "")), "Gene", "Tissue", "star",
                gene_order$Gene, gcols, fill = "")
  pi_ <- match(gene_order$Gene, gene_pooled_tbl$Gene)
  gt[, "All tissues"]  <- gene_pooled_tbl$t_pooled[pi_]
  gsh[, "All tissues"] <- gene_pooled_tbl$mean_shift_Z[pi_]
  gst[, "All tissues"] <- star_fun(gene_pooled_tbl$q_pooled[pi_])
  glab <- matrix("", nrow(gsh), ncol(gsh), dimnames = dimnames(gsh))
  ok <- !is.na(gsh); glab[ok] <- paste0(sprintf("%+.1f", gsh[ok]), gst[ok])
  gt_c <- pmax(pmin(gt, gene_t_lim), -gene_t_lim)
  
  ann_row <- data.frame(Set = factor(clean_lab(as.character(gene_order$Set)), levels = unique(clean_lab(gene_levels))),
                        row.names = gene_order$Gene)
  ann_col <- list(Set = setNames(scales::hue_pal()(length(unique(clean_lab(gene_levels)))), unique(clean_lab(gene_levels))))
  gene_gaps <- cumsum(table(gene_order$Set)); gene_gaps <- unique(gene_gaps[gene_gaps > 0 & gene_gaps < nrow(gt)])
  
  hA <- pheatmap::pheatmap(gt_c, display_numbers = glab, fontsize_number = fs_gene$number, number_color = "black",
                           cluster_rows = FALSE, cluster_cols = FALSE, gaps_col = length(tissue_order),
                           annotation_row = ann_row, annotation_colors = ann_col, gaps_row = gene_gaps,
                           labels_col = c("Proximal\nGut", "Medial\nGut", "Distal\nGut", "Lungs", "Heart", "Liver", "Pectoral", pooled_label_gene),
                           angle_col = "0", color = my_colors, breaks = seq(-gene_t_lim, gene_t_lim, length.out = 101),
                           na_col = "grey85", border_color = NA, fontsize = fs_gene$base, fontsize_row = fs_gene$row,
                           fontsize_col = fs_gene$col, silent = TRUE)
  
  # ---------------- PANELS B and C: set-level heatmaps ------------------------------------------
  col_names <- c(tissue_order, "All tissues")
  g_mat     <- to_mat(set_tbl, "Set", "Tissue", "g",           sets_ord, col_names)
  shift_mat <- to_mat(set_tbl, "Set", "Tissue", "shift_set_Z", sets_ord, col_names)
  g_mat[pooled_sc$Set, "All tissues"]     <- pooled_sc$g_pooled
  shift_mat[pooled_sc$Set, "All tissues"] <- pooled_sc$pooled_shift_Z
  
  set_heatmap <- function(star_mat, pooled_label) {
    lab <- matrix("", nrow(shift_mat), ncol(shift_mat), dimnames = dimnames(shift_mat))
    okk <- !is.na(shift_mat)
    lab[okk] <- paste0(sprintf("%+.2f", shift_mat[okk]), "\n", star_mat[okk])
    gm <- pmax(pmin(g_mat, set_g_lim), -set_g_lim)
    rownames(gm) <- rownames(lab) <- clean_lab(rownames(gm))
    pheatmap::pheatmap(gm, display_numbers = lab, fontsize_number = fs_set$number, number_color = "black",
                       cluster_rows = FALSE, cluster_cols = FALSE, gaps_col = length(tissue_order),
                       labels_col = c(tissue_order, pooled_label), color = my_colors,
                       breaks = seq(-set_g_lim, set_g_lim, length.out = 101), na_col = "grey85",
                       border_color = "white", fontsize = fs_set$base, fontsize_row = fs_set$row,
                       fontsize_col = fs_set$col, silent = TRUE)
  }
  
  sc_plot <- set_tbl
  star_sc <- to_mat(mutate(sc_plot, star = paste0(star_fun(q_selfcontained),
                                                  ifelse(!is.na(q_mixed) & q_mixed < fdr_cut, "+", ""))),
                    "Set", "Tissue", "star", sets_ord, col_names, fill = "")
  qm_pooled <- cam_pooled$q_mixed[match(pooled_sc$Set, cam_pooled$Set)]
  star_sc[pooled_sc$Set, "All tissues"] <- paste0(star_fun(pooled_sc$q_pooled_sc),
                                                  ifelse(!is.na(qm_pooled) & qm_pooled < fdr_cut, "+", ""))
  hB <- set_heatmap(star_sc, "All tissues\n(pooled lme)")
  
  star_comp <- to_mat(mutate(set_tbl, star = star_fun(q_camera)), "Set", "Tissue", "star", sets_ord, col_names, fill = "")
  star_comp[cam_pooled$Set, "All tissues"] <- star_fun(cam_pooled$q_camera)
  hC <- set_heatmap(star_comp, "All tissues\n(pooled camera)")
  
  # ---------------- ASSEMBLE ---------------------------------------------------------------------
  pA <- wrap_elements(full = hA$gtable); pB <- wrap_elements(full = hB$gtable); pC <- wrap_elements(full = hC$gtable)
  
  panel_cm <- label_cm + length(sets_ord) * set_row_cm + col_label_cm          # height of B (and of C), in cm
  right_col <- wrap_plots(labelled("B", pB), labelled("C", pC), plot_spacer(), ncol = 1,
                          heights = unit(c(panel_cm, panel_cm, 1), c("cm", "cm", "null")))   # B and C top-aligned; spacer takes any slack
  fig <- wrap_plots(labelled("A", pA), right_col, ncol = 2, widths = panel_widths) +
    plot_annotation(theme = theme(plot.margin = margin(40, 60, 40, 40)))
  
  fig_height <- max(15, 0.24 * nrow(gt) + 4, 2 * panel_cm / 2.54 + 3)
  out_base <- file.path(main_fig_dir, paste0(FILE_PREFIX, "_Composite_Figure"))
  if (save_png)  ggsave(paste0(out_base, ".png"),  fig, width = fig_width, height = fig_height, dpi = 300, bg = "white", limitsize = FALSE)
  if (save_pdf)  ggsave(paste0(out_base, ".pdf"),  fig, width = fig_width, height = fig_height, bg = "white", limitsize = FALSE)
  if (save_tiff) ggsave(paste0(out_base, ".tiff"), fig, width = fig_width, height = fig_height, dpi = 300, bg = "white",
                        compression = "lzw", limitsize = FALSE)
  tibble(pathway = FILE_PREFIX, n_genes = nrow(gt), n_sets = length(sets_ord), width_in = fig_width,
         height_in = round(fig_height, 1), status = "ok")
}

# ---- run all pathways ---------------------------------------------------------------------------------
summary_tbl <- map_dfr(pathways, function(p) {
  message("=== ", p, " ===")
  tryCatch(make_composite(p), error = function(e) {
    message("  FAILED: ", conditionMessage(e))
    tibble(pathway = p, n_genes = NA_integer_, n_sets = NA_integer_, width_in = NA_real_,
           height_in = NA_real_, status = paste("FAILED:", conditionMessage(e)))
  })
})
cat("\n"); print(as.data.frame(summary_tbl), right = FALSE)
cat("\nComposite figures are in:", main_fig_dir, "\n")