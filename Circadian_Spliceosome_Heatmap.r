# ==========================================
# GLOBAL SETUP & DIRECTORIES
# ==========================================
library(KEGGREST)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(tidyverse)
library(pheatmap)
library(here)
library(pathview)

# Define target output directory
output_dir <- "/Users/notarcha/Desktop/PNAS_code/Outputs/Pathways_Heatmaps"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create local genes folder if missing
if (!dir.exists("./genes")) {
  dir.create("./genes")
}

# ==========================================
# GLOBAL HEATMAP STYLING CONFIGURATION
# ==========================================
GLOBAL_FS_ROW  <- 20  # Gene names / Row labels
GLOBAL_FS_COL  <- 20  # Column / Sample labels
GLOBAL_FS_NUM  <- 20  # Significance stars size
GLOBAL_BORDER  <- "white" # Border color

# Color Palette Setup
purple_orange_colors <- colorRampPalette(c("#b35806", "#f1a340", "#fee0b6", "#f7f7f7", "#d8daeb", "#998ec3", "#542788"))(100)

# Custom Heatmap Wrapper Function
plot_custom_heatmap <- function(mat, ...) {
  args <- list(mat = mat, ...)
  
  args$fontsize_row   <- GLOBAL_FS_ROW
  args$fontsize_col   <- GLOBAL_FS_COL
  args$border_color   <- GLOBAL_BORDER
  
  if (is.null(args$fontsize)) args$fontsize <- 12
  
  do.call(pheatmap, args)
}

tissue_order <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Lungs", "Heart", "Liver", "Pectoral")
state_order  <- c("N", "D")

calc_z_score <- function(x) {
  if (sd(x) == 0) return(rep(0, length(x)))
  return((x - mean(x)) / sd(x))
}

# Shared Metadata & Expression Counts (Loaded directly from Results)
RNASeq_metadata  <- read.csv(here("Results", "RNASeq_metadata.csv"), row.names = 1)
Normalizedcounts <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"), row.names = 1)

RNASeq_metadata_clean <- RNASeq_metadata %>% rownames_to_column(var = "SampleID")
RNASeq_metadata_clean$SampleID <- RNASeq_metadata_clean$SampleID %>% str_trim() %>% str_replace_all("-", "_")
colnames(Normalizedcounts) <- colnames(Normalizedcounts) %>% str_trim() %>% str_replace_all("-", "_")

common_samples_global <- intersect(RNASeq_metadata_clean$SampleID, colnames(Normalizedcounts))

################################################################################
## 1. CIRCADIAN / TORPOR ENTRY PATHWAY (Yin & Lazar Model)
################################################################################
cat("\n--- STARTING ANALYSIS: TORPOR ENTRY (Yin & Lazar Model) ---\n")

FILE_PREFIX_CIRC <- "TORPOR_ENTRY"

pillar_info_circ <- data.frame(
  Gene = c(
    "ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "CRY1", "CRY2",
    "BHLHE40", "NR1D1", "RORA", "DBP", "NAMPT", "GSK3B",
    "PDK4", "TXNIP", "FGF21", "PNPLA2", "LIPE", "PFKFB3", "HK2", "PPARGC1A",
    "PCK1", "G6PC", "PPARG", "ELOVL3", "APOC3",
    "PRKAA1", "PRKAA2", "SIRT1", "SIRT3", "HIF1A", "SOD2", "GPX1", "CAT"
  ),
  Function = c(
    rep("Master Timekeepers", 7),
    rep("Circadian-Metabolic Linkers", 6),
    rep("The Fuel Selection Lock", 13),
    rep("Energy & Stress Monitoring", 8)
  )
)

metadata_circ <- RNASeq_metadata_clean %>% 
  filter(SampleID %in% common_samples_global, Metabolic_State != "T") %>%
  mutate(Tissue = recode(Tissue, "Gut1" = "Proximal Gut", "Gut2" = "Medial Gut", "Gut3" = "Distal Gut", "Pect" = "Pectoral"),
         Tissue = factor(Tissue, levels = tissue_order),
         Metabolic_State = factor(Metabolic_State, levels = state_order)) %>%
  arrange(Metabolic_State, Tissue)

counts_circ <- Normalizedcounts[, metadata_circ$SampleID]
ordered_samples_circ <- metadata_circ$SampleID

genes_circ_available <- intersect(rownames(counts_circ), pillar_info_circ$Gene)

scaled_list_circ <- list()
for (tis in levels(metadata_circ$Tissue)) {
  samps <- metadata_circ %>% filter(Tissue == tis) %>% pull(SampleID)
  if (length(samps) == 0) next
  sub_mat <- counts_circ[genes_circ_available, samps]
  scaled_list_circ[[tis]] <- t(apply(sub_mat, 1, calc_z_score))
}
counts_circ_scaled <- do.call(cbind, scaled_list_circ)[, ordered_samples_circ]

pillar_genes_circ <- pillar_info_circ$Gene[pillar_info_circ$Gene %in% rownames(counts_circ_scaled)]

counts_circ_long <- as.data.frame(counts_circ_scaled[pillar_genes_circ, ]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata_circ %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

diff_matrix_circ <- counts_circ_long %>%
  group_by(Gene, Tissue, Metabolic_State) %>%
  summarize(Mean_Z = mean(Z_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Metabolic_State, values_from = Mean_Z) %>%
  mutate(Shift = D - N) %>%
  dplyr::select(Gene, Tissue, Shift) %>%
  pivot_wider(names_from = Tissue, values_from = Shift) %>%
  column_to_rownames("Gene")
diff_matrix_circ <- diff_matrix_circ[pillar_genes_circ, tissue_order]

stats_circ <- counts_circ_long %>%
  group_by(Gene, Tissue) %>%
  summarize(
    p_val = {
      n_N <- sum(Metabolic_State == "N", na.rm = TRUE)
      n_D <- sum(Metabolic_State == "D", na.rm = TRUE)
      if (n_N >= 2 && n_D >= 2) {
        tryCatch(t.test(Z_Score ~ Metabolic_State)$p.value, error = function(e) NA)
      } else { NA }
    },
    .groups = 'drop'
  ) %>%
  mutate(sig_star = ifelse(!is.na(p_val) & p_val < 0.05, "*", ""))

sig_matrix_circ <- stats_circ %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")
sig_matrix_circ <- sig_matrix_circ[rownames(diff_matrix_circ), colnames(diff_matrix_circ)]
sig_matrix_circ[is.na(sig_matrix_circ)] <- ""

annotation_row_circ <- data.frame(
  Function = factor(pillar_info_circ$Function[pillar_info_circ$Gene %in% pillar_genes_circ], 
                    levels = unique(pillar_info_circ$Function)),
  row.names = pillar_genes_circ
)

pillar_colors_circ <- list(
  Function = c(
    "Master Timekeepers"          = "#543005", 
    "Circadian-Metabolic Linkers" = "#bf812d", 
    "The Fuel Selection Lock"     = "#01665e", 
    "Energy & Stress Monitoring"  = "#35978f"
  )
)

bold_rows_circ <- lapply(rownames(diff_matrix_circ), function(x) bquote(bold(.(x))))
bold_cols_circ <- lapply(colnames(diff_matrix_circ), function(x) bquote(bold(.(x))))
gap_locations_circ <- cumsum(table(droplevels(annotation_row_circ$Function))[unique(annotation_row_circ$Function)])

tiff(filename = file.path(output_dir, paste0(FILE_PREFIX_CIRC, "_praw_Heatmap_Mechanistic.tiff")), 
     width = 10, height = 12, units = "in", res = 300, compression = "lzw")

plot_custom_heatmap(
  mat               = as.matrix(diff_matrix_circ),
  display_numbers   = as.matrix(sig_matrix_circ), 
  fontsize_number   = GLOBAL_FS_NUM, 
  labels_row        = as.expression(bold_rows_circ),
  labels_col        = as.expression(bold_cols_circ),
  cluster_rows      = FALSE, 
  cluster_cols      = FALSE,
  annotation_row    = annotation_row_circ,
  annotation_colors = pillar_colors_circ,
  scale             = "none",
  color             = purple_orange_colors,
  breaks            = seq(-2, 2, length.out = 101),
  main              = "Torpor Entry (Yin & Lazar): Torpor Shift (D-N)\n(* = Raw p < 0.05, Purple = Up in Torpor)",
  gaps_row          = gap_locations_circ
)
dev.off()


################################################################################
## 2. SPLICEOSOME PATHWAY (hsa03040)
################################################################################
cat("\n--- STARTING ANALYSIS: SPLICEOSOME (hsa03040) ---\n")

FILE_PREFIX_SPLICE <- "SPLICEOSOME"

pillar_info_splice <- data.frame(
  Gene = c(
    "SNRPB", "SNRPD1", "SNRPD3", "SNRPE", "SNRPF", "SNRPG", "SNRPC", "SNRPB2",
    "LSM3", "LSM4", "LSM5", "LSM6", "LSM7", "LSM8", "SNRNP27", "SNRNP40", "SNU13",
    "SF3A1", "SF3A2", "SF3A3", "SF3B1", "SF3B3", "SF3B4", "SF3B5", "SF3B6", 
    "PHF5A", "U2AF1", "U2SURP", "DDX42", "DDX46", "PUF60",
    "PRPF8", "PRPF18", "PRPF3", "PRPF4", "PRPF6", "PRPF38A", "PRPF38B", "PRPF40A",
    "CDC5L", "PLRG1", "BCAS2", "CTNNBL1", "CWC15", "CWC25", "PPIL1", "SYF2", "CRNKL1",
    "EFTUD2", "SNRNP200", "TXNL4A", "SLU7", "DHX8", "DHX15", "DHX38",
    "AQR", "DDX5", "DDX23", "USP39", "HSPA8", "SUI1", "ZMAT2", "YJU2", "WBP11", "SMNDC1",
    "SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF6", "SRSF7", "SRSF10", "TRA2A", "TRA2B",
    "HNRNPK", "HNRNPU", "HNRNPM", "HNRNPA3", "RBMX", "RBM17", "RBM22", "RBM25", "RBMX2", "CLK1", 
    "CLK2", "CLK3", "CLK4", "RSRP1"
  ),
  Function = c(
    rep("Core snRNPs & Sm Proteins", 17),
    rep("U2 Complex & Branch Point Recognition", 14),
    rep("Catalytic Step & Prp19 Complex", 24),
    rep("RNA Helicases & Assembly Factors", 10),
    rep("SR Proteins & Splicing Regulators", 23)
  )
) %>% distinct(Gene, .keep_all = TRUE)

metadata_splice <- RNASeq_metadata_clean %>% 
  filter(SampleID %in% common_samples_global, Metabolic_State != "T") %>%
  mutate(Tissue = recode(Tissue, "Gut1" = "Proximal Gut", "Gut2" = "Medial Gut", "Gut3" = "Distal Gut", "Pect" = "Pectoral"),
         Tissue = factor(Tissue, levels = tissue_order),
         Metabolic_State = factor(Metabolic_State, levels = state_order)) %>%
  arrange(Metabolic_State, Tissue)

counts_splice <- Normalizedcounts[, metadata_splice$SampleID]
ordered_samples_splice <- metadata_splice$SampleID

genes_splice_available <- intersect(rownames(counts_splice), pillar_info_splice$Gene)

scaled_list_splice <- list()
for (tis in levels(metadata_splice$Tissue)) {
  samps <- metadata_splice %>% filter(Tissue == tis) %>% pull(SampleID)
  if (length(samps) == 0) next
  sub_mat <- counts_splice[genes_splice_available, samps]
  scaled_list_splice[[tis]] <- t(apply(sub_mat, 1, calc_z_score))
}
counts_splice_scaled <- do.call(cbind, scaled_list_splice)[, ordered_samples_splice]

pillar_genes_splice <- pillar_info_splice$Gene[pillar_info_splice$Gene %in% rownames(counts_splice_scaled)]

counts_splice_long <- as.data.frame(counts_splice_scaled[pillar_genes_splice, ]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata_splice %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

diff_matrix_splice <- counts_splice_long %>%
  group_by(Gene, Tissue, Metabolic_State) %>%
  summarize(Mean_Z = mean(Z_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Metabolic_State, values_from = Mean_Z) %>%
  mutate(Shift = D - N) %>%
  dplyr::select(Gene, Tissue, Shift) %>%
  pivot_wider(names_from = Tissue, values_from = Shift) %>%
  column_to_rownames("Gene")
diff_matrix_splice <- diff_matrix_splice[pillar_genes_splice, tissue_order]

stats_splice <- counts_splice_long %>%
  group_by(Gene, Tissue) %>%
  summarize(
    p_val = {
      n_N <- sum(Metabolic_State == "N", na.rm = TRUE)
      n_D <- sum(Metabolic_State == "D", na.rm = TRUE)
      if (n_N >= 2 && n_D >= 2) {
        tryCatch(t.test(Z_Score ~ Metabolic_State)$p.value, error = function(e) NA)
      } else { NA }
    },
    .groups = 'drop'
  ) %>%
  mutate(sig_star = ifelse(!is.na(p_val) & p_val < 0.05, "*", ""))

sig_matrix_splice <- stats_splice %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")
sig_matrix_splice <- sig_matrix_splice[rownames(diff_matrix_splice), colnames(diff_matrix_splice)]
sig_matrix_splice[is.na(sig_matrix_splice)] <- ""

annotation_row_splice <- data.frame(
  Function = factor(pillar_info_splice$Function[pillar_info_splice$Gene %in% pillar_genes_splice], 
                    levels = unique(pillar_info_splice$Function)),
  row.names = pillar_genes_splice
)

pillar_colors_splice <- list(
  Function = c(
    "Core snRNPs & Sm Proteins"              = "#543005", 
    "U2 Complex & Branch Point Recognition" = "#8c510a", 
    "Catalytic Step & Prp19 Complex"        = "#bf812d", 
    "RNA Helicases & Assembly Factors"      = "#01665e",
    "SR Proteins & Splicing Regulators"     = "#35978f"
  )
)

bold_rows_splice <- lapply(rownames(diff_matrix_splice), function(x) bquote(bold(.(x))))
bold_cols_splice <- lapply(colnames(diff_matrix_splice), function(x) bquote(bold(.(x))))
gap_locations_splice <- cumsum(table(droplevels(annotation_row_splice$Function))[unique(annotation_row_splice$Function)])

# --- Output Configuration with Expanded Canvas ---
tiff(filename = file.path(output_dir, paste0(FILE_PREFIX_SPLICE, "_praw_Heatmap_Spliceosome.tiff")), 
     width = 12, height = 24, units = "in", res = 300, compression = "lzw")

plot_custom_heatmap(
  mat               = as.matrix(diff_matrix_splice),
  display_numbers   = as.matrix(sig_matrix_splice), 
  fontsize_number   = 12,                # Adjusted star size to fit tighter rows
  fontsize_row      = 10,                # Reduced row font size to eliminate text overlap
  labels_row        = as.expression(bold_rows_splice), 
  labels_col        = as.expression(bold_cols_splice),
  cluster_rows      = FALSE, 
  cluster_cols      = FALSE,
  annotation_row    = annotation_row_splice,
  annotation_colors = pillar_colors_splice,
  scale             = "none",
  color             = purple_orange_colors,
  breaks            = seq(-2, 2, length.out = 101),
  main              = "Spliceosome: Torpor Shift (D-N)\n(* = Raw p < 0.05, Purple = Up in Torpor)",
  gaps_row          = gap_locations_splice
)
dev.off()

cat("\n--- ALL HEATMAP ANALYSES EXECUTED SUCCESSFULLY ---\n")