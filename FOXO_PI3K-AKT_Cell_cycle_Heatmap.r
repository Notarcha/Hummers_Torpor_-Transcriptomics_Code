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
library(ggpubr)

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

get_symbols_from_ids <- function(id_string) {
  if (is.na(id_string) || id_string == "") return(NA)
  ids <- unlist(strsplit(as.character(id_string), ","))
  syms <- tryCatch({
    mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  }, error = function(e) return(ids))
  syms[is.na(syms)] <- ids[is.na(syms)]
  paste(syms, collapse = ", ")
}

# Load Shared Metadata and Counts from Results folder
RNASeq_metadata  <- read.csv(here("Results", "RNASeq_metadata.csv"), row.names = 1)
Normalizedcounts <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"), row.names = 1)

# Format IDs universally
RNASeq_metadata_clean <- RNASeq_metadata %>% rownames_to_column(var = "SampleID")
RNASeq_metadata_clean$SampleID <- RNASeq_metadata_clean$SampleID %>% str_trim() %>% str_replace_all("-", "_")
colnames(Normalizedcounts) <- colnames(Normalizedcounts) %>% str_trim() %>% str_replace_all("-", "_")

common_samples_global <- intersect(RNASeq_metadata_clean$SampleID, colnames(Normalizedcounts))

my_colors <- colorRampPalette(c("#b35806", "#f1a340", "#fee0b6", "#f7f7f7", "#d8daeb", "#998ec3", "#542788"))(100)

################################################################################
## PATHWAY 1: PI3K-Akt Signaling Pathway (hsa04151)
################################################################################
cat("\n=================== PROCESSING: PI3K-Akt Pathway ===================\n")

PI3K_PATHWAY_ID  <- "path:hsa04151"
FILE_PREFIX_PI3K <- "PI3K_AKT"

# 1. KEGG Gene Extraction
gene_map_pi3k <- keggLink("hsa", PI3K_PATHWAY_ID)
entrez_ids_pi3k <- sub("hsa:", "", unique(gene_map_pi3k))
gene_symbols_pi3k <- mapIds(org.Hs.eg.db, keys = entrez_ids_pi3k, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
kegg_pi3k_genes <- data.frame(symbol = gene_symbols_pi3k, stringsAsFactors = FALSE) %>% na.omit()

# 2. Metadata Alignment & Filtering
metadata_pi3k <- RNASeq_metadata_clean %>% 
  filter(SampleID %in% common_samples_global, Metabolic_State != "T")
counts_pi3k <- Normalizedcounts[, metadata_pi3k$SampleID]

metadata_pi3k <- metadata_pi3k %>%
  mutate(Tissue = recode(Tissue, "Gut1" = "Proximal Gut", "Gut2" = "Medial Gut", "Gut3" = "Distal Gut", "Pect" = "Pectoral"))

# 3. Define PI3K Target Genes & Order
pi3k_pillar_info <- data.frame(
  Gene = c(
    "INSR", "IGF1R", "EGFR", "FGFR1", "MET", "ITGA1", "ITGB1", "PTK2", "IL6R", "JAK1", "STAT3", "TLR2", "TLR4", "CD19", "SYK",
    "PIK3CA", "AKT1", "PDK1", "PTEN", "TSC2",
    "MTOR", "RHEB", "RPS6KB1", "EIF4EBP1", "EIF4E", "EIF4B",
    "BCL2", "MCL1", "BAD", "BAX", "MDM2", "CASP9",
    "GSK3B", "FOXO1", "FOXO3", "PCK1", "G6PC", "SREBF1",
    "MYC", "CCND1", "CDKN1A", "CDKN1B", "BRCA1"
  ),
  Function = c(
    rep("Sensors & Receptors", 15),
    rep("Core Hub & Brakes", 5),
    rep("Protein Synthesis (mTOR)", 6),
    rep("Survival & Apoptosis", 6),
    rep("Metabolic Control (FOXO/GSK)", 6),
    rep("Cell Cycle & Repair", 5)
  ),
  stringsAsFactors = FALSE
)

combined_pi3k_genes <- unique(c(kegg_pi3k_genes$symbol, pi3k_pillar_info$Gene))
genes_of_interest_pi3k <- intersect(rownames(counts_pi3k), combined_pi3k_genes)

metadata_pi3k <- metadata_pi3k %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_order),
    Metabolic_State = factor(Metabolic_State, levels = state_order)
  ) %>% arrange(Metabolic_State, Tissue)

ordered_samples_pi3k <- metadata_pi3k$SampleID

# 4. Per-Tissue Scaling
scaled_list_pi3k <- list()
for (tis in levels(metadata_pi3k$Tissue)) {
  samps <- metadata_pi3k %>% filter(Tissue == tis) %>% pull(SampleID)
  if (length(samps) == 0) next
  sub_mat <- counts_pi3k[genes_of_interest_pi3k, samps]
  scaled_list_pi3k[[tis]] <- t(apply(sub_mat, 1, calc_z_score))
}
counts_pi3k_scaled <- do.call(cbind, scaled_list_pi3k)[, ordered_samples_pi3k]

# 5. Long-Format Conversion & Statistics
pi3k_available <- pi3k_pillar_info$Gene[pi3k_pillar_info$Gene %in% rownames(counts_pi3k_scaled)]

pi3k_long <- as.data.frame(counts_pi3k_scaled[pi3k_available, ]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata_pi3k %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

pi3k_diff_matrix <- pi3k_long %>%
  group_by(Gene, Tissue, Metabolic_State) %>%
  summarize(Mean_Z = mean(Z_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Metabolic_State, values_from = Mean_Z) %>%
  mutate(Shift = D - N) %>%
  dplyr::select(Gene, Tissue, Shift) %>%
  pivot_wider(names_from = Tissue, values_from = Shift) %>%
  column_to_rownames("Gene")
pi3k_diff_matrix <- pi3k_diff_matrix[pi3k_available, tissue_order]

pi3k_stats <- pi3k_long %>%
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

pi3k_sig_matrix <- pi3k_stats %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")
pi3k_sig_matrix <- pi3k_sig_matrix[rownames(pi3k_diff_matrix), colnames(pi3k_diff_matrix)]
pi3k_sig_matrix[is.na(pi3k_sig_matrix)] <- ""

# 6. Heatmap Output
annotation_row_pi3k <- data.frame(
  Function = factor(pi3k_pillar_info$Function[pi3k_pillar_info$Gene %in% pi3k_available],
                    levels = unique(pi3k_pillar_info$Function)),
  row.names = pi3k_available
)

cat_colors_pi3k <- list(
  Function = c(
    "Sensors & Receptors"          = "#43a2ca", 
    "Core Hub & Brakes"            = "#636363", 
    "Protein Synthesis (mTOR)"     = "#e6550d", 
    "Survival & Apoptosis"         = "#de2d26", 
    "Metabolic Control (FOXO/GSK)" = "#756bb1", 
    "Cell Cycle & Repair"          = "#31a354"  
  )
)

bold_rows_pi3k <- lapply(rownames(pi3k_diff_matrix), function(x) bquote(bold(.(x))))
bold_cols_pi3k <- lapply(colnames(pi3k_diff_matrix), function(x) bquote(bold(.(x))))
gap_indices_pi3k <- cumsum(table(droplevels(annotation_row_pi3k$Function))[unique(annotation_row_pi3k$Function)])

out_file_pi3k <- file.path(output_dir, paste0(FILE_PREFIX_PI3K, "_Waterfall_praw_Stars.tiff"))
tiff(filename = out_file_pi3k, width = 11, height = 13, units = "in", res = 300, compression = "lzw")

plot_custom_heatmap(
  mat               = as.matrix(pi3k_diff_matrix),
  cluster_rows      = FALSE, 
  cluster_cols      = FALSE,
  display_numbers   = as.matrix(pi3k_sig_matrix), 
  fontsize_number   = GLOBAL_FS_NUM, 
  number_color      = "black",
  labels_row        = as.expression(bold_rows_pi3k),
  labels_col        = as.expression(bold_cols_pi3k),
  annotation_row    = annotation_row_pi3k, 
  annotation_colors = cat_colors_pi3k,
  gaps_row          = gap_indices_pi3k,
  scale             = "none", 
  color             = my_colors, 
  breaks            = seq(-1.5, 1.5, length.out = 101),
  main              = "PI3K-AKT Pathway: Torpor Shift (D-N)\n(* = Raw p < 0.05)"
)

dev.off()

################################################################################
## PATHWAY 2: Cell Cycle Pathway (hsa04110)
################################################################################
cat("\n=================== PROCESSING: Cell Cycle Pathway ===================\n")

CELL_CYCLE_PATHWAY_ID <- "path:hsa04110"
FILE_PREFIX_CC        <- "CELL_CYCLE"

# 1. KEGG Gene Extraction
gene_map_cc <- keggLink("hsa", CELL_CYCLE_PATHWAY_ID)
entrez_ids_cc <- sub("hsa:", "", unique(gene_map_cc))
gene_symbols_cc <- mapIds(org.Hs.eg.db, keys = entrez_ids_cc, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
kegg_cc_genes <- data.frame(symbol = gene_symbols_cc, stringsAsFactors = FALSE) %>% na.omit()

# 2. Metadata Alignment & Filtering
metadata_cc <- RNASeq_metadata_clean %>% 
  filter(SampleID %in% common_samples_global, Metabolic_State != "T")
counts_cc   <- Normalizedcounts[, metadata_cc$SampleID]

metadata_cc <- metadata_cc %>%
  mutate(Tissue = recode(Tissue, "Gut1" = "Proximal Gut", "Gut2" = "Medial Gut", "Gut3" = "Distal Gut", "Pect" = "Pectoral"))

# 3. Define Cell Cycle Pillar Information
cell_cycle_pillar_info <- data.frame(
  Gene = c(
    "RB1", "E2F1", "TFDP1", "CCND1", "CCNE1", "CDK4", "CDK6", "CDK2",
    "E2F4", "E2F5", "RBL1", "HDAC1",
    "CDKN1A", "CDKN1B", "GSK3B", "TGFB2", "SMAD4",
    "PCNA", "CDC6", "CCNA2", "CCNB1", "CDK1",
    "TP53", "ATM", "CHEK1", "GADD45A", "SFN",
    "CDC20", "FZR1", "MAD2L1", "BUB1"
  ),
  Function = c(
    rep("G1/S Phase Entry", 8),
    rep("DREAM Complex (G0)", 4),
    rep("CKIs & TGF-b Signaling", 5),
    rep("S-Phase & G2/M Transition", 5),
    rep("DDR & DNA Surveillance", 5),
    rep("APC/C & Mitotic Exit", 4)
  ),
  stringsAsFactors = FALSE
)

combined_cc_genes <- unique(c(kegg_cc_genes$symbol, cell_cycle_pillar_info$Gene))
genes_of_interest_cc <- intersect(rownames(counts_cc), combined_cc_genes)

metadata_cc <- metadata_cc %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_order),
    Metabolic_State = factor(Metabolic_State, levels = state_order)
  ) %>% arrange(Metabolic_State, Tissue)

ordered_samples_cc <- metadata_cc$SampleID

# 4. Per-Tissue Scaling
scaled_list_cc <- list()
for (tis in levels(metadata_cc$Tissue)) {
  samps <- metadata_cc %>% filter(Tissue == tis) %>% pull(SampleID)
  if (length(samps) == 0) next
  sub_mat <- counts_cc[genes_of_interest_cc, samps]
  scaled_list_cc[[tis]] <- t(apply(sub_mat, 1, calc_z_score))
}
counts_cc_scaled <- do.call(cbind, scaled_list_cc)[, ordered_samples_cc]

# 5. Long-Format Conversion & Statistics
cc_available <- cell_cycle_pillar_info$Gene[cell_cycle_pillar_info$Gene %in% rownames(counts_cc_scaled)]

cc_long <- as.data.frame(counts_cc_scaled[cc_available, ]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata_cc %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

cc_diff_matrix <- cc_long %>%
  group_by(Gene, Tissue, Metabolic_State) %>%
  summarize(Mean_Z = mean(Z_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Metabolic_State, values_from = Mean_Z) %>%
  mutate(Shift = D - N) %>%
  dplyr::select(Gene, Tissue, Shift) %>%
  pivot_wider(names_from = Tissue, values_from = Shift) %>%
  column_to_rownames("Gene")
cc_diff_matrix <- cc_diff_matrix[cc_available, tissue_order]

cc_stats <- cc_long %>%
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

cc_sig_matrix <- cc_stats %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")
cc_sig_matrix <- cc_sig_matrix[rownames(cc_diff_matrix), colnames(cc_diff_matrix)]
cc_sig_matrix[is.na(cc_sig_matrix)] <- ""

# 6. Heatmap Output
annotation_row_cc <- data.frame(
  Function = factor(cell_cycle_pillar_info$Function[cell_cycle_pillar_info$Gene %in% cc_available],
                    levels = unique(cell_cycle_pillar_info$Function)),
  row.names = cc_available
)

cat_colors_cc <- list(
  Function = c(
    "G1/S Phase Entry"          = "#e41a1c",
    "DREAM Complex (G0)"        = "#377eb8",
    "CKIs & TGF-b Signaling"   = "#4daf4a",
    "S-Phase & G2/M Transition" = "#984ea3",
    "DDR & DNA Surveillance"    = "#ff7f00",
    "APC/C & Mitotic Exit"      = "#ffff33"
  )
)

bold_rows_cc <- lapply(rownames(cc_diff_matrix), function(x) bquote(bold(.(x))))
bold_cols_cc <- lapply(colnames(cc_diff_matrix), function(x) bquote(bold(.(x))))
gap_indices_cc <- cumsum(table(droplevels(annotation_row_cc$Function))[unique(annotation_row_cc$Function)])

out_file_cc <- file.path(output_dir, paste0(FILE_PREFIX_CC, "_Waterfall_praw_Stars.tiff"))
tiff(filename = out_file_cc, width = 11, height = 11, units = "in", res = 300, compression = "lzw")

plot_custom_heatmap(
  mat               = as.matrix(cc_diff_matrix),
  cluster_rows      = FALSE, 
  cluster_cols      = FALSE,
  display_numbers   = as.matrix(cc_sig_matrix), 
  fontsize_number   = GLOBAL_FS_NUM, 
  number_color      = "black",
  labels_row        = as.expression(bold_rows_cc),
  labels_col        = as.expression(bold_cols_cc),
  annotation_row    = annotation_row_cc, 
  annotation_colors = cat_colors_cc,
  gaps_row          = gap_indices_cc,
  scale             = "none", 
  color             = my_colors, 
  breaks            = seq(-1.5, 1.5, length.out = 101),
  main              = "Cell Cycle Pathway: Torpor Shift (D-N)\n(* = Raw p < 0.05)"
)

dev.off()

################################################################################
## PATHWAY 3: FOXO SIGNALING PATHWAY
################################################################################
cat("\n=================== PROCESSING: FOXO Signaling Pathway ===================\n")

FILE_PREFIX_FOXO <- "FOXO_SIGNALING"

# 1. Functional Pillars
pillar_info_foxo <- data.frame(
  Gene = c(
    "IGF1R", "INSR", "PIK3CA", "PIK3CB", "AKT1", "AKT2", "AKT3", "SGK1",
    "PRKAA1", "PRKAA2", "STK4", "MAPK8", "MAPK9",
    "FOXO1", "FOXO3", "FOXO4",
    "CDKN1B", "CDKN1A", "CCNG2",
    "PDK4", "PCK1", "G6PC", "PPARG", "PPARGC1A",
    "BNIP3", "GABARAPL1", "ATG12",
    "SOD2", "CAT", "GADD45A"
  ),
  Function = c(
    rep("Upstream Inhibitors (Growth)", 8), 
    rep("Upstream Activators (Stress)", 5),
    rep("The FoxO Transcription Factors", 3),
    rep("Cell Cycle Arrest", 3),
    rep("Metabolic Switching", 5),
    rep("Autophagy & Cleanup", 3),
    rep("Antioxidant Defense", 3)
  )
)

# 2. Metadata & Counts
metadata_foxo <- RNASeq_metadata_clean %>% 
  filter(SampleID %in% common_samples_global, Metabolic_State != "T") %>%
  mutate(Tissue = recode(Tissue, "Gut1" = "Proximal Gut", "Gut2" = "Medial Gut", "Gut3" = "Distal Gut", "Pect" = "Pectoral"),
         Tissue = factor(Tissue, levels = tissue_order),
         Metabolic_State = factor(Metabolic_State, levels = state_order)) %>%
  arrange(Metabolic_State, Tissue)

counts_foxo <- Normalizedcounts[, metadata_foxo$SampleID]
ordered_samples_foxo <- metadata_foxo$SampleID

# 3. Z-Score Scaling
genes_found_foxo <- intersect(pillar_info_foxo$Gene, rownames(counts_foxo))
pillar_info_foxo_filtered <- pillar_info_foxo %>% filter(Gene %in% genes_found_foxo)

scaled_list_foxo <- list()
for (tis in levels(metadata_foxo$Tissue)) {
  samps <- metadata_foxo %>% filter(Tissue == tis) %>% pull(SampleID)
  if (length(samps) == 0) next
  sub_mat <- counts_foxo[genes_found_foxo, samps]
  scaled_list_foxo[[tis]] <- t(apply(sub_mat, 1, calc_z_score))
}
counts_foxo_scaled <- do.call(cbind, scaled_list_foxo)[, ordered_samples_foxo]

# 4. Long Format & Metrics
counts_foxo_long <- as.data.frame(counts_foxo_scaled) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata_foxo %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

foxo_summary <- counts_foxo_long %>%
  group_by(Gene, Tissue) %>%
  summarize(
    Shift = mean(Z_Score[Metabolic_State == "D"], na.rm = TRUE) - 
      mean(Z_Score[Metabolic_State == "N"], na.rm = TRUE),
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

diff_matrix_foxo <- foxo_summary %>%
  dplyr::select(Gene, Tissue, Shift) %>%
  pivot_wider(names_from = Tissue, values_from = Shift) %>%
  column_to_rownames("Gene")

sig_matrix_foxo <- foxo_summary %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")

diff_matrix_foxo <- diff_matrix_foxo[pillar_info_foxo_filtered$Gene, tissue_order]
sig_matrix_foxo  <- sig_matrix_foxo[rownames(diff_matrix_foxo), colnames(diff_matrix_foxo)]
sig_matrix_foxo[is.na(sig_matrix_foxo)] <- ""

# 5. Annotations & Plotting
ann_row_foxo <- data.frame(
  Function = factor(pillar_info_foxo_filtered$Function, levels = unique(pillar_info_foxo_filtered$Function)),
  row.names = pillar_info_foxo_filtered$Gene
)

pillar_colors_foxo <- list(
  Function = c(
    "Upstream Inhibitors (Growth)"      = "#543005", "Upstream Activators (Stress)"     = "#8c510a",
    "The FoxO Transcription Factors"    = "#01665e", "Cell Cycle Arrest"                = "#35978f",
    "Metabolic Switching"               = "#80cdc1", "Autophagy & Cleanup"              = "#4d4d4d",
    "Antioxidant Defense"               = "#1a1a1a"
  )
)

gap_locations_foxo <- cumsum(table(droplevels(ann_row_foxo$Function))[unique(ann_row_foxo$Function)])
bold_rows_foxo <- lapply(rownames(diff_matrix_foxo), function(x) bquote(bold(.(x))))
bold_cols_foxo <- lapply(colnames(diff_matrix_foxo), function(x) bquote(bold(.(x))))

# --- FOXO Heatmap 1: Red/Blue ---
red_blue_colors <- colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7", "#ef8a62", "#b2182b"))(100)

tiff(filename = file.path(output_dir, paste0(FILE_PREFIX_FOXO, "_praw_RedBlue.tiff")), 
     width = 11, height = 12, units = "in", res = 300, compression = "lzw")

plot_custom_heatmap(
  mat               = as.matrix(diff_matrix_foxo),
  display_numbers   = as.matrix(sig_matrix_foxo),
  fontsize_number   = GLOBAL_FS_NUM,
  labels_row        = as.expression(bold_rows_foxo), 
  labels_col        = as.expression(bold_cols_foxo),
  cluster_rows      = FALSE, 
  cluster_cols      = FALSE,
  annotation_row    = ann_row_foxo, 
  annotation_colors = pillar_colors_foxo,
  scale             = "none", 
  color             = red_blue_colors,
  breaks            = seq(-2, 2, length.out = 101),
  gaps_row          = gap_locations_foxo,
  main              = "FOXO Signaling: Torpor Shift (D-N)\n(* = Raw p < 0.05)",
  angle_col         = "45"
)
dev.off()

# --- FOXO Heatmap 2: Purple/Orange ---
purple_orange_colors <- colorRampPalette(c("#b35806", "#f1a340","#fee0b6","#f7f7f7","#d8daeb","#998ec3","#542788"))(100)

tiff(filename = file.path(output_dir, paste0(FILE_PREFIX_FOXO, "_praw_PurpleOrange.tiff")), 
     width = 11, height = 12, units = "in", res = 300, compression = "lzw")

plot_custom_heatmap(
  mat               = as.matrix(diff_matrix_foxo),
  display_numbers   = as.matrix(sig_matrix_foxo),
  fontsize_number   = GLOBAL_FS_NUM,
  labels_row        = as.expression(bold_rows_foxo), 
  labels_col        = as.expression(bold_cols_foxo),
  cluster_rows      = FALSE, 
  cluster_cols      = FALSE,
  annotation_row    = ann_row_foxo, 
  annotation_colors = pillar_colors_foxo,
  scale             = "none", 
  color             = purple_orange_colors,
  breaks            = seq(-2, 2, length.out = 101),
  gaps_row          = gap_locations_foxo,
  main              = "FOXO Signaling: Torpor Shift (D-N)\n(* = Raw p < 0.05, Purple = Up in Torpor)",
  angle_col         = "45"
)
dev.off()

################################################################################
## LINEAR REGRESSION ANALYSIS: AKT/FOXO/SGK PARALOG PAIRS
################################################################################
cat("\n=================== PROCESSING: Linear Regressions ===================\n")

target_genes <- c("AKT1", "AKT3", "FOXO1", "FOXO3", "SGK1")

foxo_sgk_sub <- Normalizedcounts[intersect(target_genes, rownames(Normalizedcounts)), common_samples_global]

metadata_reg <- RNASeq_metadata_clean %>% 
  filter(SampleID %in% common_samples_global, Metabolic_State != "T") %>%
  mutate(
    Tissue = recode(Tissue, "Gut1" = "Proximal Gut", "Gut2" = "Medial Gut", "Gut3" = "Distal Gut", "Pect" = "Pectoral"),
    Tissue = factor(Tissue, levels = tissue_order),
    Metabolic_State = factor(Metabolic_State, levels = state_order)
  )

# Per-tissue Z-score calculation
scaled_reg_list <- list()
for (tis in levels(metadata_reg$Tissue)) {
  samps <- metadata_reg %>% filter(Tissue == tis) %>% pull(SampleID)
  if (length(samps) == 0) next
  sub_mat <- foxo_sgk_sub[, samps]
  scaled_reg_list[[tis]] <- t(apply(sub_mat, 1, calc_z_score))
}

reg_scaled_mat <- do.call(cbind, scaled_reg_list)[, metadata_reg$SampleID]

# Reshape to wide format with one column per gene
reg_data_wide <- as.data.frame(t(reg_scaled_mat)) %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata_reg %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

gene_pairs <- list(
  c(x = "AKT1", y = "FOXO1"),
  c(x = "AKT1", y = "FOXO3"),
  c(x = "AKT3", y = "SGK1")
)

reg_plots <- list()

for (pair in gene_pairs) {
  gene_x <- pair["x"]
  gene_y <- pair["y"]
  
  plot_title <- paste(gene_y, "vs.", gene_x)
  
  p <- ggplot(reg_data_wide, aes_string(x = gene_x, y = gene_y)) +
    geom_smooth(method = "lm", color = "black", fill = "grey80", se = TRUE, linewidth = 0.8) +
    geom_point(aes(color = Tissue, shape = Metabolic_State), size = 3, alpha = 0.85) +
    stat_cor(
      aes(label = paste(..rr.label.., ..p.label.., sep = "~`,\n`~")),
      label.x.npc = "left", 
      label.y.npc = "top",
      size = 4
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size = 12) +
    labs(
      title = plot_title,
      x = paste(gene_x, "(Z-score)"),
      y = paste(gene_y, "(Z-score)"),
      color = "Tissue",
      shape = "Metabolic State"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  reg_plots[[paste(gene_x, gene_y, sep = "_")]] <- p
}

# Save Individual Plots and Combined Panel TIFF
for (pair_name in names(reg_plots)) {
  out_path <- file.path(output_dir, paste0("Regression_", pair_name, ".tiff"))
  ggsave(out_path, plot = reg_plots[[pair_name]], width = 6, height = 5, dpi = 300, compression = "lzw")
}

combined_reg_plot <- ggarrange(
  plotlist = reg_plots, 
  ncol = 3, 
  nrow = 1, 
  common.legend = TRUE, 
  legend = "right"
)

combined_out_path <- file.path(output_dir, "Combined_AKT_FOXO_SGK_Regressions.tiff")
ggsave(combined_out_path, plot = combined_reg_plot, width = 15, height = 5, dpi = 300, compression = "lzw")

cat("\n--- SUCCESS: ALL HEATMAPS AND REGRESSION PLOTS SAVED TO OUTPUTS DIRECTORY ---\n")