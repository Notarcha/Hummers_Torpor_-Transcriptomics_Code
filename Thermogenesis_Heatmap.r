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
  
  # Forces legends, annotations, etc., to default base size
  if (is.null(args$fontsize)) args$fontsize <- 12
  
  do.call(pheatmap, args)
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
  
  # Forces legends, annotations, etc., to default base size
  if (is.null(args$fontsize)) args$fontsize <- 12
  
  do.call(pheatmap, args)
}

######## THERMOGENESIS ###########

#########################################
## 1. LOAD LIBRARIES & INITIALIZE
#########################################
library(KEGGREST)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(tidyverse)
library(pheatmap)
library(here)

# Define main figures relative output directory
main_fig_dir <- here("Outputs", "Pathways_Heatmaps")
if (!dir.exists(main_fig_dir)) {
  dir.create(main_fig_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created directory:", main_fig_dir, "\n")
}

cat("--- STARTING INTUITION-DRIVEN THERMOGENESIS ANALYSIS (N vs D ONLY) ---\n")

genes_dir <- here("genes")
if (!dir.exists(genes_dir)) {
  dir.create(genes_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created directory:", genes_dir, "\n")
}

FILE_PREFIX <- "THERMO_INTUITION"

#########################################
## 2. LOAD & CLEAN DATA (STRICT LOGIC)
#########################################
RNASeq_metadata <- read.csv(here("Results", "RNASeq_metadata.csv"), row.names = 1)
Normalizedcounts <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"), row.names = 1)

metadata <- RNASeq_metadata %>% rownames_to_column(var = "SampleID")
metadata$SampleID <- metadata$SampleID %>% str_trim() %>% str_replace_all("-", "_")
colnames(Normalizedcounts) <- colnames(Normalizedcounts) %>% str_trim() %>% str_replace_all("-", "_")

common_samples <- intersect(metadata$SampleID, colnames(Normalizedcounts))
metadata <- metadata %>% filter(SampleID %in% common_samples)
Normalizedcounts <- Normalizedcounts[, common_samples]

# Filter Step: Remove "T"
metadata <- metadata %>% filter(Metabolic_State != "T")
Normalizedcounts <- Normalizedcounts[, metadata$SampleID]

metadata <- metadata %>%
  mutate(Tissue = recode(Tissue, 
                         "Gut1" = "Proximal Gut", 
                         "Gut2" = "Medial Gut", 
                         "Gut3" = "Distal Gut", 
                         "Pect" = "Pectoral"))

#########################################
## 3. DEFINE GENES & PILLARS (INTUITION-BASED)
#########################################
pillar_info <- data.frame(
  Gene = c(
    # 1. The Hidden Furnace (Calcium Shunts & Uncoupling)
    "SLN", "PLN", "CASQ1", "CASQ2", "ATP2A1", "ATP2A2", "RYR1", "RYR2", "UCP2", "UCP3",
    
    # 2. Redox Shield (Protecting the High-Flux Liver)
    "SOD1", "SOD2", "GPX1", "GPX4", "CAT", "PRDX1", "TXN", "TXNRD1",
    
    # 3. Mitochondrial Grooming & Quality Control
    "PINK1", "PRKN", "DNM1L", "MFN1", "MFN2", "TFAM", "PPARGC1A",
    
    # 4. Adrenergic Signaling (The Master Switch)
    "ADRB3", "GNAS", "ADCY1", "ADCY3", "ADCY6", "PRKACB", "CREB1", "ATF2",
    
    # 5. Mitochondrial ETC Hub (Expanded for Liver Patterns)
    "NDUFA1", "NDUFA9", "NDUFB8", "NDUFS1", "SDHA", "SDHB", "SDHC", "SDHD", 
    "COX4I1", "ATP5F1A", "CYC1", "UQCRC1",
    
    # 6. mTOR & Growth Suppression (Energy Saving)
    "MTOR", "RPTOR", "TSC1", "TSC2", "RPS6", "EIF4EBP1", "MAPK14"
  ),
  Function = c(
    rep("The Hidden Furnace (Ca2+/Uncoupling)", 10),
    rep("Redox Shield (Antioxidants)", 8),
    rep("Mito-Grooming & Quality Control", 7),
    rep("Adrenergic Signaling", 8),
    rep("Mitochondrial ETC Hub", 12),
    rep("mTOR & Growth Suppression", 7)
  )
)

pillar_info <- pillar_info %>% distinct(Gene, .keep_all = TRUE)
genes_of_interest <- intersect(rownames(Normalizedcounts), pillar_info$Gene)

tissue_order <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Lungs", "Heart", "Liver", "Pectoral")
state_order  <- c("N", "D")

metadata <- metadata %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order),
         Metabolic_State = factor(Metabolic_State, levels = state_order))

#########################################
## 4 & 5. SORT & Z-SCORE
#########################################
metadata <- metadata %>% arrange(Metabolic_State, Tissue)
ordered_samples_state_first <- metadata$SampleID

calc_z_score <- function(x) {
  if(sd(x) == 0) return(rep(0, length(x)))
  return((x - mean(x)) / sd(x))
}

tissues <- levels(metadata$Tissue)
scaled_list <- list()

for (tis in tissues) {
  samps_in_tissue <- metadata %>% filter(Tissue == tis) %>% pull(SampleID)
  if(length(samps_in_tissue) == 0) next
  sub_mat <- Normalizedcounts[genes_of_interest, samps_in_tissue]
  scaled_sub <- t(apply(sub_mat, 1, calc_z_score))
  scaled_list[[tis]] <- scaled_sub
}

counts_per_tissue_scaled <- do.call(cbind, scaled_list)
counts_final_sorted <- counts_per_tissue_scaled[, ordered_samples_state_first]

#########################################
## 6. CALCULATE SHIFT (D - N) & RAW P-VALUES
#########################################
pillar_genes_available <- pillar_info$Gene[pillar_info$Gene %in% rownames(counts_final_sorted)]

# 1. Prepare Long Format Data
counts_long <- as.data.frame(counts_final_sorted[pillar_genes_available, ]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

# 2. Calculate Mean Shifts (The Heatmap Colors)
diff_matrix <- counts_long %>%
  group_by(Gene, Tissue, Metabolic_State) %>%
  summarize(Mean_Z = mean(Z_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Metabolic_State, values_from = Mean_Z) %>%
  mutate(Shift = D - N) %>%
  dplyr::select(Gene, Tissue, Shift) %>%
  pivot_wider(names_from = Tissue, values_from = Shift) %>%
  column_to_rownames("Gene")

diff_matrix <- diff_matrix[pillar_genes_available, tissue_order]

# 3. Calculate Stats and apply RAW p-value analysis mapping
stats_summary <- counts_long %>%
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

# 4. Generate Significance Matrix
sig_matrix <- stats_summary %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")

# Re-align sig_matrix with diff_matrix
sig_matrix <- sig_matrix[rownames(diff_matrix), colnames(diff_matrix)]
sig_matrix[is.na(sig_matrix)] <- ""

#########################################
## 7. HEATMAP PLOTTING
#########################################
annotation_row <- data.frame(
  Function = factor(pillar_info$Function[pillar_info$Gene %in% pillar_genes_available], 
                    levels = unique(pillar_info$Function))
)
rownames(annotation_row) <- pillar_genes_available

pillar_colors <- list(
  Function = c(
    "The Hidden Furnace (Ca2+/Uncoupling)" = "#543005", 
    "Redox Shield (Antioxidants)"          = "#8c510a", 
    "Mito-Grooming & Quality Control"      = "#bf812d", 
    "Adrenergic Signaling"                 = "#01665e", 
    "Mitochondrial ETC Hub"                = "#35978f", 
    "mTOR & Growth Suppression"            = "#1a1a1a"  
  )
)

my_colors <- colorRampPalette(c("#b35806", "#f1a340","#fee0b6","#f7f7f7","#d8daeb","#998ec3","#542788"))(100)

# Save to relative figures folder using high-resolution tiff with LZW compression
tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_praw_Heatmap.tiff")), 
     width = 11, height = 18, units = "in", res = 300, compression = "lzw")

bold_rows <- lapply(rownames(diff_matrix), function(x) bquote(bold(.(x))))
bold_cols <- lapply(colnames(diff_matrix), function(x) bquote(bold(.(x))))

plot_custom_heatmap(diff_matrix,
                    display_numbers = sig_matrix, 
                    fontsize_number = GLOBAL_FS_NUM, 
                    labels_row = as.expression(bold_rows), 
                    labels_col = as.expression(bold_cols),
                    cluster_rows = FALSE, cluster_cols = FALSE,
                    annotation_row = annotation_row,
                    annotation_colors = pillar_colors,
                    scale = "none",
                    color = my_colors,
                    breaks = seq(-2, 2, length.out = 101),
                    main = "Thermogenesis: Torpor Shift (D-N)\n(* = Raw p < 0.05, Purple = Up in Torpor)",
                    gaps_row = cumsum(table(annotation_row$Function)[unique(annotation_row$Function)]))

dev.off()

cat(paste0("Success! Thermogenesis raw p-value heatmap saved to: ", main_fig_dir, "\n"))