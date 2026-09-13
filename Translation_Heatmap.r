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


###### Translational machinery ##### 
#########################################
## LOAD LIBRARIES
#########################################
library(KEGGREST)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(tidyverse)
library(pheatmap)
library(here)
library(pathview)

# Define relative output directory (creates full path if missing)
main_fig_dir <- here("Outputs", "Pathways_Heatmaps")
if (!dir.exists(main_fig_dir)) {
  dir.create(main_fig_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created directory:", main_fig_dir, "\n")
}

cat("--- STARTING ANALYSIS (N vs D ONLY) FOR RIBOSOME/TRANSLATION PATHWAY (hsa03010) ---\n")

# Ensure the relative './genes' output directory exists
genes_dir <- here("genes")
if (!dir.exists(genes_dir)) {
  dir.create(genes_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created directory:", genes_dir, "\n")
}

# Define Translation (Ribosome) Pathway ID
TRANSLATION_PATHWAY_ID <- "path:hsa03010"
TRANSLATION_HSA_ID <- "hsa03010"
FILE_PREFIX <- "Translation"

#########################################
## 1. PREPARE KEGG GENE LIST (TRANSLATION)
#########################################
pathway_id <- TRANSLATION_PATHWAY_ID 
gene_map <- keggLink("hsa", pathway_id)
genes_in_pathway <- unique(gene_map)

# Map to Symbols
entrez_ids <- sub("hsa:", "", genes_in_pathway)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

kegg_pathway_genes <- data.frame(
  kegg_id = genes_in_pathway,
  entrez_id = entrez_ids,
  symbol = gene_symbols,
  stringsAsFactors = FALSE
)
kegg_pathway_genes <- na.omit(kegg_pathway_genes)

cat("Step 1: Mapped", nrow(kegg_pathway_genes), "Ribosome pathway genes.\n")

#########################################
## 2. LOAD & CLEAN DATA
#########################################
RNASeq_metadata <- read.csv(here("Results", "RNASeq_metadata.csv"), row.names = 1)
Normalizedcounts <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"), row.names = 1)

metadata <- RNASeq_metadata %>% rownames_to_column(var = "SampleID")
metadata$SampleID <- metadata$SampleID %>% str_trim() %>% str_replace_all("-", "_")
colnames(Normalizedcounts) <- colnames(Normalizedcounts) %>% str_trim() %>% str_replace_all("-", "_")

common_samples <- intersect(metadata$SampleID, colnames(Normalizedcounts))
metadata <- metadata %>% filter(SampleID %in% common_samples)
Normalizedcounts <- Normalizedcounts[, common_samples]

# Filter "T" state
metadata <- metadata %>% filter(Metabolic_State != "T")
Normalizedcounts <- Normalizedcounts[, metadata$SampleID]

# Recode Tissue Names
metadata <- metadata %>%
  mutate(Tissue = recode(Tissue, 
                         "Gut1" = "Proximal Gut", 
                         "Gut2" = "Medial Gut", 
                         "Gut3" = "Distal Gut", 
                         "Pect" = "Pectoral"))

#########################################
## 3. DEFINE TRANSLATION PILLARS (NEW)
#########################################

# Defining logical groupings for the Ribosomal machinery
rpl_genes <- grep("^RPL", rownames(Normalizedcounts), value = TRUE) # Large Subunit
rps_genes <- grep("^RPS", rownames(Normalizedcounts), value = TRUE) # Small Subunit
mrpl_genes <- grep("^MRPL", rownames(Normalizedcounts), value = TRUE) # Mito Large
mrps_genes <- grep("^MRPS", rownames(Normalizedcounts), value = TRUE) # Mito Small
eif_genes <- grep("^EIF", rownames(Normalizedcounts), value = TRUE)  # Initiation Factors

all_custom_pillars <- c(rpl_genes, rps_genes, mrpl_genes, mrps_genes, eif_genes)

combined_target_genes <- unique(c(kegg_pathway_genes$symbol, all_custom_pillars))
genes_of_interest <- intersect(rownames(Normalizedcounts), combined_target_genes)

tissue_order <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Lungs", "Heart", "Liver", "Pectoral")
state_order  <- c("N", "D")

metadata <- metadata %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_order),
    Metabolic_State = factor(Metabolic_State, levels = state_order)
  )

#########################################
## 4. SORT METADATA
#########################################
metadata <- metadata %>% arrange(Metabolic_State, Tissue)
ordered_samples_state_first <- metadata$SampleID

#########################################
## 5. PER-TISSUE Z-SCORE SCALING
#########################################
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
## 7. PREPARE HEATMAP ANNOTATION
#########################################
annotation_col <- metadata %>%
  dplyr::select(SampleID, Metabolic_State, Tissue) %>%
  column_to_rownames("SampleID")

ann_colors <- list(
  Metabolic_State = c(N = "gray", D = "gray30"),
  Tissue = c(
    "Proximal Gut" = "#C2E7D9", "Medial Gut" = "#A8D5BA", "Distal Gut" = "#8EC39C",
    "Lungs" = "#D2B6E0", "Heart" = "#E8A4A4", "Liver" = "#F1C27D", "Pectoral" = "#AEC6CF"
  )
)

gaps_col <- cumsum(table(metadata$Metabolic_State))

#########################################
## 8. PLOT HEATMAP (ALL GENES)
#########################################
my_colors <- colorRampPalette(c("#b35806", "#f1a340","#fee0b6","#f7f7f7","#d8daeb","#998ec3","#542788"))(100)
my_breaks <- seq(-2.5, 2.5, length.out = 101)

tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_All_Genes_Heatmap.tiff")), 
     width = 8, height = 30, units = "in", res = 300, compression = "lzw")

# UPDATED: Replaced pheatmap with custom wrapper template
plot_custom_heatmap(counts_final_sorted,
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE,
                    annotation_col = annotation_col,
                    annotation_colors = ann_colors,
                    gaps_col = gaps_col,
                    scale = "none",
                    color = my_colors,
                    breaks = my_breaks,
                    main = paste0(FILE_PREFIX, " pathway: N vs D (Normalized within Tissues)"),
                    show_rownames = FALSE)
dev.off()

#########################################
## A. Filter for Top 10% Variable Genes
#########################################
gene_vars <- apply(counts_final_sorted, 1, var)
variance_cutoff <- quantile(gene_vars, 0.90)
top_10_genes <- names(gene_vars[gene_vars >= variance_cutoff])
counts_top_10pct <- counts_final_sorted[top_10_genes, ]

tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_Top_10_Percent_Heatmap.tiff")), 
     width = 8, height = 10, units = "in", res = 300, compression = "lzw")

# UPDATED: Replaced pheatmap with custom wrapper template
plot_custom_heatmap(counts_top_10pct,
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE, 
                    annotation_col = annotation_col,
                    annotation_colors = ann_colors,
                    gaps_col = gaps_col,
                    scale = "none",
                    color = my_colors,
                    breaks = my_breaks,
                    main = paste0("Top 10% Most Variable Genes ", FILE_PREFIX, " (N vs D)"),
                    show_rownames = TRUE)
dev.off()

############################################
### Variance Plots
############################################
samples_N <- metadata %>% filter(Metabolic_State == "N") %>% pull(SampleID)
samples_D <- metadata %>% filter(Metabolic_State == "D") %>% pull(SampleID)

valid_genes <- rownames(counts_final_sorted) 
matrix_N <- Normalizedcounts[valid_genes, samples_N]
matrix_D <- Normalizedcounts[valid_genes, samples_D]

var_N <- apply(matrix_N, 1, var)
var_D <- apply(matrix_D, 1, var)

var_df <- data.frame(Gene = names(var_N), Variance_N = var_N, Variance_D = var_D)

plot_scatter <- ggplot(var_df, aes(x = Variance_N, y = Variance_D)) +
  geom_point(alpha = 0.5, color = "#542788") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size=1) +
  theme_minimal() +
  labs(title = paste0("Gene Variance: Normothermy vs Torpor (D) - ", FILE_PREFIX),
       x = "Variance in Normothermy (N)",
       y = "Variance in Torpor (D)")

ggsave(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_Variance_Scatter.tiff")), 
       plot = plot_scatter, width = 6, height = 6, units = "in", dpi = 300, compression = "lzw")

# --- DEFINITION FOR AUDIT STEP ---
get_symbols_from_ids <- function(id_string) {
  if (is.na(id_string) || id_string == "") return(NA)
  ids <- unlist(strsplit(as.character(id_string), ","))
  syms <- tryCatch({
    mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  }, error = function(e) return(ids))
  syms[is.na(syms)] <- ids[is.na(syms)]
  paste(syms, collapse = ", ")
}

################################################################
## PATHVIEW & AUDIT (ITERATE BY TISSUE)
################################################################
TISSUE_TYPES <- levels(metadata$Tissue)
original_wd <- getwd()
setwd(genes_dir)

for (TISSUE_NAME in TISSUE_TYPES) {
  samps_N <- metadata %>% filter(Tissue == TISSUE_NAME, Metabolic_State == "N") %>% pull(SampleID)
  samps_D <- metadata %>% filter(Tissue == TISSUE_NAME, Metabolic_State == "D") %>% pull(SampleID)
  
  if(length(samps_N) == 0 || length(samps_D) == 0) next
  
  SUBSET_MATRIX_TISSUE <- counts_final_sorted[, c(samps_N, samps_D), drop=FALSE]
  mean_N <- rowMeans(SUBSET_MATRIX_TISSUE[, samps_N, drop=FALSE])
  mean_D <- rowMeans(SUBSET_MATRIX_TISSUE[, samps_D, drop=FALSE])
  diff_vector <- mean_D - mean_N
  
  entrez_map <- AnnotationDbi::select(org.Hs.eg.db, keys = names(diff_vector), columns = "ENTREZID", keytype = "SYMBOL")
  pathview_input_df <- merge(data.frame(SYMBOL = names(diff_vector), value = diff_vector), entrez_map, by="SYMBOL")
  pathview_input_df <- pathview_input_df[!is.na(pathview_input_df$ENTREZID), ]
  pathview_input_df <- aggregate(value ~ ENTREZID, data=pathview_input_df, mean)
  
  gene_data_vector <- pathview_input_df$value
  names(gene_data_vector) <- pathview_input_df$ENTREZID
  
  PV_SUFFIX <- paste0(FILE_PREFIX, "_kegg_native_", TISSUE_NAME)
  pv_out <- pathview(gene.data = gene_data_vector, pathway.id = TRANSLATION_HSA_ID, species = "hsa", 
                     out.suffix = PV_SUFFIX, limit = list(gene=2, cpd=1), low = "goldenrod2", 
                     mid = "gray", high = "purple", kegg.native = TRUE)
  
  # Redirect Audit Data CSV exports straight into main_fig_dir
  if (!is.null(pv_out$plot.data.gene)) {
    node_data <- pv_out$plot.data.gene
    bin_audit <- data.frame(Node_Label = node_data$labels, My_Mapped_Entrez = node_data$all.mapped, stringsAsFactors = FALSE)
    bin_audit$Tissue <- TISSUE_NAME
    bin_audit <- bin_audit[!is.na(bin_audit$My_Mapped_Entrez) & bin_audit$My_Mapped_Entrez != "", ]
    bin_audit$My_Genes_In_Bin <- sapply(bin_audit$My_Mapped_Entrez, get_symbols_from_ids)
    
    write.csv(bin_audit, file.path(main_fig_dir, paste0(FILE_PREFIX, "_Node_Audit_", TISSUE_NAME, ".csv")), row.names = FALSE)
  }
}
setwd(original_wd)

################################################################
## 8D. MOLECULAR PILLARS: THE TORPOR SHIFT (D - N) WITH RAW P-VALUE STARS
################################################################

# 1. Define the "High-Impact" Regulatory & Structural Gene List
selected_regulators <- c(
  # Global Initiation Gatekeepers
  "EIF4E", "EIF4G1", "EIF4G3", "EIF4A1", "EIF4B", "EIF3A", "EIF5",
  
  # The "Brakes" (Stress kinases that inhibit translation)
  "EIF2AK1", "EIF2AK2", "EIF2AK3", "EIF2AK4", "EIF4EBP1",
  
  # Cytosolic Ribosome
  "RPL7L1", "RPL29", "RPL34", "RPS7",
  
  # Mitochondrial Translation
  "MRPS2", "MRPS5", "MRPS18A", "MRPS28", "MRPL12", "MRPL22", "MRPL46", "MRPL47",
  
  # Ribosomal Kinases/Signaling
  "RPS6KA1", "RPS6KA3", "RPS6KB1"
)

# 2. Create the annotation dataframe
pillar_info <- data.frame(
  Gene = selected_regulators,
  Function = c(rep("Initiation Gatekeepers", 7),
               rep("Stress Response (Brakes)", 5),
               rep("Cytosolic Ribosome (Selected)", 4),
               rep("Mito-Translation Control", 8),
               rep("mTOR/Signaling Proxies", 3))
)

# 3. Filter for available genes
pillar_genes_available <- pillar_info$Gene[pillar_info$Gene %in% rownames(counts_final_sorted)]

# 4. Prepare Long Format Data for Calculation
counts_long <- as.data.frame(counts_final_sorted[pillar_genes_available, ]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

# 5. Calculate Mean Shifts (D - N)
diff_matrix <- counts_long %>%
  group_by(Gene, Tissue, Metabolic_State) %>%
  summarize(Mean_Z = mean(Z_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Metabolic_State, values_from = Mean_Z) %>%
  mutate(Shift = D - N) %>%
  dplyr::select(Gene, Tissue, Shift) %>%
  pivot_wider(names_from = Tissue, values_from = Shift) %>%
  column_to_rownames("Gene")

# Reorder columns to tissue order
diff_matrix <- diff_matrix[pillar_genes_available, tissue_order]

# 6. Calculate Stats using RAW p-values for Significance Stars
stats_summary <- counts_long %>%
  group_by(Gene, Tissue) %>%
  summarize(
    p_val = {
      n_N <- sum(Metabolic_State == "N", na.rm = TRUE)
      n_D <- sum(Metabolic_State == "D", na.rm = TRUE)
      if (n_N >= 2 && n_D >= 2) {
        tryCatch(t.test(Z_Score ~ Metabolic_State)$p.value, error = function(e) NA)
      } else { NA }
    }, .groups = 'drop'
  ) %>%
  # Evaluate RAW p-values directly for significance asterisks mapping
  mutate(sig_star = ifelse(!is.na(p_val) & p_val < 0.05, "*", ""))

# 7. Generate Significance Matrix
sig_matrix <- stats_summary %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")

# Align sig_matrix with diff_matrix
sig_matrix <- sig_matrix[rownames(diff_matrix), colnames(diff_matrix)]
sig_matrix[is.na(sig_matrix)] <- ""

# 8. Setup Visual Annotations
annotation_row <- data.frame(
  Function = factor(pillar_info$Function[pillar_info$Gene %in% pillar_genes_available], 
                    levels = unique(pillar_info$Function))
)
rownames(annotation_row) <- pillar_genes_available

pillar_colors <- list(
  Function = c(
    "Initiation Gatekeepers"        = "#2166ac", 
    "Stress Response (Brakes)"      = "#b2182b", 
    "Cytosolic Ribosome (Selected)" = "#1b7837", 
    "Mito-Translation Control"      = "#762a83", 
    "mTOR/Signaling Proxies"        = "#00441b"
  )
)

# 9. Plot Final Heatmap using raw p-value configuration straight to main figures target
tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_praw_Filtered_Regulators.tiff")), 
     width = 10, height = 13, units = "in", res = 300, compression = "lzw")

bold_rows <- lapply(rownames(diff_matrix), function(x) bquote(bold(.(x))))
bold_cols <- lapply(colnames(diff_matrix), function(x) bquote(bold(.(x))))

# UPDATED: Replaced pheatmap with custom wrapper template & mapped GLOBAL_FS_NUM dynamically
plot_custom_heatmap(diff_matrix,
                    display_numbers = sig_matrix,
                    fontsize_number = GLOBAL_FS_NUM,
                    labels_row = as.expression(bold_rows),
                    labels_col = as.expression(bold_cols),
                    cluster_rows = FALSE, 
                    cluster_cols = FALSE,
                    annotation_row = annotation_row,
                    annotation_colors = pillar_colors,
                    scale = "none",
                    color = my_colors,
                    breaks = seq(-2, 2, length.out = 101),
                    main = "Translation Control: Torpor Shift (D-N)\n(* = Raw p < 0.05, Purple = Up in Torpor)")
dev.off()

cat(paste0("Success! Translation raw p-value heatmap and audit CSV data files saved to: ", main_fig_dir, "\n"))