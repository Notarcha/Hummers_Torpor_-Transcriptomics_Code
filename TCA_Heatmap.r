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

#### TCA cycle ######
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

# --- STARTING ANALYSIS FOR TCA CYCLE PATHWAY ---
cat("--- STARTING ANALYSIS (N vs D ONLY) FOR TCA CYCLE (hsa00020) ---\n")

# Ensure the relative './genes' output directory exists
genes_dir <- here("genes")
if (!dir.exists(genes_dir)) {
  dir.create(genes_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created directory:", genes_dir, "\n")
}

FILE_PREFIX <- "TCA_CYCLE"
PATHWAY_ID <- "hsa00020"

#########################################
## 2. LOAD & CLEAN DATA (EXACT LOGIC)
#########################################

# Load Metadata & Counts
RNASeq_metadata <- read.csv(here("Results", "RNASeq_metadata.csv"), row.names = 1)
Normalizedcounts <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"), row.names = 1)

# Create Metadata Object
metadata <- RNASeq_metadata %>% rownames_to_column(var = "SampleID")

# Clean Names
metadata$SampleID <- metadata$SampleID %>% str_trim() %>% str_replace_all("-", "_")
colnames(Normalizedcounts) <- colnames(Normalizedcounts) %>% str_trim() %>% str_replace_all("-", "_")

# Intersect Samples
common_samples <- intersect(metadata$SampleID, colnames(Normalizedcounts))
metadata <- metadata %>% filter(SampleID %in% common_samples)
Normalizedcounts <- Normalizedcounts[, common_samples]

# --- FILTER STEP: REMOVE "T" ---
cat("Samples before filtering:", nrow(metadata), "\n")
metadata <- metadata %>% filter(Metabolic_State != "T")
Normalizedcounts <- Normalizedcounts[, metadata$SampleID]
cat("Samples after removing 'T':", nrow(metadata), "\n")

# *** RECODE TISSUE NAMES ***
metadata <- metadata %>%
  mutate(Tissue = recode(Tissue, 
                         "Gut1" = "Proximal Gut", 
                         "Gut2" = "Medial Gut", 
                         "Gut3" = "Distal Gut", 
                         "Pect" = "Pectoral"))

#########################################
## 3. DEFINE GENES & PILLARS
#########################################

# Define the TCA Cycle Pillars (User List + Critical Additions)
pillar_info <- data.frame(
  Gene = c(
    # 1. The Pyruvate Gate (Entry & Regulation)
    "MPC1", "MPC2",                 # The Door
    "PDHA1", "PDHB", "DLAT", "DLD", # The PDH Complex (The Converter)
    "PDK1", "PDK2", "PDK3", "PDK4", # The Regulators (The Lock)
    
    # 2. Citrate & Isocitrate Steps (Early Oxidation)
    "CS", "ACLY",                   # Entry & Exit (Lipogenesis)
    "ACO1", "ACO2",                 # Isomerization
    "IDH1", "IDH2", "IDH3A", "IDH3B", # Decarboxylation 1
    
    # 3. The Alpha-KG & Succinate Node (Energy Harvest)
    "OGDH", "OGDHL", "DLST",       # Alpha-KG Complex
    "SUCLG1", "SUCLG2", "SUCLA2",  # Succinyl-CoA Synthetase
    
    # 4. The Final Stretch (Succinate -> OAA)
    "SDHA", "SDHB", "SDHC", "SDHD", # Complex II
    "FH",                           # Fumarase
    "MDH1", "MDH2",                 # Malate Dehydrogenase
    
    # 5. Anaplerosis & Gluconeogenic Exit
    "PC",                           # Pyruvate Carboxylase (Refill OAA)
    "PCK1"                          # PEPCK (Drain OAA to Sugar)
  ),
  Function = c(
    rep("The Pyruvate Gate & Regulation", 10),
    rep("Citrate & Isocitrate Steps", 8),
    rep("The Alpha-KG & Succinate Node", 6),
    rep("Succinate -> OAA", 7),
    rep("Anaplerosis & Gluconeogenic Exit", 2)
  )
)

# Filter for genes actually in your dataset
genes_of_interest <- intersect(rownames(Normalizedcounts), pillar_info$Gene)

# *** TISSUE ORDER ***
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
## 5. PER-TISSUE Z-SCORE SCALING (EXACT LOGIC)
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

#########################################
## 6. REORDER MATRIX COLUMNS
#########################################
counts_final_sorted <- counts_per_tissue_scaled[, ordered_samples_state_first]

cat("Final Z-Score Matrix Ready:", dim(counts_final_sorted)[1], "genes x", dim(counts_final_sorted)[2], "samples.\n")

#########################################
## 7. PATHVIEW & AUDIT BY TISSUE
#########################################
get_symbols_from_ids <- function(id_string) {
  if (is.na(id_string) || id_string == "") return(NA)
  ids <- unlist(strsplit(as.character(id_string), ","))
  syms <- tryCatch({
    mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  }, error = function(e) return(ids))
  syms[is.na(syms)] <- ids[is.na(syms)]
  paste(syms, collapse = ", ")
}

original_wd <- getwd()
setwd(genes_dir)

for (TISSUE_NAME in levels(metadata$Tissue)) {
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
  
  pv_out <- pathview(gene.data = gene_data_vector, pathway.id = PATHWAY_ID, species = "hsa", 
                     out.suffix = paste0(FILE_PREFIX, "_", TISSUE_NAME), limit = list(gene=2, cpd=1), 
                     low = "goldenrod2", mid = "gray", high = "purple", kegg.native = TRUE)
  
  # Structural node tracking audit records exported to Figures
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
## 8D. MOLECULAR PILLARS: TCA SHIFT WITH RAW P-VALS
################################################################

# Filter for genes actually in your dataset
pillar_genes_available <- pillar_info$Gene[pillar_info$Gene %in% rownames(counts_final_sorted)]

# Transform Z-Scores to Long Format
counts_long <- as.data.frame(counts_final_sorted[pillar_genes_available, ]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata %>% dplyr::select(SampleID, Tissue, Metabolic_State), by = "SampleID")

# 1. Calculate Shift (Mean D Z-Score - Mean N Z-Score)
diff_matrix <- counts_long %>%
  group_by(Gene, Tissue, Metabolic_State) %>%
  summarize(Mean_Z = mean(Z_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Metabolic_State, values_from = Mean_Z) %>%
  mutate(Shift = D - N) %>%
  dplyr::select(Gene, Tissue, Shift) %>%
  pivot_wider(names_from = Tissue, values_from = Shift) %>%
  column_to_rownames("Gene")

# Ensure column order matches tissue order
diff_matrix <- diff_matrix[pillar_genes_available, tissue_order]

# 2. Calculate P-Values using RAW values for significance checks
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
  # Apply raw significance threshold filtering directly
  mutate(sig_star = ifelse(!is.na(p_val) & p_val < 0.05, "*", ""))

# 3. Create Significance Matrix for Plotting
sig_matrix <- stats_summary %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")

# Re-align sig_matrix with diff_matrix
sig_matrix <- sig_matrix[rownames(diff_matrix), colnames(diff_matrix)]
sig_matrix[is.na(sig_matrix)] <- ""

# 4. Prepare Annotations
annotation_row <- data.frame(
  Function = factor(pillar_info$Function[pillar_info$Gene %in% pillar_genes_available], 
                    levels = unique(pillar_info$Function))
)
rownames(annotation_row) <- pillar_genes_available

pillar_colors <- list(
  Function = c(
    "The Pyruvate Gate & Regulation"       = "#543005", 
    "Citrate & Isocitrate Steps"           = "#8c510a", 
    "The Alpha-KG & Succinate Node"        = "#01665e", 
    "Succinate -> OAA"                     = "#35978f", 
    "Anaplerosis & Gluconeogenic Exit"     = "#1a1a1a"  
  )
)

my_colors <- colorRampPalette(c("#b35806", "#f1a340","#fee0b6","#f7f7f7","#d8daeb","#998ec3","#542788"))(100)

# 5. Plot Heatmap to global figures output folder
tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_praw_Heatmap_TCA.tiff")), 
     width = 10, height = 12, units = "in", res = 300, compression = "lzw")

bold_rows <- lapply(rownames(diff_matrix), function(x) bquote(bold(.(x))))
bold_cols <- lapply(colnames(diff_matrix), function(x) bquote(bold(.(x))))

# UPDATED: Replaced pheatmap with custom wrapper template & mapped GLOBAL_FS_NUM dynamically
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
                    main = "TCA Cycle: Torpor Shift (D-N)\n(* = Raw p < 0.05, Purple = Up in Torpor)",
                    gaps_row = cumsum(table(annotation_row$Function)[unique(annotation_row$Function)]))
dev.off()

cat(paste0("Success! TCA Cycle raw p-value heatmap and audit data documents saved to: ", main_fig_dir, "\n"))