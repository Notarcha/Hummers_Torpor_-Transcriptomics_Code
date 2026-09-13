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

####### Transcription machinery ######

## LOAD LIBRARIES
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

# Ensure the relative './genes' output directory exists
genes_dir <- here("genes")
if (!dir.exists(genes_dir)) {
  dir.create(genes_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created directory:", genes_dir, "\n")
}

cat("--- STARTING ANALYSIS (N vs D ONLY) FOR BASAL TRANSCRIPTION PATHWAY (hsa03022) ---\n")

# Define Transcription Pathway ID
TRANSCRIPTION_PATHWAY_ID <- "path:hsa03022"
TRANSCRIPTION_HSA_ID <- "hsa03022"
FILE_PREFIX <- "Transcription"


## 1. PREPARE KEGG GENE LIST (TRANSCRIPTION)
pathway_id <- TRANSCRIPTION_PATHWAY_ID 
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

cat("Step 1: Mapped", nrow(kegg_pathway_genes), "Basal Transcription pathway genes.\n")


## 2. LOAD & CLEAN DATA
# Load Metadata & Counts using project-relative paths
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

# --- NEW FILTER STEP: REMOVE "T" ---
cat("Samples before filtering:", nrow(metadata), "\n")
metadata <- metadata %>% filter(Metabolic_State != "T")
Normalizedcounts <- Normalizedcounts[, metadata$SampleID]
cat("Samples after removing 'T':", nrow(metadata), "\n")

# RECODE TISSUE NAMES
metadata <- metadata %>%
  mutate(Tissue = recode(Tissue, 
                         "Gut1" = "Proximal Gut", 
                         "Gut2" = "Medial Gut", 
                         "Gut3" = "Distal Gut", 
                         "Pect" = "Pectoral"))

## 3. DEFINE GENES & ORDER (UPDATED FILTER)
# Define custom pillar genes early so they survive the KEGG filter
pol1_genes <- c("POLR1A", "POLR1B", "POLR1C", "POLR1D", "POLR1E")
pol2_genes <- c("POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E", "POLR2F", 
                "POLR2G", "POLR2H", "POLR2I", "POLR2J", "POLR2K", "POLR2L")
pol3_genes <- c("POLR3A", "POLR3B", "POLR3C", "POLR3D", "POLR3E", "POLR3F", 
                "POLR3G", "POLR3H", "POLR3K", "POLR3GL")
gtf_genes <- c("TBP", "GTF2A1", "GTF2A2", "GTF2B", "GTF2E1", "GTF2E2", 
               "GTF2F1", "GTF2F2", "GTF2H1", "GTF2H2", "GTF2H3", "GTF2H4", "GTF2H5")
mediator_genes <- c("MED1", "MED4", "MED6", "MED7", "MED8", "MED9", "MED10", 
                    "MED11", "MED12", "MED13", "MED14", "MED15", "MED16", "MED17")

all_custom_pillars <- c(pol1_genes, pol2_genes, pol3_genes, gtf_genes, mediator_genes)

# Combine KEGG genes with Custom genes, then intersect with your matrix
combined_target_genes <- unique(c(kegg_pathway_genes$symbol, all_custom_pillars))
genes_of_interest <- intersect(rownames(Normalizedcounts), combined_target_genes)

# UPDATED TISSUE ORDER NAMES
tissue_order <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Lungs", "Heart", "Liver", "Pectoral")
state_order  <- c("N", "D")

metadata <- metadata %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_order),
    Metabolic_State = factor(Metabolic_State, levels = state_order)
  )

## 4. SORT METADATA
metadata <- metadata %>% arrange(Metabolic_State, Tissue)
ordered_samples_state_first <- metadata$SampleID


## 5. PER-TISSUE Z-SCORE SCALING
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


## 6. REORDER MATRIX COLUMNS
counts_final_sorted <- counts_per_tissue_scaled[, ordered_samples_state_first]

cat("Final Matrix Ready:", dim(counts_final_sorted)[1], "genes x", dim(counts_final_sorted)[2], "samples.\n")


## 7. PREPARE HEATMAP ANNOTATION
annotation_col <- metadata %>%
  dplyr::select(SampleID, Metabolic_State, Tissue) %>%
  column_to_rownames("SampleID")

# UPDATED ANNOTATION COLORS NAMES
ann_colors <- list(
  Metabolic_State = c(N = "gray", D = "gray30"),
  Tissue = c(
    "Proximal Gut" = "#C2E7D9", "Medial Gut" = "#A8D5BA", "Distal Gut" = "#8EC39C",
    "Lungs" = "#D2B6E0", "Heart" = "#E8A4A4", "Liver" = "#F1C27D", "Pectoral" = "#AEC6CF"
  )
)

gaps_col <- cumsum(table(metadata$Metabolic_State))

## 8. PLOT HEATMAP (ALL GENES)
my_colors <- colorRampPalette(c("#b35806", "#f1a340","#fee0b6","#f7f7f7","#d8daeb","#998ec3","#542788"))(100)
my_breaks <- seq(-2.5, 2.5, length.out = 101)

tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_All_Genes_Heatmap.tiff")), 
     width = 8, height = 30, units = "in", res = 300, compression = "lzw")

# FIXED: cluster_rows = FALSE to prevent C stack overflow on massive matrices
plot_custom_heatmap(counts_final_sorted,
                    cluster_rows = FALSE, 
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


## A. Filter for Top 10% Variable Genes
gene_vars <- apply(counts_final_sorted, 1, var)
variance_cutoff <- quantile(gene_vars, 0.90)
top_10_genes <- names(gene_vars[gene_vars >= variance_cutoff])
counts_top_10pct <- counts_final_sorted[top_10_genes, ]

tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_Top_10_Percent_Heatmap.tiff")), 
     width = 8, height = 10, units = "in", res = 300, compression = "lzw")

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


### Is Variance Lower in D when compared to N 
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

var_long <- var_df %>%
  pivot_longer(cols = c("Variance_N", "Variance_D"), names_to = "State", values_to = "Variance")

plot_box <- ggplot(var_long, aes(x = State, y = Variance, fill = State)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + 
  geom_jitter()+coord_cartesian(ylim = quantile(var_long$Variance, c(0, 0.95))) +
  scale_fill_manual(values = c("Variance_N" = "#998ec3", "Variance_D" = "#f1a340")) +
  theme_minimal() +
  labs(title = paste0("Distribution of Variance (Zoomed to 95%)_", FILE_PREFIX), y = "Variance (rlog)")

ggsave(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_Variance_Boxplot.tiff")), 
       plot = plot_box, width = 5, height = 7, units = "in", dpi = 300, compression = "lzw")


## 9. BOXPLOT OF TOP 10% GENES BY STATE
top_genes_df <- as.data.frame(counts_top_10pct) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(metadata %>% dplyr::select(SampleID, Metabolic_State), by = "SampleID")

tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_Top_10_PercentVariableGenes_Boxplot.tiff")), 
     width = 8, height = 10, units = "in", res = 300, compression = "lzw")

print(
  ggplot(top_genes_df, aes(x = Metabolic_State, y = Z_Score, fill = Metabolic_State)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
    facet_wrap(~ Gene, scales = "free_y") +
    scale_fill_manual(values = c("N" = "#998ec3", "D" = "#f1a340")) + 
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 10), legend.position = "none") +
    labs(title = paste0("Expression of Top 10% Most Variable ", FILE_PREFIX, " Genes"),
         subtitle = "Comparison of Normalized Z-Scores (N vs D)", y = "Z-Score Expression", x = "Physiological State")
)
dev.off()

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

## FULL WORKFLOW: PATHVIEW & AUDIT (ITERATE BY TISSUE)
TISSUE_TYPES <- levels(metadata$Tissue)
original_wd <- getwd()
setwd(genes_dir)

for (TISSUE_NAME in TISSUE_TYPES) {
  cat("\n--- STARTING ANALYSIS FOR TISSUE:", TISSUE_NAME, "---\n")
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
  pv_out <- pathview(gene.data = gene_data_vector, pathway.id = TRANSCRIPTION_HSA_ID, species = "hsa", 
                     out.suffix = PV_SUFFIX, limit = list(gene=2, cpd=1), low = "goldenrod2", 
                     mid = "gray", high = "purple", kegg.native = TRUE)
  
  # Audit
  if (!is.null(pv_out$plot.data.gene)) {
    node_data <- pv_out$plot.data.gene
    bin_audit <- data.frame(Node_Label = node_data$labels, My_Mapped_Entrez = node_data$all.mapped, stringsAsFactors = FALSE)
    bin_audit$Tissue <- TISSUE_NAME
    bin_audit <- bin_audit[!is.na(bin_audit$My_Mapped_Entrez) & bin_audit$My_Mapped_Entrez != "", ]
    bin_audit$My_Genes_In_Bin <- sapply(bin_audit$My_Mapped_Entrez, get_symbols_from_ids)
    write.csv(bin_audit, paste0(FILE_PREFIX, "_Node_Audit_", TISSUE_NAME, ".csv"), row.names = FALSE)
  }
}
setwd(original_wd)

## 8D. MOLECULAR PILLARS: THE TORPOR SHIFT (D - N) WITH RAW P-VALUE STARS
pillar_info <- data.frame(
  Gene = c(pol1_genes, pol2_genes, pol3_genes, gtf_genes, mediator_genes),
  Function = c(rep("RNA Polymerase I (rRNA)", length(pol1_genes)),
               rep("RNA Polymerase II (mRNA)", length(pol2_genes)),
               rep("RNA Polymerase III (tRNA)", length(pol3_genes)),
               rep("General Transcription Factors (GTFs)", length(gtf_genes)),
               rep("Mediator Complex", length(mediator_genes)))
)

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

# 3. Calculate Stats using RAW p-values for Significance Stars
stats_summary <- counts_long %>%
  group_by(Gene, Tissue) %>%
  summarize(
    p_val = {
      n_N <- sum(Metabolic_State == "N", na.rm = TRUE)
      n_D <- sum(Metabolic_State == "D", na.rm = TRUE)
      if (n_N >= 2 && n_D >= 2) {
        # Welch's t-test
        tryCatch(t.test(Z_Score ~ Metabolic_State)$p.value, error = function(e) NA)
      } else {
        NA
      }
    },
    .groups = 'drop'
  ) %>%
  mutate(sig_star = ifelse(!is.na(p_val) & p_val < 0.05, "*", ""))

# 4. Generate Significance Matrix
sig_matrix <- stats_summary %>%
  dplyr::select(Gene, Tissue, sig_star) %>%
  pivot_wider(names_from = Tissue, values_from = sig_star) %>%
  column_to_rownames("Gene")

sig_matrix <- sig_matrix[rownames(diff_matrix), colnames(diff_matrix)]
sig_matrix[is.na(sig_matrix)] <- ""

# 5. Annotations and Colors
annotation_row <- data.frame(
  Function = factor(pillar_info$Function[pillar_info$Gene %in% pillar_genes_available], 
                    levels = unique(pillar_info$Function))
)
rownames(annotation_row) <- pillar_genes_available

pillar_colors <- list(
  Function = c(
    "RNA Polymerase I (rRNA)"                = "#543005", 
    "RNA Polymerase II (mRNA)"               = "#8c510a", 
    "RNA Polymerase III (tRNA)"              = "#bf812d", 
    "General Transcription Factors (GTFs)" = "#01665e", 
    "Mediator Complex"                       = "#35978f"
  )
)

# 6. Plot Final Heatmap directly to global Figures location
tiff(filename = file.path(main_fig_dir, paste0(FILE_PREFIX, "_praw_Heatmap_Transcription.tiff")), 
     width = 10, height = 24, units = "in", res = 300, compression = "lzw")

bold_rows <- lapply(rownames(diff_matrix), function(x) bquote(bold(.(x))))
bold_cols <- lapply(colnames(diff_matrix), function(x) bquote(bold(.(x))))

plot_custom_heatmap(diff_matrix,
                    display_numbers = sig_matrix,          # <--- OVERLAYS THE RAW P-VALUE ASTERISKS
                    fontsize_number = GLOBAL_FS_NUM,       # <--- USES GLOBAL SIZE FOR ASTERISKS
                    labels_row = as.expression(bold_rows),
                    labels_col = as.expression(bold_cols),
                    cluster_rows = FALSE, cluster_cols = FALSE,
                    annotation_row = annotation_row,
                    annotation_colors = pillar_colors,
                    scale = "none",
                    color = my_colors,
                    breaks = seq(-2, 2, length.out = 101),
                    main = "Basal Transcription: Torpor Shift (D-N)\n(* = Raw p < 0.05, Purple = Up in Torpor)",
                    gaps_row = cumsum(table(annotation_row$Function)[unique(annotation_row$Function)]))
dev.off()
cat(paste0("Success! Basal Transcription heatmap with raw p-value stars saved to: ", main_fig_dir, "\n"))