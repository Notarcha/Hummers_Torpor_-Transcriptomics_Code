# Loading Libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(patchwork)
library(apeglm)
library(ashr)
library(tibble)
library(tidyr)
library(here)
library(RColorBrewer)
library(UpSetR)
library(purrr)
library(ggVennDiagram)

# Set main output directory relative to project root
main_fig_dir <- here("Outputs")
dir.create(main_fig_dir, showWarnings = FALSE, recursive = TRUE)

############### Read in Counts Data and Metadata ###########
data <- read.csv(here("Results", "data.csv"), row.names = 1)
metadata <- read.csv(here("Results", "RNASeq_metadata.csv"), row.names = 1)

metadata$Tissue <- recode(
  metadata$Tissue,
  "Gut1" = "Proximal Gut",
  "Gut2" = "Medial Gut",
  "Gut3" = "Distal Gut",
  "Pect" = "Pectoral"
)

# Enforce requested explicit factor ordering for the PCA legend
metadata$Tissue <- factor(
  metadata$Tissue,
  levels = c("Proximal Gut", "Medial Gut", "Distal Gut", "Lungs", "Heart", "Liver", "Pectoral")
)

metadata$Metabolic_State <- factor(
  metadata$Metabolic_State,
  levels = c("N", "T", "D")
)

data <- data[, row.names(metadata)]
stopifnot(all(row.names(metadata) == colnames(data)))

data <- data %>% filter(rowSums(.) >= 10)
counts.matrix <- as.matrix(data)

#### PCA to Identify Outliers #####
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(data), 
  colData   = metadata, 
  design    = ~ Tissue + Metabolic_State
)
rld <- vst(dds, blind = FALSE)

# Manual PCA Generation (for PC3 access)
pca_res <- prcomp(t(assay(rld)))
percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2))

# Create plotting dataframe
pcadata <- as.data.frame(pca_res$x)
pcadata$Tissue <- metadata$Tissue
pcadata$Metabolic_State <- metadata$Metabolic_State
pcadata$SampleID <- rownames(pcadata)

# Global Settings
tissue_colors <- c(
  "Proximal Gut" = "#b35806",
  "Medial Gut"   = "#f1a340",
  "Distal Gut"   = "#fee0b6",
  "Lungs"        = "#67a9cf",
  "Heart"        = "#5ab4ac",
  "Liver"        = "#998ec3",
  "Pectoral"     = "#542788"
)
state_shapes <- c("N" = 16, "T" = 17, "D" = 15)

# Shared labeling settings
repel_settings <- list(
  aes(label = row.names(pcadata)), 
  size = 2.5,
  max.overlaps = Inf,
  box.padding = 0.25,
  point.padding = 0.2,
  min.segment.length = 2, 
  segment.alpha = 0.3,
  force = 0.3,             
  bg.color = "white",
  bg.r = 0.1
)

# Create Individual Plots
p1 <- ggplot(pcadata, aes(PC1, PC2)) +
  geom_point(aes(color = Tissue, shape = Metabolic_State), size = 3, stroke = 1) +
  do.call(geom_text_repel, repel_settings) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("PC1 vs PC2") +
  scale_color_manual(values = tissue_colors) +
  scale_shape_manual(values = state_shapes) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") 

p2 <- ggplot(pcadata, aes(PC1, PC3)) +
  geom_point(aes(color = Tissue, shape = Metabolic_State), size = 3, stroke = 1) +
  do.call(geom_text_repel, repel_settings) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC3: ", percentVar[3], "%")) +
  ggtitle("PC1 vs PC3") +
  scale_color_manual(values = tissue_colors) +
  scale_shape_manual(values = state_shapes) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

p3 <- ggplot(pcadata, aes(PC2, PC3)) +
  geom_point(aes(color = Tissue, shape = Metabolic_State), size = 3, stroke = 1) +
  do.call(geom_text_repel, repel_settings) +
  xlab(paste0("PC2: ", percentVar[2], "%")) +
  ylab(paste0("PC3: ", percentVar[3], "%")) +
  ggtitle("PC2 vs PC3") +
  scale_color_manual(values = tissue_colors) +
  scale_shape_manual(values = state_shapes) +
  theme_classic(base_size = 12)

# Stitch together using Patchwork
final_plot <- (p1 | p2 | p3) + 
  plot_layout(guides = 'collect') + 
  plot_annotation(
    title = 'Multi-Component PCA Analysis',
    theme = theme(plot.title = element_text(size = 18, hjust = 0.5))
  )

ggsave(file.path(main_fig_dir, "PCA_tissue_type.tiff"), plot = final_plot, width = 14, height = 10, dpi = 300, compression = "lzw")

#### Pairwise Euclidean Distance Heatmap ##### 
graphics.off()

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

sampleDists <- dist(t(assay(rld))) 
sampleDistMatrix <- as.matrix(sampleDists)

# Create Tissue Labels
new_labels <- paste(rownames(metadata), metadata$Tissue, sep = " - ")
dist_matrix_labeled <- sampleDistMatrix
rownames(dist_matrix_labeled) <- new_labels
colnames(dist_matrix_labeled) <- new_labels

# Generate Heatmap
tiff(file.path(main_fig_dir, "Sample_Distances_Pre_Outliers.tiff"), 
     width = 18, height = 18, units = "in", res = 300, compression = "lzw")
pheatmap(
  dist_matrix_labeled,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors,
  main = "Sample-to-Sample Distances by Tissue"
)
dev.off()

### Outlier Dropping #### 
outliers <- c("AS1", "AS53", "AS6")
metadata_clean <- metadata[!(rownames(metadata) %in% outliers), ]
data_clean <- data[, rownames(metadata_clean)]

## Grouped DESeq2 Design
metadata_clean$Group <- factor(
  paste0(metadata_clean$Tissue, "_", metadata_clean$Metabolic_State)
)

dds_grouped <- DESeqDataSetFromMatrix(
  countData = as.matrix(data_clean),
  colData   = metadata_clean,
  design    = ~ Group
)

dds_grouped <- dds_grouped[rowSums(counts(dds_grouped)) >= 10, ]
dds_grouped <- DESeq(dds_grouped)

######## Defining DEG Threshold ####################
padj.cutoff <- 0.05
lfc.cutoff  <- 0.58

##### DEG Analysis #####
tissues <- unique(metadata_clean$Tissue)
all_tissue_results <- list()
deg_summary_list <- list()

for (t in tissues) {
  
  treat <- paste0(t, "_D")
  ref   <- paste0(t, "_N")
  
  if (treat %in% levels(metadata_clean$Group) &&
      ref   %in% levels(metadata_clean$Group)) {
    
    message("Processing tissue: ", t)
    
    res_raw <- results(
      dds_grouped,
      contrast = c("Group", treat, ref),
      alpha = padj.cutoff
    )
    
    res_shrunk <- lfcShrink(
      dds_grouped,
      contrast = c("Group", treat, ref),
      res = res_raw,
      type = "ashr"
    )
    
    res_df <- res_shrunk %>%
      data.frame() %>%
      rownames_to_column("gene") %>%
      as_tibble() %>%
      mutate(
        diffexpressed = case_when(
          log2FoldChange >  lfc.cutoff & padj < padj.cutoff ~ "UP",
          log2FoldChange < -lfc.cutoff & padj < padj.cutoff ~ "DOWN",
          TRUE ~ "NO"
        )
      )
    
    file_name <- paste0(t, "_D_vs_N.csv")
    out_path  <- file.path(main_fig_dir, file_name)
    
    write.csv(res_df, file = out_path, row.names = FALSE)
    
    all_tissue_results[[t]] <- res_df
    
    deg_summary_list[[t]] <- tibble(
      Tissue        = t,
      Downregulated = sum(res_df$diffexpressed == "DOWN", na.rm = TRUE),
      Upregulated   = sum(res_df$diffexpressed == "UP",   na.rm = TRUE),
      Total_DEGs    = sum(res_df$diffexpressed != "NO",   na.rm = TRUE)
    )
  }
}

## Final DEG Summary Table
deg_summary <- bind_rows(deg_summary_list) %>%
  mutate(Down_minus_Up = Downregulated - Upregulated) %>%
  arrange(desc(Total_DEGs))

print(deg_summary)

write.csv(deg_summary, file.path(main_fig_dir, "DEG_summary_table.csv"), row.names = FALSE)

## Save Normalized Counts
normalized_counts <- counts(dds_grouped, normalized = TRUE)
write.csv(
  normalized_counts,
  file.path(main_fig_dir, "normalized_counts_STAR_master.csv"),
  row.names = TRUE
)

## Boxplots for Top 30 UP/DOWN Genes Per Tissue
metadata_clean$State <- recode(
  metadata_clean$Metabolic_State,
  "N" = "Normothermy",
  "T" = "Transition",
  "D" = "Deep Torpor"
)

metadata_clean$State <- factor(
  metadata_clean$State,
  levels = c("Normothermy", "Transition", "Deep Torpor")
)

state_colors <- c(
  "Normothermy" = "goldenrod",
  "Transition"  = "skyblue",
  "Deep Torpor" = "purple"
)

norm_counts <- counts(dds_grouped, normalized = TRUE)
norm_counts <- norm_counts[, rownames(metadata_clean)]

plot_top_gene_boxplots <- function(tissue, res_df, norm_counts, metadata_clean) {
  
  samples_tissue <- rownames(metadata_clean)[metadata_clean$Tissue == tissue]
  meta_tissue <- metadata_clean[samples_tissue, ]
  
  top_up <- res_df %>%
    filter(diffexpressed == "UP") %>%
    arrange(padj) %>%
    slice_head(n = 30) %>%
    pull(gene)
  
  top_down <- res_df %>%
    filter(diffexpressed == "DOWN") %>%
    arrange(padj) %>%
    slice_head(n = 30) %>%
    pull(gene)
  
  top_genes <- unique(c(top_up, top_down))
  
  if (length(top_genes) == 0) return(NULL)
  
  df_long <- norm_counts[top_genes, samples_tissue, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "SampleID",
      values_to = "Expression"
    ) %>%
    left_join(
      meta_tissue %>% rownames_to_column("SampleID"),
      by = "SampleID"
    )
  
  ggplot(df_long, aes(x = State, y = Expression, fill = State)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
    facet_wrap(~ gene, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = state_colors) +
    theme_classic(base_size = 12) +
    theme(
      strip.text = element_text(size = 7),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = paste("Top 30 UP & DOWN Genes:", tissue),
      x = NULL,
      y = "Normalized expression"
    )
}

for (t in tissues) {
  message("Plotting tissue: ", t)
  
  p <- plot_top_gene_boxplots(
    tissue = t,
    res_df = all_tissue_results[[t]],
    norm_counts = norm_counts,
    metadata_clean = metadata_clean
  )
  
  if (!is.null(p)) {
    ggsave(
      file.path(main_fig_dir, paste0("Boxplots_Top30_UP_DOWN_D_Vs_N_", t, ".tiff")),
      plot = p,
      width = 16,
      height = 12,
      dpi = 300,
      compression = "lzw"
    )
  }
}

## UpSet Plot of DEGs
deg_list_clean <- all_tissue_results %>%
  map(~ .x %>%
        filter(diffexpressed != "NO") %>%   
        pull(gene) %>%                       
        as.character() %>%                  
        unique()
  )

upset_data <- UpSetR::fromList(deg_list_clean)

tiff(file.path(main_fig_dir, "UpSet_DEGs.tiff"), width = 10, height = 7, units = "in", res = 300, compression = "lzw")
UpSetR::upset(
  upset_data,
  sets = names(deg_list_clean),
  nsets = length(deg_list_clean),
  order.by = "freq",
  decreasing = TRUE,
  point.size = 3.5,
  line.size = 1,
  main.bar.color = "midnightblue",
  sets.bar.color = "goldenrod3",
  text.scale = c(1.6, 1.3, 1.6, 1.3, 1.6, 1.3)
)
dev.off()

###### Shared Genes #########
gene_counts <- deg_list_clean %>%
  enframe(name = "Tissue", value = "Gene") %>%  
  unnest(cols = Gene) %>%                       
  group_by(Gene) %>%                            
  summarize(
    Count = n(),                                
    Shared_By = paste(sort(unique(Tissue)), collapse = ", "),  
    .groups = "drop"
  ) %>%
  filter(Count >= 3) %>%                        
  arrange(desc(Count))

write.csv(
  gene_counts,
  file.path(main_fig_dir, "Shared_DEGs_3plus_Tissues_D_vs_N.csv"),
  row.names = FALSE
)

p_venn <- ggVennDiagram(
  deg_list_clean,
  label = "count",
  label_alpha = 0
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )

ggsave(file.path(main_fig_dir, "VennDiagram_DEGs.tiff"), plot = p_venn, width = 8, height = 6, dpi = 300, compression = "lzw")

################ Volcano Plots ################
for (t in names(all_tissue_results)) {
  
  res_table_dn_df <- all_tissue_results[[t]]
  
  if (!is.null(res_table_dn_df) && nrow(res_table_dn_df) > 0) {
    message("Generating Volcano Plot for: ", t)
    
    top_up_fc <- res_table_dn_df %>% 
      filter(diffexpressed == "UP") %>% 
      arrange(desc(log2FoldChange)) %>% 
      head(20)
    
    top_down_fc <- res_table_dn_df %>% 
      filter(diffexpressed == "DOWN") %>% 
      arrange(log2FoldChange) %>% 
      head(20)
    
    plot_labels <- rbind(top_up_fc, top_down_fc)
    
    max_lfc <- max(abs(res_table_dn_df$log2FoldChange), na.rm = TRUE)
    x_limit <- max_lfc * 1.1 
    
    cols_points <- c("DOWN" = "goldenrod3", "NO" = "grey30", "UP" = "mediumpurple1")
    cols_text   <- c("DOWN" = "goldenrod4", "NO" = "grey30",  "UP" = "purple3")
    
    p_volcano <- ggplot(data = res_table_dn_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(
        aes(fill = diffexpressed), 
        shape = 21, 
        size = 2.5, 
        color = "transparent", 
        alpha = 0.7
      ) + 
      scale_fill_manual(values = cols_points) +
      xlim(-x_limit, x_limit) +
      theme_minimal() +
      theme(
        text = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 16),
        axis.text  = element_text(color = "black", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        legend.position = "none"
      ) +
      geom_vline(xintercept = c(-0.58, 0.58), col = "grey30", linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), col = "grey30", linetype = "dashed") +
      labs(
        title = paste("Volcano Plot:", t),
        subtitle = "Top 20 Up/Down Ranked by Fold Change",
        x = "log2 Fold Change",
        y = "-log10 Adjusted P-value"
      ) +
      geom_text_repel(
        data = plot_labels,
        aes(label = gene, color = diffexpressed), 
        size = 3.5,                  
        force = 10,                  
        max.overlaps = Inf,           
        segment.color = 'grey50',
        segment.alpha = 0.6,
        box.padding = 0.5,
        point.padding = 0.5,
        fontface = "italic",
        min.segment.length = 0
      ) +
      scale_color_manual(values = cols_text)
    
    safe_name <- gsub(" ", "_", t)
    save_path <- file.path(main_fig_dir, paste0("Volcano_Dynamic_D_vs_N_", safe_name, ".tiff"))
    
    ggsave(save_path, plot = p_volcano, width = 8, height = 7, dpi = 300, compression = "lzw")
  }
}

################ Torpor – Normothermy Heatmap ################
genes_to_plot <- gene_counts$Gene

norm_counts_df <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% genes_to_plot)

sample_metadata <- metadata_clean %>%
  rownames_to_column("SampleID") %>%
  select(SampleID, Tissue, Metabolic_State)

norm_long <- norm_counts_df %>%
  pivot_longer(
    cols = -gene,
    names_to = "SampleID",
    values_to = "expression"
  ) %>%
  left_join(sample_metadata, by = "SampleID")

gene_tissue_state_means <- norm_long %>%
  filter(Metabolic_State %in% c("D", "N")) %>%
  group_by(gene, Tissue, Metabolic_State) %>%
  summarize(
    mean_expression = mean(expression),
    .groups = "drop"
  )

torpor_vs_norm <- gene_tissue_state_means %>%
  pivot_wider(
    names_from = Metabolic_State,
    values_from = mean_expression
  ) %>%
  mutate(
    log2FC_Torpor_vs_Norm = log2(D + 1) - log2(N + 1)
  )

heatmap_mat <- torpor_vs_norm %>%
  select(gene, Tissue, log2FC_Torpor_vs_Norm) %>%
  pivot_wider(
    names_from = Tissue,
    values_from = log2FC_Torpor_vs_Norm
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

tissue_order <- c(
  "Proximal Gut", "Medial Gut", "Distal Gut",
  "Liver", "Heart", "Lungs", "Pectoral"
)
heatmap_mat <- heatmap_mat[, tissue_order, drop = FALSE]
heatmap_mat <- heatmap_mat[complete.cases(heatmap_mat), ]

heat_colors <- colorRampPalette(c("goldenrod3", "white", "purple"))(255)

tiff(file.path(main_fig_dir, "Heatmap_Torpor_vs_Normothermy_TissueDifferences.tiff"), 
     width = 10, height = 12, units = "in", res = 300, compression = "lzw")
pheatmap(
  heatmap_mat,
  color = heat_colors,
  breaks = seq(-2, 2, length.out = 256),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_cols = FALSE,
  scale = "none",
  border_color = NA,
  fontsize_row = 12,
  fontsize_col = 12,
  main = "Torpor – Normothermy Expression Differences by Tissue"
)
dev.off()

################ Sample-Level Heatmap ################
genes_to_plot <- gene_counts$Gene

norm_counts_df <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% genes_to_plot)

sample_metadata <- metadata_clean %>%
  rownames_to_column("SampleID") %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = c(
        "Proximal Gut", "Medial Gut", "Distal Gut",
        "Heart", "Lungs", "Liver", "Pectoral"
      )
    ),
    Metabolic_State = recode(
      metadata_clean$Metabolic_State,
      "N" = "Normothermy",
      "T" = "Transition",
      "D" = "Deep Torpor"
    ),
    Metabolic_State = factor(
      Metabolic_State,
      levels = c("Normothermy", "Transition", "Deep Torpor")
    )
  ) %>%
  select(SampleID, Tissue, Metabolic_State)

heatmap_mat_samples <- norm_counts_df %>%
  column_to_rownames("gene") %>%
  as.matrix()

sample_order <- sample_metadata %>%
  arrange(Metabolic_State, Tissue) %>%
  pull(SampleID)

heatmap_mat_samples <- heatmap_mat_samples[, sample_order]

annotation_col <- sample_metadata %>%
  column_to_rownames("SampleID") %>%
  select(Tissue, Metabolic_State)

heat_colors_samples <- colorRampPalette(c("goldenrod", "white", "purple"))(255)

ann_colors <- list(
  Tissue = c(
    "Proximal Gut" = "#b35806",
    "Medial Gut"   = "#f1a340",
    "Distal Gut"   = "#fee0b6",
    "Heart"        = "#5ab4ac",
    "Lungs"        = "#67a9cf",
    "Liver"        = "#998ec3",
    "Pectoral"     = "#542788"
  ),
  Metabolic_State = c(
    "Normothermy" = "lightgray",
    "Transition"  = "gray60",
    "Deep Torpor" = "black"
  )
)

tiff(file.path(main_fig_dir, "Heatmap_Sample_Level_Shared_DEGs.tiff"), 
     width = 12, height = 12, units = "in", res = 300, compression = "lzw")
pheatmap(
  heatmap_mat_samples,
  color = heat_colors_samples,
  breaks = seq(-2, 2, length.out = 256),
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  scale = "row",          
  cluster_cols = FALSE,  
  border_color = NA,
  fontsize_row = 12,
  fontsize_col = 8,
  main = "Sample-level Expression of Shared DEGs",
  annotation_col = annotation_col,
  annotation_colors = ann_colors
)
dev.off()

###### Get DEG Panels #####
tissues_list <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Liver", "Heart", "Lungs", "Pectoral")
file_names <- paste0(tissues_list, "_D_vs_N.csv")
file_paths <- file.path(main_fig_dir, file_names)
names(file_paths) <- tissues_list

master_deg_table <- purrr::map_dfr(names(file_paths), function(tissue) {
  path <- file_paths[[tissue]]
  if (file.exists(path)) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    if("X" %in% colnames(df)) {
      df <- rename(df, gene = X)
    }
    deg_filtered <- df %>%
      filter(!is.na(padj)) %>% 
      filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
      mutate(
        Tissue = tissue,
        diffexpressed = ifelse(log2FoldChange > 0, "UP", "DOWN")
      ) %>%
      dplyr::select(Tissue, gene, diffexpressed, log2FoldChange, padj, baseMean)
    return(deg_filtered)
  } else {
    warning(paste("Could not find file:", path))
    return(NULL)
  }
})

master_deg_table <- master_deg_table %>%
  arrange(Tissue, diffexpressed, padj)

output_file <- file.path(main_fig_dir, "Actual_DEGs_List_All_Tissues.csv")
write.csv(master_deg_table, output_file, row.names = FALSE)
message("Saved the complete list of actual genes to: ", output_file)

print_top_genes <- master_deg_table %>%
  group_by(Tissue, diffexpressed) %>%
  slice_head(n = 5) %>%
  dplyr::select(Tissue, diffexpressed, gene, log2FoldChange)

message("\nHere is a preview of the top actual genes per tissue:")
print(as.data.frame(print_top_genes))

top_3_down_per_tissue <- master_deg_table %>%
  filter(diffexpressed == "DOWN") %>%
  group_by(Tissue) %>%
  arrange(log2FoldChange, .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  dplyr::select(Tissue, gene, log2FoldChange, padj)

message("\n--- Top 3 Most Downregulated Genes Per Tissue ---")
print(as.data.frame(top_3_down_per_tissue))

highly_shared_genes <- master_deg_table %>%
  group_by(gene) %>%
  summarize(
    Num_Tissues = n_distinct(Tissue),
    Tissues = paste(Tissue, collapse = ", "),
    Directions = paste(diffexpressed, collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(Num_Tissues >= 5) %>%
  arrange(desc(Num_Tissues))

message("\n--- Genes Differentially Expressed in 5 or More Tissues ---")
print(as.data.frame(highly_shared_genes))

write.csv(
  highly_shared_genes, 
  file.path(main_fig_dir, "Highly_Shared_DEGs_5plus_Tissues.csv"), 
  row.names = FALSE
)