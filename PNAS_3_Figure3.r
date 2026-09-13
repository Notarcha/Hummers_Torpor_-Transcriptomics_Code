############################################################
## FULL MULTI-TISSUE KEGG ENRICHMENT PIPELINE
## Species: Calypte anna
## Tissues: 7 Tissues (Gut, Heart, Lungs, Liver, Pectoral)
############################################################   

## 0. Load Libraries
suppressPackageStartupMessages({
  library(rtracklayer)    # GTF import
  library(dplyr)          # Data wrangling
  library(stringr)        # String handling
  library(readr)          # Reading tabular text
  library(clusterProfiler) # Enrichment analysis
  library(ggplot2)        # Plotting
  library(purrr)
  library(tibble)
  library(tidyr)
  library(here)
})

# Create main output directory
main_fig_dir <- here("Outputs")
dir.create(main_fig_dir, showWarnings = FALSE, recursive = TRUE)

## 1. Define Input Directories & Core Files
kegg_files_dir <- "Files for KEGG enrichments"

gtf_file  <- here(kegg_files_dir, "GCF_003957555.1_bCalAnn1_v1.p_genomic.gtf")
kaas_file <- here(kegg_files_dir, "Anna's Hummingbird_KAAS_run_birdrefs.txt")

cat("--> Loading Genome Annotation (GTF)...\n")
gtf <- import(gtf_file)

gene2protein <- gtf %>%
  as.data.frame() %>%
  filter(type == "CDS") %>%
  select(gene_id, protein_id) %>%
  distinct() %>%
  filter(!is.na(gene_id), !is.na(protein_id))

cat("Total genes in GTF:", n_distinct(gene2protein$gene_id), "\n")
cat("Total proteins in GTF:", n_distinct(gene2protein$protein_id), "\n")

## 2. Parse & Clean KAAS Output
cat("--> Processing KAAS Annotations...\n")
kaas_lines <- readLines(kaas_file)
kaas_lines <- kaas_lines[kaas_lines != "" & !grepl("^#", kaas_lines)]
kaas_split <- strsplit(kaas_lines, "\\s+")

kaas <- data.frame(
  raw_id = sapply(kaas_split, `[`, 1),
  KO     = sapply(kaas_split, function(x) if (length(x) >= 2) x[2] else NA),
  stringsAsFactors = FALSE
)

kaas$KO[kaas$KO == "-" | kaas$KO == ""] <- NA

kaas <- kaas %>%
  mutate(protein_id = str_extract(raw_id, "XP_[0-9]+\\.[0-9]+"))

# Map GTF Genes -> KAAS Proteins -> KO terms
gene2ko <- gene2protein %>%
  left_join(kaas %>% select(protein_id, KO), by = "protein_id") %>%
  filter(!is.na(KO)) %>%
  select(gene_id, KO) %>%
  distinct()

bg_ko_vec <- unique(gene2ko$KO)
cat("Background KO count:", length(bg_ko_vec), "\n")

## 3. Load KEGG Mappings (Pathways & Modules)
cat("--> Reading KEGG Reference Files...\n")

# Pathway Mapping
ko2pathway <- read.delim(here(kegg_files_dir, "ko_pathway.list"), header = FALSE, stringsAsFactors = FALSE)
colnames(ko2pathway) <- c("KO", "pathway")
ko2pathway$KO <- sub("ko:", "", ko2pathway$KO)
ko2pathway$pathway <- sub("path:map", "map", ko2pathway$pathway)

pathway_names <- read.delim(here(kegg_files_dir, "pathway_names.list"), header = FALSE, stringsAsFactors = FALSE)
colnames(pathway_names) <- c("pathway", "pathway_name")
pathway_names$pathway <- sub("path:map", "map", pathway_names$pathway)

ko2pathway <- ko2pathway %>% left_join(pathway_names, by = "pathway")

# Module Mapping
ko2module <- read.delim(here(kegg_files_dir, "ko_module_raw.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(ko2module) <- c("module", "KO")
ko2module$module <- sub("^md:", "", ko2module$module)
ko2module$KO     <- sub("^ko:", "", ko2module$KO)

module_names <- read.delim(here(kegg_files_dir, "ko_module_names.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(module_names) <- c("module", "module_name")
module_names$module <- sub("^md:", "", module_names$module)

ko2module <- ko2module %>% left_join(module_names, by = "module")

## 4. Define Tissue Order & Run Loop
tissue_order <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Heart", "Lungs", "Liver", "Pectoral")

cat("\n============================================================\n")
cat("STARTING ENRICHMENT LOOP FOR 7 TISSUES\n")
cat("============================================================\n")

for (tissue in tissue_order) {
  
  deg_file_path <- here(kegg_files_dir, paste0("Full_Results_", tissue, "_D_vs_N.csv"))
  
  if (!file.exists(deg_file_path)) {
    warning(paste("File not found for tissue:", tissue, "- Skipping."))
    next
  }
  
  cat("\n--> Processing:", tissue, "\n")
  
  deg_data <- read.csv(deg_file_path)
  
  # Extract DEGs (UP and DOWN)
  deg_genes <- deg_data %>%
    dplyr::filter(diffexpressed %in% c("UP", "DOWN")) %>%
    pull(gene)
  
  deg_genes <- deg_genes[deg_genes %in% gene2ko$gene_id]
  deg_ko_vec <- unique(gene2ko %>% filter(gene_id %in% deg_genes) %>% pull(KO))
  
  cat("    DEGs mapped to KO space:", length(deg_ko_vec), "\n")
  
  clean_tissue_name <- gsub(" ", "", tissue)
  
  if (length(deg_ko_vec) == 0) {
    warning(paste("No DEGs with KO terms found for:", tissue))
    next
  }
  
  # --- KEGG PATHWAY ENRICHMENT ---
  kk <- enricher(
    gene          = deg_ko_vec,
    universe      = bg_ko_vec,
    TERM2GENE     = ko2pathway %>% select(pathway, KO),
    TERM2NAME     = ko2pathway %>% select(pathway, pathway_name),
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  kk_results <- as.data.frame(kk)
  pathway_out <- file.path(main_fig_dir, paste0("KEGG_pathway_Enrichment_1_", clean_tissue_name, "_DvsN.csv"))
  write.csv(kk_results, pathway_out, row.names = FALSE)
  cat("    Saved Pathway Results:", nrow(kk_results), "terms found.\n")
  
  # Save Gene-to-Pathway mapping table if enriched terms exist
  if (nrow(kk_results) > 0) {
    genes_in_paths <- gene2ko %>%
      filter(gene_id %in% deg_genes) %>%
      inner_join(ko2pathway %>% filter(pathway %in% kk_results$ID), by = "KO", relationship = "many-to-many") %>%
      distinct()
    write.csv(genes_in_paths, file.path(main_fig_dir, paste0("KEGG_pathway_Enrichment_2_", clean_tissue_name, "_DvsN.csv")), row.names = FALSE)
  }
  
  # --- KEGG MODULE ENRICHMENT ---
  km <- enricher(
    gene          = deg_ko_vec,
    universe      = bg_ko_vec,
    TERM2GENE     = ko2module %>% select(module, KO),
    TERM2NAME     = ko2module %>% select(module, module_name),
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  km_results <- as.data.frame(km)
  module_out <- file.path(main_fig_dir, paste0("KEGG_module_Enrichment_1_", clean_tissue_name, "_DvsN.csv"))
  write.csv(km_results, module_out, row.names = FALSE)
  cat("    Saved Module Results:", nrow(km_results), "terms found.\n")
  
  # Save Gene-to-Module mapping table if enriched terms exist
  if (nrow(km_results) > 0) {
    genes_in_mods <- gene2ko %>%
      filter(gene_id %in% deg_genes) %>%
      inner_join(ko2module %>% filter(module %in% km_results$ID), by = "KO", relationship = "many-to-many") %>%
      distinct()
    write.csv(genes_in_mods, file.path(main_fig_dir, paste0("KEGG_module_Enrichment_2_", clean_tissue_name, "_DvsN.csv")), row.names = FALSE)
  }
}

## 5. Robust Plot Generator Function
cat("\n============================================================\n")
cat("GENERATING MULTI-TISSUE DOT PLOTS\n")
cat("============================================================\n")

process_and_plot_kegg <- function(search_keyword, plot_title, output_filename) {
  
  # Scan Outputs folder for files matching keyword (e.g. "pathway_Enrichment_1" or "module_Enrichment_1")
  all_files <- list.files(main_fig_dir, pattern = "\\.csv$", full.names = TRUE)
  files <- all_files[grepl(search_keyword, basename(all_files), ignore.case = TRUE) & grepl("_1_", basename(all_files))]
  
  if (length(files) == 0) {
    warning(paste("No matching CSV files found for keyword:", search_keyword))
    return(NULL)
  }
  
  all_data <- list()
  
  for (file in files) {
    fname <- basename(file)
    matched_tissue <- NA
    
    for (t in tissue_order) {
      clean_t <- gsub(" ", "", t)
      if (grepl(clean_t, fname, ignore.case = TRUE) || grepl(t, fname, ignore.case = TRUE)) {
        matched_tissue <- t
        break
      }
    }
    
    if (is.na(matched_tissue)) next
    
    df <- tryCatch({
      read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) return(data.frame()))
    
    # Must have rows and essential columns
    if (nrow(df) == 0 || !"Description" %in% colnames(df) || !"Count" %in% colnames(df)) next
    
    sort_col <- NULL
    if ("p.adjust" %in% colnames(df)) {
      sort_col <- "p.adjust"
    } else if ("pvalue" %in% colnames(df)) {
      sort_col <- "pvalue"
    }
    
    if (is.null(sort_col)) next
    
    df$p.adjust <- df[[sort_col]]
    df <- df %>% filter(!str_starts(Description, "path:"))
    
    if (nrow(df) > 0) {
      df$Tissue <- matched_tissue
      all_data[[fname]] <- df[, c("Description", "p.adjust", "Count", "Tissue")]
    }
  }
  
  if (length(all_data) == 0) {
    message(paste("No significant enriched terms found across any tissues for:", search_keyword, "- Skipping plot."))
    return(NULL)
  }
  
  big_df <- bind_rows(all_data)
  big_df$Tissue <- factor(big_df$Tissue, levels = tissue_order)
  
  # Select Top 5 terms per tissue ordered by p.adjust
  top_terms <- big_df %>%
    group_by(Tissue) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    pull(Description) %>%
    unique()
  
  plot_df <- big_df %>% filter(Description %in% top_terms)
  
  if (nrow(plot_df) == 0) return(NULL)
  
  p <- ggplot(plot_df, aes(x = Tissue, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    scale_x_discrete(drop = FALSE) +
    labs(title = plot_title, x = "", y = "", color = "FDR", size = "Genes") +
    theme_bw(base_size = 14) + 
    theme(
      plot.title       = element_text(color = "black", face = "bold", size = 18, hjust = 0.5, margin = margin(b=15)),
      axis.text.x      = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 14),
      axis.text.y      = element_text(color = "black", face = "bold", size = 14),
      legend.title     = element_text(color = "black", face = "bold", size = 14),
      legend.text      = element_text(color = "black", face = "bold", size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(main_fig_dir, output_filename), 
    plot = p, width = 13, height = 9, units = "in", dpi = 300, compression = "lzw"
  )
  
  message(paste("--> Successfully generated and saved:", output_filename))
}

# Generate Pathway Dot Plot
process_and_plot_kegg(
  search_keyword = "pathway",
  plot_title = "KEGG Pathway Enrichment Analysis",
  output_filename = "KEGG_Enrichment_Pathway_DotPlot.tiff"
)

# Generate Module Dot Plot
process_and_plot_kegg(
  search_keyword = "module",
  plot_title = "KEGG Module Enrichment Analysis",
  output_filename = "KEGG_Enrichment_Module_DotPlot.tiff"
)

cat("\n============================================================\n")
cat("PIPELINE COMPLETE! ALL OUTPUTS SAVED TO: ", main_fig_dir, "\n")
cat("============================================================\n")