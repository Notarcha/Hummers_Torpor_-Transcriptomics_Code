library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(here)

# ==============================================================================
# 0. SETUP OUTPUT DIRECTORY
# ==============================================================================

destination_dir <- "/Users/notarcha/Desktop/PNAS_code/Outputs/WGCNA/KEGG_Enrichments"

if (!dir.exists(destination_dir)) {
  dir.create(destination_dir, recursive = TRUE)
  message(paste("Created output directory:", destination_dir))
}

# Helper function to rename tissues to full display names
rename_tissue <- function(t_name) {
  case_when(
    t_name == "Gut1" ~ "Proximal Gut",
    t_name == "Gut2" ~ "Medial Gut",
    t_name == "Gut3" ~ "Distal Gut",
    t_name == "Pect" ~ "Pectoral",
    TRUE ~ t_name
  )
}

# ==============================================================================
# 1. BUILD GENE SYMBOL TO KO MAPPING VIA GTF
# ==============================================================================

annot_dir <- here("Files for KEGG enrichments")
gtf_file <- file.path(annot_dir, "GCF_003957555.1_bCalAnn1_v1.p_genomic.gtf")
kaas_file <- file.path(annot_dir, "Anna's Hummingbird_KAAS_run_birdrefs.txt")

message("Parsing GTF and KAAS annotations...")

# Parse GTF
gtf_lines <- readLines(gtf_file)
gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]

gtf_df <- data.frame(info = gtf_lines[grepl("protein_id", gtf_lines)], stringsAsFactors = FALSE) %>%
  mutate(
    Protein_ID = str_match(info, 'protein_id "([^"]+)"')[,2],
    Gene_Symbol = str_match(info, 'gene "([^"]+)"')[,2],
    Gene_Synonym = str_match(info, 'gene_synonym "([^"]+)"')[,2]
  )

gtf_map <- bind_rows(
  gtf_df %>% dplyr::select(Protein_ID, Gene = Gene_Symbol),
  gtf_df %>% dplyr::select(Protein_ID, Gene = Gene_Synonym)
) %>%
  filter(!is.na(Protein_ID) & !is.na(Gene) & Gene != "") %>%
  mutate(Gene = trimws(Gene), Protein_ID = trimws(Protein_ID)) %>%
  distinct()

# Parse KAAS
kaas_raw <- read.table(kaas_file, header = FALSE, sep = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = "")[, 1:2]
colnames(kaas_raw) <- c("Raw_Header", "KO")

kaas_gene_ko <- kaas_raw %>%
  filter(!is.na(KO) & KO != "") %>%
  mutate(
    Protein_ID = str_match(Raw_Header, "([XNA][PR]_\\d+\\.\\d+)")[,2],
    KO = trimws(gsub("^ko:", "", KO))
  ) %>%
  filter(!is.na(Protein_ID)) %>%
  inner_join(gtf_map, by = "Protein_ID", relationship = "many-to-many") %>%
  distinct(Gene, KO)

# ==============================================================================
# 2. LOAD KEGG PATHWAYS & MODULES
# ==============================================================================

# --- A. PATHWAYS ---
ko_to_path <- read_tsv(file.path(annot_dir, "ko_pathway.list"), col_names = c("KO", "Pathway_ID"), show_col_types = FALSE) %>% 
  mutate(
    KO = trimws(gsub("^ko:", "", KO)),
    Pathway_ID = trimws(gsub("^path:", "", Pathway_ID)),
    Path_Num = str_extract(Pathway_ID, "\\d{5}")
  ) %>%
  filter(!is.na(Path_Num))

path_names <- read_tsv(file.path(annot_dir, "pathway_names.list"), col_names = c("Pathway_ID", "Description"), show_col_types = FALSE) %>% 
  mutate(
    Pathway_ID = trimws(gsub("^path:", "", Pathway_ID)),
    Path_Num = str_extract(Pathway_ID, "\\d{5}")
  ) %>%
  filter(!is.na(Path_Num)) %>%
  distinct(Path_Num, Description)

gene_pathway_map <- kaas_gene_ko %>%
  inner_join(ko_to_path, by = "KO", relationship = "many-to-many") %>%
  inner_join(path_names, by = "Path_Num", relationship = "many-to-many") %>%
  dplyr::select(Gene, Pathway_ID = Path_Num, Description) %>%
  distinct()

# --- B. MODULES ---
ko_to_mod <- read_tsv(file.path(annot_dir, "ko_module_raw.txt"), col_names = c("Module_ID", "KO"), show_col_types = FALSE) %>% 
  mutate(
    KO = trimws(gsub("^ko:", "", KO)),
    Module_ID = trimws(gsub("^md:", "", Module_ID))
  )

mod_names <- read_tsv(file.path(annot_dir, "ko_module_names.txt"), col_names = c("Module_ID", "Description"), show_col_types = FALSE) %>% 
  mutate(Module_ID = trimws(gsub("^md:", "", Module_ID)))

gene_module_map <- kaas_gene_ko %>%
  inner_join(ko_to_mod, by = "KO", relationship = "many-to-many") %>%
  inner_join(mod_names, by = "Module_ID", relationship = "many-to-many") %>%
  dplyr::select(Gene, Pathway_ID = Module_ID, Description) %>%
  distinct()

N_pathway <- length(unique(gene_pathway_map$Gene))
N_module <- length(unique(gene_module_map$Gene))

message(paste("Master database built | Pathways mapped genes:", N_pathway, "| Modules mapped genes:", N_module))

# ==============================================================================
# 3. HYPERGEOMETRIC ENRICHMENT FUNCTION
# ==============================================================================

run_kegg_ora <- function(target_genes, annotation_map, N) {
  target_genes <- unique(target_genes[!is.na(target_genes)])
  overlap_map <- annotation_map %>% filter(Gene %in% target_genes)
  
  if (nrow(overlap_map) == 0) return(NULL)
  
  n <- length(unique(overlap_map$Gene))
  
  results <- overlap_map %>%
    group_by(Pathway_ID, Description) %>%
    summarise(
      Count = n_distinct(Gene),
      GeneIDs = paste(unique(Gene), collapse = "/"),
      .groups = "drop"
    ) %>%
    mutate(
      M = map_int(Pathway_ID, ~ n_distinct(annotation_map$Gene[annotation_map$Pathway_ID == .x])),
      pvalue = phyper(Count - 1, M, N - M, n, lower.tail = FALSE)
    ) %>%
    mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue)
  
  return(results)
}

# ==============================================================================
# 4. RECURSIVE DIRECTORY TRAVERSAL & WRITING TO KEGG_ENRICHMENTS
# ==============================================================================

wgcna_dir <- here("Outputs", "WGCNA")

all_lm_files <- list.files(wgcna_dir, pattern = "^lm_results.*\\.csv$", recursive = TRUE, full.names = TRUE)

message("\n--------------------------------------------------")
message(paste("Processing", length(all_lm_files), "WGCNA result files..."))
message("--------------------------------------------------")

saved_pathway_count <- 0
saved_module_count  <- 0

for (lm_file in all_lm_files) {
  
  lm_data <- read.csv(lm_file, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(lm_data) <- trimws(colnames(lm_data))
  
  gene_col <- colnames(lm_data)[grep("gene|id|symbol", colnames(lm_data), ignore.case = TRUE)][1]
  
  if (!is.na(gene_col)) {
    module_genes <- trimws(lm_data[[gene_col]])
    
    trait_name <- gsub("^lm_results_|_\\d+\\.csv$|\\.csv$", "", basename(lm_file))
    path_components <- strsplit(lm_file, "/")[[1]]
    m_name <- path_components[length(path_components) - 1]
    
    # Extract tissue name and map to full label
    t_raw <- path_components[length(path_components) - 2]
    t_name <- rename_tissue(t_raw)
    
    # 1. Pathway ORA
    res_pathway <- run_kegg_ora(module_genes, gene_pathway_map, N_pathway)
    if (!is.null(res_pathway) && nrow(res_pathway) > 0) {
      out_name <- paste0(t_name, "_", m_name, "_", trait_name, "_KEGGpathway_1.csv")
      write.csv(res_pathway, file.path(destination_dir, out_name), row.names = FALSE)
      message(paste("SUCCESS -> Saved Pathway:", out_name))
      saved_pathway_count <- saved_pathway_count + 1
    }
    
    # 2. Module ORA
    res_module <- run_kegg_ora(module_genes, gene_module_map, N_module)
    if (!is.null(res_module) && nrow(res_module) > 0) {
      out_name <- paste0(t_name, "_", m_name, "_", trait_name, "_KEGGmodule_1.csv")
      write.csv(res_module, file.path(destination_dir, out_name), row.names = FALSE)
      message(paste("SUCCESS -> Saved Module:", out_name))
      saved_module_count <- saved_module_count + 1
    }
  }
}

message("\n==================================================")
message(paste("FINISHED! Saved to:", destination_dir))
message(paste("Total Pathway CSVs:", saved_pathway_count, "| Total Module CSVs:", saved_module_count))
message("==================================================")