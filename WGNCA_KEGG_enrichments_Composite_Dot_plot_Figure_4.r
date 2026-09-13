library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(scales)

# ==============================================================================
# SET WORKING DIRECTORIES
# ==============================================================================

input_dir <- "/Users/notarcha/Desktop/PNAS_code/Outputs/WGCNA/KEGG_Enrichments"

if (!dir.exists(input_dir)) {
  stop(paste("Target folder does not exist:", input_dir))
}

all_tissues <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Heart", "Liver", "Lungs", "Pectoral")

# ==============================================================================
# FUNCTION TO GENERATE TISSUE-VS-TERM COMPOSITE DOT PLOT
# ==============================================================================

generate_composite_dotplot <- function(file_pattern, plot_title, output_prefix) {
  
  files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    warning(paste("No files matching pattern", file_pattern, "found in:", input_dir))
    return(NULL)
  }
  
  message(paste("Processing", length(files), "files for", plot_title, "..."))
  
  # Load and aggregate all files
  combined_df <- lapply(files, function(f) {
    file_name <- basename(f)
    parts <- unlist(strsplit(file_name, "_"))
    tissue_name <- parts[1]
    
    df <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      df$Tissue <- tissue_name
      return(df)
    }
    return(NULL)
  }) %>% bind_rows()
  
  if (nrow(combined_df) == 0) {
    warning(paste("No valid data rows found for", plot_title))
    return(NULL)
  }
  
  # Pool duplicate terms per tissue
  pooled_df <- combined_df %>%
    group_by(Tissue, Description) %>%
    summarise(
      p.adjust = min(p.adjust, na.rm = TRUE),
      pvalue = min(pvalue, na.rm = TRUE),
      Count = max(Count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(NegLogP = -log10(p.adjust))
  
  # Filter top 10 per tissue
  filtered_df <- pooled_df %>%
    filter(Count >= 2 & p.adjust < 0.05) %>%
    group_by(Tissue) %>%
    slice_max(order_by = NegLogP, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  if (nrow(filtered_df) < 3) {
    message("Including top terms by raw p-value...")
    filtered_df <- pooled_df %>%
      filter(Count >= 2) %>%
      mutate(p.adjust = pvalue, NegLogP = -log10(pvalue)) %>%
      group_by(Tissue) %>%
      slice_min(order_by = pvalue, n = 10, with_ties = FALSE) %>%
      ungroup()
  }
  
  filtered_df <- filtered_df %>%
    mutate(Tissue = factor(Tissue, levels = all_tissues))
  
  # Generate Dot Plot with regular log10 scale ticks
  p <- ggplot(filtered_df, aes(x = Tissue, y = reorder(Description, NegLogP))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_x_discrete(drop = FALSE) +
    
    # --- CLEAN & REGULAR P-ADJUST SCALE ---
    scale_color_continuous(
      low = "firebrick1", 
      high = "royalblue", 
      trans = "log10",
      name = "p.adjust",
      breaks = 10^(-seq(2, 10, by = 2)),
      labels = label_log(digits = 2)
    ) +
    
    scale_size_continuous(range = c(4, 11), name = "Gene Count (≥2)") +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 14) +
    labs(
      title = str_wrap(plot_title, width = 45),
      subtitle = str_wrap("Pooled across modules & traits per tissue (Filtered for Count ≥ 2 & p.adj < 0.05)", width = 60),
      x = "Tissue Type",
      y = "Enriched Term"
    ) +
    theme(
      # Title & Subtitle settings with high padding
      plot.title = element_text(face = "bold", size = 18, margin = margin(b = 8), hjust = 0),
      plot.subtitle = element_text(size = 13, margin = margin(b = 15), hjust = 0),
      
      # Axis Labels
      axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 12)),
      axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
      
      # Axis Tick Text
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 13, color = "black"),
      axis.text.y = element_text(face = "bold", size = 11, color = "black"),
      
      # Legend Formatting
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      
      # Grid and Large Outer Plot Margins
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      plot.margin = margin(t = 40, r = 25, b = 20, l = 20)
    )
  
  out_pdf <- file.path(input_dir, paste0(output_prefix, "_Composite_DotPlot_Min2Genes.pdf"))
  out_png <- file.path(input_dir, paste0(output_prefix, "_Composite_DotPlot_Min2Genes.png"))
  
  ggsave(out_pdf, plot = p, width = 12.5, height = 13, units = "in")
  ggsave(out_png, plot = p, width = 12.5, height = 13, units = "in", dpi = 300)
  
  message(paste("Saved plot:", basename(out_png)))
}

# ==============================================================================
# EXECUTE FOR PATHWAYS AND MODULES
# ==============================================================================

generate_composite_dotplot(
  file_pattern = ".*_KEGGpathway_1\\.csv$",
  plot_title = "KEGG Pathways — Cross-Tissue Composite",
  output_prefix = "AllTissues_KEGGpathway"
)

generate_composite_dotplot(
  file_pattern = ".*_KEGGmodule_1\\.csv$",
  plot_title = "KEGG Modules — Cross-Tissue Composite",
  output_prefix = "AllTissues_KEGGmodule"
)

message("\n==================================================")
message(paste("Finished! Dot plots with formatted p.adjust legend saved to:\n", input_dir))
message("==================================================")