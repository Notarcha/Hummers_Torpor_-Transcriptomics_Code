library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(forcats)
library(here)
library(matrixStats)
library(reshape2)
library(ggplot2)
library(e1071)
library(patchwork)
library(relaimpo)
library(RColorBrewer)
library(lme4)
library(MuMIn)
library(broom.mixed)
library(lubridate)
library(pheatmap)

allowWGCNAThreads()

# Define project output directory for Lungs
out_dir <- here("Outputs", "WGCNA", "Lungs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

######## All tissue types all three states ########## 
########## Organismal level #############
### Fetch data and metadata
data <- read.csv(here("WGCNA_data", "data.csv"), row.names = 1)
metadata <- read.csv(here("WGCNA_data", "RNASeq_metadata.csv"), row.names = 1)

#### Remove outliers ####
outliers <- c("AS1", "AS53", "AS6")
metadata <- metadata[!(rownames(metadata) %in% outliers), ]
data <- data[, rownames(metadata)]

metadata <- metadata %>% filter(Tissue %in% c("Lungs"))
data <- data[, row.names(metadata)]
all(rownames(metadata) == colnames(data))
metadata_surfTemp_mlO2_min_euthTime <- read.csv(here("WGCNA_data", "metadata_surfTemp_mlO2_min_euthTime.csv"))

### detecting outliers
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detected as outliers
data <- data[gsg$goodGenes == TRUE, ]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# pca - method 2
pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var / sum(pca.var) * 100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

data <- data[, row.names(metadata)]
all(rownames(metadata) == colnames(data))

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ 1)

##### Keep genes where 75% of samples have a read count larger than 10
keep <- rowSums(counts(dds) >= 10) >= (0.75 * ncol(dds))
dds_filtered <- dds[keep, ]
nrow(dds_filtered)

### Normalize counts using variance stabilizing transformation
dds_rlog <- rlog(dds_filtered, fitType = "local")
norm.counts <- assay(dds_rlog) %>% 
  t()

write_csv(as.data.frame(norm.counts), file.path(out_dir, "Lungs_norm_counts_rlog.csv"))

############## Reducing noise in data prior to WGCNA ############################
gene_means <- colMeans(norm.counts)
gene_vars <- colVars(as.matrix(norm.counts))

jpeg(file.path(out_dir, "gene_means_histogram.jpg"), width = 1200, height = 900, res = 150) 
hist(gene_means, breaks = 100, main = "Distribution of Gene Means (rlog)", xlab = "Mean expression")
abline(v = 4, col = "blue", lty = 2)
abline(v = 5, col = "green", lty = 2)
abline(v = 5.5, col = "orange", lty = 2)
abline(v = 6, col = "red", lty = 2)
sapply(c(4, 5, 5.5, 6), function(cut) sum(gene_means > cut))
dev.off()

jpeg(file.path(out_dir, "gene_vars_histogram.jpg"), width = 1200, height = 900, res = 150) 
hist(gene_vars, breaks = 100)
dev.off()

# Step 2: Find the variance threshold for the top 80% most variable genes
var_threshold <- quantile(gene_vars, probs = 0.20)

# Step 3: Apply combined filter
keep_genes <- (gene_vars > var_threshold) & (gene_means > 6)

# Step 4: Subset filtered expression matrix
norm.counts.filtered <- norm.counts[, keep_genes]

gene_means_1 <- colMeans(norm.counts.filtered)
gene_vars_1 <- colVars(as.matrix(norm.counts.filtered))
hist(gene_means_1, breaks = 100)
hist(gene_vars_1, breaks = 100)

cat("Total genes before filtering:", ncol(norm.counts), "\n")
cat("Genes retained after filtering:", ncol(norm.counts.filtered), "\n")

# Plot to visualize retained genes
jpeg(file.path(out_dir, "genes_retained_postFiltering.jpg"), width = 1200, height = 900, res = 150) 
plot(gene_means, gene_vars,
     pch = 20,
     col = ifelse(keep_genes, "blue", "gray"),
     xlab = "Mean expression (VST)",
     ylab = "Variance",
     main = "Filtered Genes for WGCNA")
abline(h = var_threshold, col = "red", lty = 2)
abline(v = 6, col = "orange", lty = 2)
legend("topright", legend = c("Retained", "Filtered out"),
       col = c("blue", "gray"), pch = 20)
dev.off()

norm.counts.filtered <- as.data.frame(norm.counts.filtered)

## Network Construction ---------------------------------------------------
power <- c(1:25)

sft <- pickSoftThreshold(norm.counts.filtered,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

jpeg(file.path(out_dir, "softPowerThresholding.jpg"), width = 1200, height = 900, res = 150) 
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.75, color = 'black') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
dev.off()

norm.counts.filtered <- as.matrix(norm.counts.filtered)

cor <- WGCNA::cor
bicor <- WGCNA::bicor 

bwnet <- blockwiseModules(norm.counts.filtered,
                          maxBlockSize = 15000,
                          TOMType = "signed",
                          power = 17,
                          minModuleSize = 20,
                          mergeCutHeight = 0.3,
                          numericLabels = FALSE,
                          corType = "pearson",
                          randomSeed = 1234,
                          verbose = 3)

save(bwnet, file = file.path(out_dir, "Lungs_bwnet.RData"))

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs
moduleColors <- bwnet$colors

head(module_eigengenes)
table(bwnet$colors)
write.csv(table(bwnet$colors), file.path(out_dir, "Lungs_module_counts.csv"))

jpeg(file.path(out_dir, "ClusterDendrogram.jpg"), width = 1200, height = 900, res = 150) 
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)
dev.off()

# 6A. Relate modules to traits --------------------------------------------------
traits1 <- metadata %>% 
  mutate(Torpor = ifelse(grepl('D', Metabolic_State), 1, 0)) %>% 
  mutate(Transition = ifelse(grepl('T', Metabolic_State), 1, 0)) %>%
  mutate(Normothermy = ifelse(grepl('N', Metabolic_State), 1, 0))

traits1 <- traits1[, -c(1:4)]

nSamples <- nrow(norm.counts.filtered)
nGenes <- ncol(norm.counts.filtered)
module.trait.corr <- cor(module_eigengenes, traits1, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

melted_cor <- melt(module.trait.corr)
melted_pval <- melt(module.trait.corr.pvals)

plot_df <- cbind(melted_cor, p = melted_pval$value)
colnames(plot_df) <- c("Module", "Trait", "Correlation", "pvalue")

plot_df$label <- ifelse(plot_df$pvalue < 0.05,
                        paste0("p=", round(plot_df$pvalue, 2), "*"),
                        paste0("p=", round(plot_df$pvalue, 2)))

jpeg(file.path(out_dir, "modulesVsCategoricalTraits.jpg"), width = 1500, height = 1200, res = 150) 
ggplot(plot_df, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 4.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module–Metabolic Trait Correlation (p-values with significance)_Lungs", fill = "Correlation")
dev.off()

####### continuous traits of surfaceTemps and Metabolic rates
metadata <- metadata %>% rownames_to_column(var = "SampleID")

combined_metadata <- metadata %>%
  left_join(metadata_surfTemp_mlO2_min_euthTime, by = "BirdID")

traits1 <- combined_metadata %>%
  dplyr::filter(SampleID %in% rownames(module_eigengenes)) %>%
  dplyr::select(SampleID, mean_VO2, Surf_Temp) %>%
  drop_na() %>%
  column_to_rownames("SampleID")

if ("SampleID" %in% colnames(module_eigengenes)) {
  module_eigengenes <- module_eigengenes %>% select(-SampleID)
}

traits1 <- traits1[rownames(module_eigengenes), ]

nSamples <- nrow(traits1)
module.trait.corr <- cor(module_eigengenes, traits1, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

plot_df <- reshape2::melt(
  module.trait.corr,
  varnames = c("Module", "Trait"),
  value.name = "Correlation"
) %>%
  mutate(
    pvalue = reshape2::melt(
      module.trait.corr.pvals,
      varnames = c("Module", "Trait"),
      value.name = "pvalue"
    )$pvalue,
    label = ifelse(
      pvalue < 0.05,
      paste0("p=", round(pvalue, 2), "*"),
      paste0("p=", round(pvalue, 2))
    )
  )

jpeg(file.path(out_dir, "modulesVsContinuousTraits.jpg"), width = 1500, height = 1200, res = 150) 
ggplot(plot_df, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 4.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module–Trait Correlation with Metabolic Rate & Surface Temp",
       x = "Trait", y = "Module", fill = "Correlation")
dev.off()

### === EIGENGENE PLOT vs Continuous Trait === ###
if (!"SampleID" %in% colnames(metadata)) {
  metadata <- metadata %>% rownames_to_column(var = "SampleID")
}

traits1 <- metadata %>%
  left_join(metadata_surfTemp_mlO2_min_euthTime, by = "BirdID")

module_eigengenes$SampleID <- rownames(module_eigengenes)

plot_data <- left_join(
  module_eigengenes,
  traits1[, c("SampleID", "mean_VO2", "Surf_Temp", "Metabolic_State.y", "end_time", "BirdID")],
  by = "SampleID"
)

eigengene_cols <- grep("^ME", names(plot_data), value = TRUE)

long_df <- plot_data %>%
  pivot_longer(
    cols = all_of(eigengene_cols),
    names_to = "Module",
    values_to = "Eigengene"
  )

jpeg(file.path(out_dir, "modulesVsMeanVO2.jpg"), width = 1500, height = 1200, res = 150) 
ggplot(long_df, aes(x = mean_VO2, y = Eigengene, color = Metabolic_State.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  facet_wrap(~ Module, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Module Eigengene Expression vs Mean VO2",
       x = "Mean VO2 (ml/min)", y = "Module Eigengene Expression") +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

jpeg(file.path(out_dir, "modulesVsSurf_Temp.jpg"), width = 1500, height = 1200, res = 150)
ggplot(long_df, aes(x = Surf_Temp, y = Eigengene, color = Metabolic_State.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  facet_wrap(~ Module, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Module Eigengene Expression vs SurfTemp",
       x = "Surf Temp", y = "Module Eigengene Expression") +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

####################### GLMMs ###########################
eigengene_cols <- names(module_eigengenes) %>% str_subset("^ME[A-Za-z]+")

module_eigengenes <- module_eigengenes %>%
  { if (!"SampleID" %in% names(.)) tibble::rownames_to_column(., "SampleID") else . }

metadata <- metadata %>%
  { if (!"SampleID" %in% names(.)) tibble::rownames_to_column(., "SampleID") else . }

metadata_clean <- metadata %>% dplyr::select(-any_of(eigengene_cols))

plot_data <- module_eigengenes %>%
  left_join(metadata_clean, by = "SampleID") %>%
  left_join(metadata_surfTemp_mlO2_min_euthTime, by = "BirdID")

plot_data <- plot_data %>%
  mutate(
    end_time = ymd_hms(end_time, tz = "UTC"),
    h = lubridate::hour(end_time),
    m = lubridate::minute(end_time),
    s = lubridate::second(end_time),
    total_seconds = h * 3600 + m * 60 + s,
    shifted_seconds = (total_seconds - 22*3600) %% (24*3600),
    night_clock = sprintf("%02d:%02d:%02d",
                          shifted_seconds %/% 3600,
                          (shifted_seconds %% 3600) %/% 60,
                          shifted_seconds %% 60),
    night_hour = shifted_seconds / 3600,
    Surf_Temp = as.numeric(Surf_Temp),
    mean_VO2 = as.numeric(mean_VO2)
  )

lmm_data <- plot_data %>%
  dplyr::select(SampleID, mean_VO2, Surf_Temp, night_hour, all_of(eigengene_cols)) %>%
  drop_na(mean_VO2, Surf_Temp, night_hour)

fit_lmm <- function(module) {
  formula <- as.formula(paste0(module, " ~ mean_VO2 + Surf_Temp + night_hour"))
  model <- lm(formula, data = lmm_data)
  tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(Module = module)
}

lmm_results <- map_dfr(eigengene_cols, fit_lmm)
lmm_results$tissue <- "Lungs"
write.csv(lmm_results, file.path(out_dir, "Effects_predictors_moduleEigenGenes.csv"))

COLOR_VO2 <- "purple"
COLOR_SURF_TEMP <- "goldenrod2"
COLOR_NIGHT_HOUR <- "midnightblue"

order_and_plot <- function(data, predictor_term) {
  plot_data <- data %>%
    dplyr::filter(term == predictor_term) %>%
    dplyr::mutate(
      group_id = dplyr::case_when(
        conf.low > 0  ~ "3_Pos",
        conf.high < 0 ~ "2_Neg",
        TRUE          ~ "1_NS"
      )
    ) %>%
    dplyr::arrange(group_id, estimate) %>%
    dplyr::mutate(
      Module = factor(Module, levels = unique(Module))
    )
  
  p <- ggplot(plot_data,
              aes(x = estimate, y = Module,
                  xmin = conf.low, xmax = conf.high)) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1.0) +
    geom_pointrange(aes(color = term), fatten = 2, size = 0.8) + 
    labs(
      title = paste(predictor_term),
      x = "Effect Estimate (± 95% CI)",
      y = "Module"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_line(color = "grey90")
    )
  return(p)
}

p_vo2 <- order_and_plot(lmm_results, "mean_VO2") + scale_color_manual(values = c("mean_VO2" = COLOR_VO2))
p_surf_temp <- order_and_plot(lmm_results, "Surf_Temp") + scale_color_manual(values = c("Surf_Temp" = COLOR_SURF_TEMP))
p_night_hour <- order_and_plot(lmm_results, "night_hour") + scale_color_manual(values = c("night_hour" = COLOR_NIGHT_HOUR))

jpeg(file.path(out_dir, "Effects_predictors_moduleEigenGenes.jpg"), width = 1500, height = 1200, res = 150)
combined_plot <- (p_vo2 + p_surf_temp + p_night_hour) +
  plot_layout(nrow = 1, widths = c(1, 1, 1)) +
  plot_annotation(
    title = 'Effect of Predictors on Module Eigengenes',
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
print(combined_plot)
dev.off()

lm1 <- lm(mean_VO2 ~ Surf_Temp, lmm_data)
cor(lmm_data$mean_VO2, lmm_data$Surf_Temp, use = "complete.obs")
jpeg(file.path(out_dir, "VO2 vs Surf_Temp.jpg"), width = 1500, height = 1200, res = 150)
plot(lmm_data$Surf_Temp, lmm_data$mean_VO2, xlab = "Surf_Temp", ylab = "Mean VO2")
abline(lm1)
dev.off()

### Variance partitioning
var_part_results <- list()
for (mod in eigengene_cols) {
  fit <- lm(as.formula(paste(mod, "~ mean_VO2 + Surf_Temp + night_hour")), data = lmm_data)
  rel_imp <- calc.relimp(fit, type = "lmg", rela = TRUE)
  var_df <- tibble(
    Predictor = names(rel_imp$lmg),
    RelativeImportance = as.numeric(rel_imp$lmg),
    Module = mod
  )
  var_part_results[[mod]] <- var_df
}
var_part_df <- bind_rows(var_part_results)

jpeg(file.path(out_dir, "FixedEffects_variancePartitioning_EigenGenes.jpg"), width = 3000, height = 2000, res = 150)
ggplot(var_part_df, aes(x = Predictor, y = RelativeImportance, fill = Predictor)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Module, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Variance Partitioning of Module Eigengenes",
    x = "Predictor",
    y = "Proportion of Variance Explained"
  ) +
  theme(legend.position = "none")
dev.off()

run_module_pipeline <- function(
    module_name,
    predictor,
    base_dir = out_dir
) { 
  message("Running pipeline for module: ", module_name, " with predictor: ", predictor)
  genes <- sort(names(bwnet$colors[bwnet$colors == tolower(module_name)]))
  
  expr_long <- norm.counts.filtered %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SampleID") %>%
    dplyr::select(SampleID, all_of(genes)) %>%
    pivot_longer(-SampleID, names_to = "Gene", values_to = "Expression")
  
  predictor_source <- if (predictor %in% colnames(plot_data)) plot_data else combined_metadata
  
  expr_long <- expr_long %>%
    left_join(
      dplyr::select(predictor_source, SampleID, all_of(predictor)),
      by = "SampleID"
    )
  
  lm_results <- expr_long %>%
    group_by(Gene) %>%
    summarise(model = list(lm(Expression ~ .data[[predictor]])), .groups = "drop")
  
  tidy_res <- lm_results %>%
    mutate(tidy_res = map(model, broom::tidy)) %>%
    unnest(tidy_res) %>%
    filter(term != "(Intercept)") %>%
    transmute(Gene, estimate, std.error, statistic, p.value)
  
  glance_res <- lm_results %>%
    mutate(glance_res = map(model, broom::glance)) %>%
    unnest(glance_res) %>%
    transmute(Gene, r.squared, adj.r.squared)
  
  lm_results <- tidy_res %>%
    left_join(glance_res, by = "Gene") %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr"), Predictor = predictor)
  
  annotation_df <- lm_results %>%
    mutate(label = paste0("Slope=", signif(estimate, 3), ", R²=", signif(r.squared, 2), ", FDR=", signif(p.adj, 3)))
  
  plot_module_in_chunks <- function(data, predictor, chunk_size = 20) {
    genes <- unique(data$Gene)
    split_genes <- split(genes, ceiling(seq_along(genes) / chunk_size))
    plots <- map(split_genes, function(gene_subset) {
      ann_subset <- annotation_df %>% filter(Gene %in% gene_subset)
      ggplot(data %>% filter(Gene %in% gene_subset),
             aes_string(x = predictor, y = "Expression")) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        facet_wrap(~ Gene, scales = "free_y") +
        geom_text(
          data = ann_subset,
          aes(x = -Inf, y = Inf, label = label),
          hjust = -0.1, vjust = 1.2,
          inherit.aes = FALSE, size = 3.5, color = "blue"
        ) +
        theme_minimal(base_size = 14) +
        labs(title = paste(module_name, "Module:", predictor, "vs Expression"), x = predictor, y = "Expression")
    })
    return(plots)
  }
  
  plots <- plot_module_in_chunks(expr_long, predictor, chunk_size = 20)
  
  module_dir <- file.path(base_dir, paste0("Module", module_name))
  if (!dir.exists(module_dir)) dir.create(module_dir, recursive = TRUE)
  
  write.csv(lm_results, file = file.path(module_dir, paste0("lm_results_", predictor, ".csv")), row.names = FALSE)
  
  for (i in seq_along(plots)) {
    jpeg(file.path(module_dir, paste0(module_name, "_set", i, "_lm_", predictor, ".jpg")),
         width = 3000, height = 2400, res = 150)
    print(plots[[i]])
    dev.off()
  }
}

run_module_pipeline("Black", "mean_VO2")
run_module_pipeline("Pink", "mean_VO2")
run_module_pipeline("Lightyellow", "mean_VO2")
run_module_pipeline("Darkgreen", "mean_VO2")

run_module_pipeline("Darkgreen", "Surf_Temp")
run_module_pipeline("Lightyellow", "Surf_Temp")
run_module_pipeline("Salmon", "Surf_Temp")
run_module_pipeline("Black", "Surf_Temp")
run_module_pipeline("Pink", "Surf_Temp")

parent_dir <- out_dir
module_dirs <- list.dirs(path = parent_dir, full.names = TRUE, recursive = FALSE) %>%
  keep(~ str_detect(basename(.x), "^Module"))

compiled_df <- module_dirs %>%
  map_dfr(function(dir) {
    file_path <- file.path(dir, "lm_results_mean_VO2.csv")
    if (!file.exists(file_path)) return(NULL)
    read_csv(file_path, show_col_types = FALSE) %>% mutate(Module = basename(dir))
  })

sig_df <- compiled_df %>% filter(p.adj < 0.05) %>% arrange(p.adj)
write.csv(sig_df, file.path(out_dir, "genes_in_metabolic_pathways_GenesCorrelatedWith_Mean_VO2_FDR_0.05.csv"))

###### VO2 correlated #######
metadata <- data.frame(
  SampleID = rownames(norm.counts),
  Metabolic_State = c("N","T","D")[sample(1:3, nrow(norm.counts), replace=TRUE)]
)

norm_long <- norm.counts %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to="Gene", values_to="Expression") %>%
  left_join(metadata, by="SampleID")

metabolic_genes <- read_csv(file.path(out_dir, "genes_in_metabolic_pathways_GenesCorrelatedWith_Mean_VO2_FDR_0.05.csv"))
genes_of_interest <- metabolic_genes$Gene

norm_long <- norm_long %>%
  mutate(Metabolic_State = recode(Metabolic_State,
                                  "N" = "Normothermy",
                                  "T" = "Transition",
                                  "D" = "Deep Torpor"),
         Metabolic_State = factor(Metabolic_State, levels = c("Normothermy", "Transition", "Deep Torpor")))

anova_results <- norm_long %>%
  group_by(Gene) %>%
  summarise(
    p_value = tryCatch({
      m <- aov(Expression ~ Metabolic_State, data = cur_data())
      summary(m)[[1]][["Pr(>F)"]][1]
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_label = paste0("p = ", signif(p_value, 3)))

n_chunks <- ceiling(length(genes_of_interest) / 20)

for (i in seq_len(n_chunks)) {
  subset_genes <- genes_of_interest[((i - 1) * 20 + 1):min(i * 20, length(genes_of_interest))]
  
  df_plot <- norm_long %>%
    filter(Gene %in% subset_genes) %>%
    left_join(anova_results, by="Gene")
  
  p <- ggplot(df_plot, aes(x = Metabolic_State, y = Expression, fill = Metabolic_State)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    facet_wrap(~Gene + p_label, scales = "free_y") +
    scale_fill_manual(values = c("Normothermy" = "#1b9e77", "Transition" = "#d95f02", "Deep Torpor" = "#7570b3")) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(filename = file.path(out_dir, paste0("metabolic_genes_chunk_VO2_", i, ".pdf")), plot = p, width = 14, height = 10)
}

###### Surf_Temp correlated #######
compiled_df <- module_dirs %>%
  map_dfr(function(dir) {
    file_path <- file.path(dir, "lm_results_Surf_Temp.csv")
    if (!file.exists(file_path)) return(NULL)
    read_csv(file_path, show_col_types = FALSE) %>% mutate(Module = basename(dir))
  })

sig_df <- compiled_df %>% filter(p.adj < 0.05) %>% arrange(p.adj)
write.csv(sig_df, file.path(out_dir, "genes_in_metabolic_pathways_GenesCorrelatedWith_Surf_Temp_FDR_0.05.csv"))

metabolic_genes <- read_csv(file.path(out_dir, "genes_in_metabolic_pathways_GenesCorrelatedWith_Surf_Temp_FDR_0.05.csv"))
genes_of_interest <- metabolic_genes$Gene

n_chunks <- ceiling(length(genes_of_interest) / 20)

for (i in seq_len(n_chunks)) {
  subset_genes <- genes_of_interest[((i - 1) * 20 + 1):min(i * 20, length(genes_of_interest))]
  
  df_plot <- norm_long %>%
    filter(Gene %in% subset_genes) %>%
    left_join(anova_results, by="Gene")
  
  p <- ggplot(df_plot, aes(x = Metabolic_State, y = Expression, fill = Metabolic_State)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    facet_wrap(~Gene + p_label, scales = "free_y") +
    scale_fill_manual(values = c("Normothermy" = "goldenrod", "Transition" = "midnightblue", "Deep Torpor" = "purple")) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(filename = file.path(out_dir, paste0("metabolic_genes_chunk_Surf_Temp_", i, ".pdf")), plot = p, width = 14, height = 10)
}

##### Module membership ######
datExpr <- norm.counts.filtered
MEs_full <- module_eigengenes

MEs <- MEs_full[, !(colnames(MEs_full) %in% c("MEgrey", "SampleID")), drop = FALSE]

geneModuleMembership <- cor(datExpr, MEs, use = "p")
MMPvalue <- corPvalueStudent(geneModuleMembership, nrow(datExpr))

geneModuleMembership <- as.data.frame(geneModuleMembership)
MMPvalue <- as.data.frame(MMPvalue)

colnames(geneModuleMembership) <- paste0("kME.", colnames(MEs))
colnames(MMPvalue) <- paste0("p.kME.", colnames(MEs))

geneInfo <- data.frame(
  GeneID = colnames(datExpr),
  Assigned_Module = moduleColors,
  geneModuleMembership,
  MMPvalue,
  row.names = colnames(datExpr),
  check.names = FALSE
)

geneInfo$kME_Assigned <- NA
geneInfo$p_kME_Assigned <- NA

ME_colors_available <- tolower(gsub("^ME", "", colnames(MEs)))
geneInfo$Assigned_Module <- tolower(geneInfo$Assigned_Module)

for (module_color in ME_colors_available) {
  kME_col_name <- paste0("kME.ME", module_color)
  p_kME_col_name <- paste0("p.kME.ME", module_color)
  
  if (!(kME_col_name %in% names(geneInfo))) next
  
  is_gene_in_module <- geneInfo$Assigned_Module == module_color
  if (sum(is_gene_in_module) == 0) next
  
  geneInfo$kME_Assigned[is_gene_in_module] <- geneInfo[[kME_col_name]][is_gene_in_module]
  geneInfo$p_kME_Assigned[is_gene_in_module] <- geneInfo[[p_kME_col_name]][is_gene_in_module]
}

TARGET_MODULES <- c("black", "lightyellow", "pink", "darkgreen", "salmon")

final_report <- geneInfo[
  geneInfo$Assigned_Module %in% TARGET_MODULES,
  c("GeneID", "Assigned_Module", "kME_Assigned", "p_kME_Assigned")
]

final_report_sorted <- final_report[order(final_report$Assigned_Module, -final_report$kME_Assigned), ]
write.csv(final_report_sorted, file.path(out_dir, "IntraModuleMembership_Lungs.csv"))

######### Heatmap per target module #######
TARGET_MODULE <- "salmon"
genes_in_module <- names(bwnet$colors)[bwnet$colors == tolower(TARGET_MODULE)]

expr_module <- t(norm.counts.filtered[, genes_in_module, drop = FALSE])
expr_module_scaled <- t(scale(t(expr_module)))

metadata_plot <- metadata %>%
  dplyr::select(SampleID, Metabolic_State) %>%
  mutate(
    Metabolic_State = recode(Metabolic_State, "N" = "Normothermy", "T" = "Transition", "D" = "Deep Torpor"),
    Metabolic_State = factor(Metabolic_State, levels = c("Normothermy", "Transition", "Deep Torpor"))
  ) %>%
  column_to_rownames("SampleID")

common_samples <- intersect(colnames(expr_module_scaled), rownames(metadata_plot))
expr_module_scaled <- expr_module_scaled[, common_samples, drop = FALSE]
metadata_plot <- metadata_plot[common_samples, , drop = FALSE]

sample_order <- order(metadata_plot$Metabolic_State)
expr_module_scaled <- expr_module_scaled[, sample_order, drop = FALSE]
metadata_plot <- metadata_plot[sample_order, , drop = FALSE]

heatmap_colors <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)

annotation_colors <- list(
  Metabolic_State = c("Normothermy" = "#1b9e77", "Transition" = "#d95f02", "Deep Torpor" = "#7570b3")
)

pheatmap(
  expr_module_scaled,
  color = heatmap_colors,
  annotation_col = metadata_plot,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize = 12,
  main = paste("Heatmap of", TARGET_MODULE, "module genes"),
  filename = file.path(out_dir, paste0("Heatmap_", TARGET_MODULE, ".pdf")),
  width = 10,
  height = 12
)