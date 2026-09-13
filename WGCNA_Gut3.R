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

allowWGCNAThreads()

# Output directories
output_root <- here("Outputs/WGCNA")
output_dir <- here("Outputs/WGCNA/Gut3")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

######## All tissue types all three states ########## 
########## Organismal level #############
### Fetch data and metadata
data<-read.csv(here("WGCNA_data/data.csv"), row.names=1)
metadata<-read.csv(here("WGCNA_data/RNASeq_metadata.csv"),row.names=1)

#### Remove outliers ####
outliers <- c("AS1", "AS53", "AS6")
metadata<- metadata[!(rownames(metadata) %in% outliers), ]
data<- data[, rownames(metadata)]

#metadata<-metadata %>% filter(!Metabolic_State =="T")
metadata<-metadata %>% filter(Tissue %in% c("Gut3"))
data<-data[,row.names(metadata)]
all(rownames(metadata) == colnames(data))
metadata_surfTemp_mlO2_min_euthTime <- read.csv(here("WGCNA_data/metadata_surfTemp_mlO2_min_euthTime.csv"))

### detecting outliers
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


#### Outliers not very apparent. We are good! 

# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead

data<-data[,row.names(metadata)]
all(rownames(metadata) == colnames(data))

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ 1) # not specifying model

##### Keep genes where 75% of samples have a read count larger than 10
keep <- rowSums(counts(dds) >= 10) >= (0.75 * ncol(dds))
dds_filtered <- dds[keep, ]
nrow(dds_filtered) # 12260 genes

### Normalize counts using variance stabilizing transformation
dds_rlog <- rlog(dds_filtered, fitType = "local")
norm.counts <- assay(dds_rlog) %>% 
  t()

write_csv(as.data.frame(norm.counts),here("Outputs/WGCNA/Gut3/Gut3_norm_counts_rlog.csv"))

############## Reducing noise in data prior to WGCNA ############################
##### Idea : Use only genes with moderate expression values and ones that are highly variable for WGCNA ######## - Hans Suggestion

# Load required packages
library(matrixStats)
library(reshape2)
library(ggplot2)
library(e1071) # Load required library
library(matrixStats)

# Input: norm.counts (normalized), transposed: rows = samples, cols = genes)

# Step 1: Calculate gene-wise statistics
gene_means <- colMeans(norm.counts)
gene_vars <- colVars(as.matrix(norm.counts))

jpeg(here("Outputs/WGCNA/Gut3/gene_means_histogram.jpg"), width = 1200, height = 900, res = 150) 
hist(gene_means, breaks = 100, main = "Distribution of Gene Means (rlog)", xlab = "Mean expression")
abline(v = 4, col = "blue", lty = 2)
abline(v = 5, col = "green", lty = 2)
abline(v = 5.5, col = "orange", lty = 2)
abline(v = 6, col = "red", lty = 2)
sapply(c(4, 5, 5.5, 6), function(cut) sum(gene_means > cut))
dev.off ()

jpeg(here("Outputs/WGCNA/Gut3/gene_vars_histogram.jpg"), width = 1200, height = 900, res = 150) 
hist(gene_vars,breaks=100)
dev.off ()

# Step 2: Find the variance threshold for the top 80% most variable genes
var_threshold <- quantile(gene_vars, probs = 0.20)

# Step 3: Apply combined filter
keep_genes <- (gene_vars > var_threshold) & (gene_means > 6)

# Step 4: Subset filtered expression matrix
norm.counts.filtered <- norm.counts[, keep_genes]

gene_means_1 <- colMeans(norm.counts.filtered)
gene_vars_1 <- colVars(as.matrix(norm.counts.filtered))
hist(gene_means_1,breaks=100)
hist(gene_vars_1,breaks=100)

# Step 5: Print summary
cat("Total genes before filtering:", ncol(norm.counts), "\n")
cat("Genes retained after filtering:", ncol(norm.counts.filtered), "\n")

######### 7860 genes #########

# Plot to visualize retained genes
jpeg(here("Outputs/WGCNA/Gut3/genes_retained_postFiltering.jpg"), width = 1200, height = 900, res = 150) 
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

# # ---------- 2. Sample genes and plot distributions ----------
# set.seed(123)
# sample_genes <- sample(colnames(norm.counts.filtered), 16)
# 
# # Histograms of expression across samples for selected genes
# par(mfrow = c(4, 4))  # 4x4 grid
# for (gene in sample_genes) {
#   hist(norm.counts.filtered[, gene],
#        main = gene, xlab = "Expression", col = "lightblue", breaks = 10)
# }
# 
# # QQ-plots for the same genes
# par(mfrow = c(4, 4))
# for (gene in sample_genes) {
#   qqnorm(norm.counts.filtered[, gene], main = paste("QQ:", gene))
#   qqline(norm.counts.filtered[, gene], col = "red")
# }
# 
# # ---------- 3. Plot overall density ----------
# df_long <- melt(norm.counts.filtered)
# 
# ggplot(df_long, aes(x = value)) +
#   geom_density(fill = "steelblue", alpha = 0.5) +
#   labs(title = "Density of Filtered Expression Values",
#        x = "Normalized Expression", y = "Density") +
#   theme_minimal()
# 
# # ---------- 4. (Optional) Skewness check ----------
# gene_skews <- apply(norm.counts.filtered, 2, skewness)
# 
# # Histogram of skewness values
# hist(gene_skews, breaks = 30,
#      main = "Skewness Across Genes", xlab = "Skewness")
# 
# # Quick summary
# summary(gene_skews)

#####
norm.counts.filtered<-as.data.frame(norm.counts.filtered)

## Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(1:25)

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts.filtered,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power
jpeg(here("Outputs/WGCNA/Gut3/softPowerThresholding.jpg"), width = 1200, height = 900, res = 150) 
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

#### choosing power of 7 and continuing ####

# convert matrix to numeric
norm.counts.filtered<-as.matrix(norm.counts.filtered)

# IMPORTANT: Resolve WGCNA namespace conflict for 'cor' and 'bicor'
# This ensures WGCNA uses its own internal correlation functions
cor <- WGCNA::cor
bicor <- WGCNA::bicor # Add this line as well, just in case bicor is used internally by WGCNA for certain checks or calculations even with corType="pearson"


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts.filtered,
                          maxBlockSize = 15000,
                          TOMType = "signed",
                          power = 7,
                          minModuleSize = 20,
                          mergeCutHeight = 0.3,
                          numericLabels = FALSE,
                          corType = "pearson", # Correctly specified corType
                          randomSeed = 1234,
                          verbose = 3)

save(bwnet, file = here("Outputs/WGCNA/Gut3/Gut3_bwnet.RData"))

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs
moduleColors <- bwnet$colors

# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)
write.csv(table(bwnet$colors),here("Outputs/WGCNA/Gut3.csv"))


# Plot the dendrogram and the module colors before and after merging underneath
jpeg(here("Outputs/WGCNA/Gut3/ClusterDendrogram.jpg"), width = 1200, height = 900, res = 150) 
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
dev.off()
# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

### Categorical traits
# create traits file - binarize categorical variables
### Metabolic states
traits1 <- metadata%>% 
  mutate(Torpor= ifelse(grepl('D', Metabolic_State), 1, 0)) %>% 
  mutate(Transition= ifelse(grepl('T', Metabolic_State), 1, 0)) %>%
  mutate(Normothermy= ifelse(grepl('N', Metabolic_State), 1, 0))

traits1<-traits1[,-c(1:4)]

# Define numbers of genes and samples
nSamples <- nrow(norm.counts.filtered)
nGenes <- ncol(norm.counts.filtered)
module.trait.corr <- cor(module_eigengenes, traits1, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

library(reshape2)
library(ggplot2)

# Melt correlation and p-value matrices
melted_cor <- melt(module.trait.corr)
melted_pval <- melt(module.trait.corr.pvals)

# Combine into a single dataframe
plot_df <- cbind(melted_cor, p = melted_pval$value)
colnames(plot_df) <- c("Module", "Trait", "Correlation", "pvalue")

# Create label: round p-value + asterisk if < 0.05
plot_df$label <- ifelse(plot_df$pvalue < 0.05,
                        paste0("p=", round(plot_df$pvalue, 2), "*"),
                        paste0("p=", round(plot_df$pvalue, 2)))

# Plot heatmap
jpeg(here("Outputs/WGCNA/Gut3/modulesVsCategoricalTraits.jpg"), width = 1500, height = 1200, res = 150) 
ggplot(plot_df, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 4.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module–Metabolic Trait Correlation (p-values with significance)_Gut3", fill = "Correlation")
dev.off()

####### continuous traits of surfaceTemps and Metabolic rates
library(WGCNA)
library(tidyverse)
library(reshape2)
library(ggplot2)

# Step 1: Move sample IDs from rownames into a column in metadata
metadata <- metadata %>%
  rownames_to_column(var = "SampleID")  # Converts rownames (AS1, AS2...) into a column

# Step 2: Merge tissue-level metadata with bird-level traits
combined_metadata <- metadata %>%
  left_join(metadata_surfTemp_mlO2_min_euthTime, by = "BirdID")

# Step 3: Prepare trait matrix using SampleID
traits1 <- combined_metadata %>%
  dplyr::filter(SampleID %in% rownames(module_eigengenes)) %>%
  dplyr::select(SampleID, mean_VO2, Surf_Temp) %>%
  drop_na() %>%
  column_to_rownames("SampleID")

# Step 4: Align order with module eigengenes
if ("SampleID" %in% colnames(module_eigengenes)) {
  module_eigengenes <- module_eigengenes %>% select(-SampleID)
}

traits1 <- traits1[rownames(module_eigengenes), ]

# Step 5: Correlation and p-values
nSamples <- nrow(traits1)
module.trait.corr <- cor(module_eigengenes, traits1, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Step 6: Format for plotting
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


# Step 7: Plot heatmap
jpeg(here("Outputs/WGCNA/Gut3//modulesVsContinuousTraits.jpg"), width = 1500, height = 1200, res = 150) 
ggplot(plot_df, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 4.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module–Trait Correlation with Metabolic Rate & Surface Temp",
       x = "Trait", y = "Module", fill = "Correlation")  # y-axis title bold

dev.off()

### so cool! Lots of modules important! 
### Much power in using quantitative trait data

### === EIGENGENE PLOT vs Continuous Trait === ###

library(tidyverse)

# 1. Ensure SampleID is present in metadata
if (!"SampleID" %in% colnames(metadata)) {
  metadata <- metadata %>% rownames_to_column(var = "SampleID")
}

# 2. Merge metadata with physiological traits
traits1 <- metadata %>%
  left_join(metadata_surfTemp_mlO2_min_euthTime, by = "BirdID")

# 3. Add SampleID to module eigengenes
module_eigengenes$SampleID <- rownames(module_eigengenes)

# 4. Merge eigengenes with physiological traits
plot_data <- left_join(
  module_eigengenes,
  traits1[, c("SampleID", "mean_VO2", "Surf_Temp", "Metabolic_State.y","end_time","BirdID")],
  by = "SampleID"
)

# 5. Extract eigengene columns
eigengene_cols <- grep("^ME", names(plot_data), value = TRUE)

# 6. Reshape for plotting
long_df <- plot_data %>%
  pivot_longer(
    cols = all_of(eigengene_cols),
    names_to = "Module",
    values_to = "Eigengene"
  )

# 7. Plot eigengene vs mean VO₂
jpeg(here("Outputs/WGCNA/Gut3/modulesVsMeanVO2.jpg"), width = 1500, height = 1200, res = 150) 
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

# 8. Plot eigengene vs surfTemp
jpeg(here("Outputs/WGCNA/Gut3/modulesVsSurf_Temp.jpg"), width = 1500, height = 1200, res = 150)
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
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(broom.mixed)
library(ggplot2)
library(purrr)
library(lubridate)

### Step 1: Identify eigengene columns
eigengene_cols <- names(module_eigengenes) %>% str_subset("^ME[A-Za-z]+")

### Step 2: Ensure SampleID exists
module_eigengenes <- module_eigengenes %>%
  { if (!"SampleID" %in% names(.)) tibble::rownames_to_column(., "SampleID") else . }

metadata <- metadata %>%
  { if (!"SampleID" %in% names(.)) tibble::rownames_to_column(., "SampleID") else . }

### Step 3: Merge metadata and physiology
metadata_clean <- metadata %>% dplyr::select(-any_of(eigengene_cols))

plot_data <- module_eigengenes %>%
  left_join(metadata_clean, by = "SampleID") %>%
  left_join(metadata_surfTemp_mlO2_min_euthTime, by = "BirdID")

### Step 4: Define linear shifted clock (22:00 = 00:00)
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

### Step 5: Prepare modeling dataset
lmm_data <- plot_data %>%
  dplyr::select(SampleID, mean_VO2, Surf_Temp, night_hour, all_of(eigengene_cols)) %>%
  drop_na(mean_VO2, Surf_Temp, night_hour)

### Step 6: Fit LM per module (since 1 row per bird)
fit_lmm <- function(module) {
  formula <- as.formula(paste0(module, " ~ mean_VO2 + Surf_Temp + night_hour"))
  model <- lm(formula, data = lmm_data)
  tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(Module = module)
}

lmm_results <- map_dfr(eigengene_cols, fit_lmm)
lmm_results$tissue<-"Gut3"
write.csv(lmm_results,here("Outputs/WGCNA/Gut3/Effects_predictors_moduleEigenGenes.csv"))

library(tidyverse)
library(patchwork) 
library(here)

# Define the colors for clarity
COLOR_VO2 <- "purple"
COLOR_SURF_TEMP <- "goldenrod2"
COLOR_NIGHT_HOUR <- "midnightblue"

# --- 1. Ordering and Plotting Function (Fixed) ---
order_and_plot <- function(data, predictor_term) {
  
  # 1. Prepare data with explicit grouping
  plot_data <- data %>%
    dplyr::filter(term == predictor_term) %>%
    dplyr::mutate(
      # Define Groups: 3=Top(Pos), 2=Middle(Neg), 1=Bottom(NS)
      group_id = dplyr::case_when(
        conf.low > 0  ~ "3_Pos",  # Significant Positive
        conf.high < 0 ~ "2_Neg",  # Significant Negative
        TRUE          ~ "1_NS"    # Not Significant
      )
    ) %>%
    # 2. Sort the Dataframe (Bottom of plot -> Top of plot)
    # Sorting by 'group_id' puts NS first (bottom), then Neg, then Pos (top).
    # Sorting by 'estimate' ensures values flow naturally (-0.5 is below -0.1).
    dplyr::arrange(group_id, estimate) %>%
    dplyr::mutate(
      # Lock the order by converting Module to a factor based on the sorted rows
      Module = factor(Module, levels = unique(Module))
    )
  
  # 3. Create the ggplot object
  p <- ggplot(plot_data,
              aes(x = estimate, y = Module,
                  xmin = conf.low, xmax = conf.high)) +
    
    # Vertical zero line
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1.0) +
    
    # Points and Error bars
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
      panel.grid.major.y = element_line(color = "grey90") # Add light lines for readability
    )
  
  return(p)
}

# --- 2. Generate and Apply Specific Colors ---

p_vo2 <- order_and_plot(lmm_results, "mean_VO2") +
  scale_color_manual(values = c("mean_VO2" = COLOR_VO2))

p_surf_temp <- order_and_plot(lmm_results, "Surf_Temp") +
  scale_color_manual(values = c("Surf_Temp" = COLOR_SURF_TEMP))

p_night_hour <- order_and_plot(lmm_results, "night_hour") +
  scale_color_manual(values = c("night_hour" = COLOR_NIGHT_HOUR))

# --- 3. Combine Plots with Patchwork ---
jpeg(here("Outputs/WGCNA/Gut3/Effects_predictors_moduleEigenGenes.jpg"), width = 1500, height = 1200, res = 150)

combined_plot <- (p_vo2 + p_surf_temp + p_night_hour) +
  plot_layout(nrow = 1, widths = c(1, 1, 1)) +
  plot_annotation(
    title = 'Effect of Predictors on Module Eigengenes',
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

print(combined_plot)
dev.off()


lm1<-lm(mean_VO2~Surf_Temp,lmm_data)
cor(lmm_data$mean_VO2, lmm_data$Surf_Temp, use = "complete.obs")
jpeg(here("Outputs/WGCNA/Gut3/VO2 vs Surf_Temp.jpg"), width = 1500, height = 1200, res = 150)
plot(lmm_data$Surf_Temp,lmm_data$mean_VO2, xlab="Surf_Temp", ylab="Mean VO2")
abline(lm1)
dev.off()

### Step 6: Plot variance partitioning
library(relaimpo)
library(tidyverse)
library(RColorBrewer)
library(lme4)
library(MuMIn)
library(here)

### Step 1 — Predictor-level variance partitioning (LMG)
var_part_results <- list()

for (mod in eigengene_cols) {
  fit <- lm(as.formula(paste(mod, "~ mean_VO2 + Surf_Temp + night_hour")), 
            data = lmm_data)
  
  rel_imp <- calc.relimp(fit, type = "lmg", rela = TRUE)
  
  var_df <- tibble(
    Predictor = names(rel_imp$lmg),
    RelativeImportance = as.numeric(rel_imp$lmg),
    Module = mod
  )
  
  var_part_results[[mod]] <- var_df
}

var_part_df <- bind_rows(var_part_results)

# Plot fixed effects variance partitioning
jpeg(here("Outputs/WGCNA/Gut3/FixedEffects_variancePartitioning_EigenGenes.jpg"), 
     width = 3000, height = 2000, res = 150)
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

###### Lesssss Go! Let's inspect Gut3 modules we want carefully. 
### VO2 --- Turquoise, Tan, Salmon, Purple, MidnightBlue, Brown  ###
### Night_Hour - Purple #### 
### Surf_Temp - Tan, Salmon, Purple, MidnightBlue, Green ######

run_module_pipeline <- function(
    module_name,
    predictor,
    base_dir = here("Outputs/WGCNA", "Gut3")
) { 
  message("Running pipeline for module: ", module_name, " with predictor: ", predictor)
  
  # ---- Genes from module ----
  genes <- sort(names(bwnet$colors[bwnet$colors == tolower(module_name)]))
  
  # ---- Expression data ----
  expr_long <- norm.counts.filtered %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SampleID") %>%
    dplyr::select(SampleID, all_of(genes)) %>%
    pivot_longer(-SampleID, names_to = "Gene", values_to = "Expression")
  
  # ---- Pick correct metadata source ----
  predictor_source <- if (predictor %in% colnames(plot_data)) plot_data else combined_metadata
  
  expr_long <- expr_long %>%
    left_join(
      dplyr::select(predictor_source, SampleID, all_of(predictor)),
      by = "SampleID"
    )
  
  # ---- Linear regression per gene ----
  lm_results <- expr_long %>%
    group_by(Gene) %>%
    summarise(model = list(lm(Expression ~ .data[[predictor]])), .groups = "drop")
  
  # Slopes and p-values
  tidy_res <- lm_results %>%
    mutate(tidy_res = map(model, broom::tidy)) %>%
    unnest(tidy_res) %>%
    filter(term != "(Intercept)") %>%
    transmute(
      Gene,
      estimate,
      std.error,
      statistic,
      p.value
    )
  
  # Model fit (R² etc.)
  glance_res <- lm_results %>%
    mutate(glance_res = map(model, broom::glance)) %>%
    unnest(glance_res) %>%
    transmute(
      Gene,
      r.squared,
      adj.r.squared
    )
  
  # Join everything
  lm_results <- tidy_res %>%
    left_join(glance_res, by = "Gene") %>%
    mutate(
      p.adj = p.adjust(p.value, method = "fdr"),
      Predictor = predictor
    )
  
  
  # ---- Annotations for plots ----
  annotation_df <- lm_results %>%
    mutate(label = paste0(
      "Slope=", signif(estimate, 3),
      ", R²=", signif(r.squared, 2),
      ", FDR=", signif(p.adj, 3)
    ))
  
  # ---- Plot function ----
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
        labs(
          title = paste(module_name, "Module:", predictor, "vs Expression"),
          x = predictor, y = "Expression"
        )
    })
    return(plots)
  }
  
  plots <- plot_module_in_chunks(expr_long, predictor, chunk_size = 20)
  
  # ---- Output directories ----
  module_dir <- file.path(base_dir, paste0("Module", module_name))
  if (!dir.exists(module_dir)) dir.create(module_dir, recursive = TRUE)
  
  # ---- Save results ----
  write.csv(lm_results,
            file = file.path(module_dir, paste0("lm_results_", predictor, ".csv")),
            row.names = FALSE)
  
  for (i in seq_along(plots)) {
    jpeg(file.path(module_dir, paste0(module_name, "_set", i, "_lm_", predictor, ".jpg")),
         width = 3000, height = 2400, res = 150)
    print(plots[[i]])
    dev.off()
  }
  
  message("Saved results to: ", module_dir)
}

### Surf_Temp - Tan, Salmon, Purple, MidnightBlue, Green ######
run_module_pipeline("Magenta", "Surf_Temp")

# ########### Looking at metabolic genes ######### 
# 
# # Set the parent directory containing Module* folders
# parent_dir <- here("Outputs/WGCNA", "Gut3")
# 
# # Identify Module folders
# module_dirs <- list.dirs(
#   path = parent_dir,
#   full.names = TRUE,
#   recursive = FALSE
# ) %>%
#   keep(~ str_detect(basename(.x), "^Module"))
# 
# # Read and compile lm_results_mean_VO2.csv across modules
# compiled_df <- module_dirs %>%
#   map_dfr(function(dir) {
#     
#     file_path <- file.path(dir, "lm_results_mean_VO2.csv")
#     
#     if (!file.exists(file_path)) {
#       return(NULL)
#     }
#     
#     read_csv(file_path, show_col_types = FALSE) %>%
#       mutate(
#         Module = basename(dir)
#       )
#   })
# 
# # Sanity check
# stopifnot(nrow(compiled_df) > 0)
# 
# # Filter significant results
# sig_df <- compiled_df %>%
#   filter(p.adj < 0.05)
# 
# # Optional: arrange for readability
# sig_df <- sig_df %>%
#   arrange(p.adj)
# 
# sig_df
# 
# write.csv(sig_df,here("Outputs/WGCNA/Gut3/genes_in_metabolic_pathways_GenesCorrelatedWith_Mean_VO2_FDR_0.05.csv"))

########## No genes have expression correlated with VO2 ########### 
##### Let go of this ######### 

# ###### VO2 correlated #######
# 
# library(tidyverse)
# 
# # Suppose norm.counts is your data frame:
# # rows = samples, cols = genes
# # rownames(norm.counts) = sample IDs like AS5, AS6, etc.
# 
# library(tidyverse)
# 
# # --------------------------
# # 1. Example: your normalized counts
# # rows = samples, cols = genes
# #rownames(norm.counts) = sample IDs like "AS5", "AS6", etc.
# 
# # --------------------------
# # 2. Your sample metadata (replace with your actual mapping!)
# metadata <- data.frame(
#   SampleID = rownames(norm.counts),
#   Metabolic_State = c("N","T","D")[sample(1:3, nrow(norm.counts), replace=TRUE)]
# )
# 
# # --------------------------
# # 1. Reshape counts + join metadata
# norm_long <- norm.counts %>%
#   as.data.frame() %>% 
#   rownames_to_column("SampleID") %>%
#   pivot_longer(-SampleID, names_to="Gene", values_to="Expression") %>%
#   left_join(metadata, by="SampleID")
# 
# # --------------------------
# # 2. Genes of interest
# metabolic_genes <- read_csv(here("Outputs/WGCNA/Gut3/genes_in_metabolic_pathways_GenesCorrelatedWith_VO2_FDR_0.05.csv"))
# genes_of_interest <- metabolic_genes$SYMBOL
# 
# # --------------------------
# # 3. Relabel and reorder states
# norm_long <- norm_long %>%
#   mutate(Metabolic_State = recode(Metabolic_State,
#                                   "N" = "Normothermy",
#                                   "T" = "Transition",
#                                   "D" = "Deep Torpor"),
#          Metabolic_State = factor(Metabolic_State,
#                                   levels = c("Normothermy", "Transition", "Deep Torpor")))
# 
# # --------------------------
# # 5. Compute ANOVA p-values per gene
# anova_results <- norm_long %>%
#   group_by(Gene) %>%
#   summarise(
#     p_value = tryCatch({
#       m <- aov(Expression ~ Metabolic_State, data = cur_data())
#       summary(m)[[1]][["Pr(>F)"]][1]
#     }, error = function(e) NA_real_),
#     .groups = "drop"
#   ) %>%
#   mutate(p_label = paste0("p = ", signif(p_value, 3)))
# 
# # --------------------------
# # 6. Break into chunks of 20 genes
# n_chunks <- ceiling(length(genes_of_interest) / 20)
# 
# for (i in seq_len(n_chunks)) {
#   subset_genes <- genes_of_interest[((i - 1) * 20 + 1):min(i * 20, length(genes_of_interest))]
#   
#   df_plot <- norm_long %>%
#     filter(Gene %in% subset_genes) %>%
#     left_join(anova_results, by="Gene")
#   
#   p <- ggplot(df_plot,
#               aes(x = Metabolic_State, y = Expression, fill = Metabolic_State)) +
#     geom_violin(trim = FALSE, alpha = 0.6) +
#     geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
#     facet_wrap(~Gene + p_label, scales = "free_y") +
#     scale_fill_manual(values = c("Normothermy" = "#1b9e77",
#                                  "Transition" = "#d95f02",
#                                  "Deep Torpor" = "#7570b3")) +
#     theme_bw(base_size = 14) +
#     theme(axis.text.x = element_text(angle = 30, hjust = 1))
#   
#   print(p)  
#   ggsave(
#     filename = here::here("Outputs/WGCNA/Gut3", paste0("metabolic_genes_chunk_VO2_", i, ".pdf")),
#     plot = p,
#     width = 14,
#     height = 10
#   )
# }

###### Surf_Temp correlated #######

# Set the parent directory containing Module* folders
parent_dir <- here("Outputs/WGCNA", "Gut3")

# Identify Module folders
module_dirs <- list.dirs(
  path = parent_dir,
  full.names = TRUE,
  recursive = FALSE
) %>%
  keep(~ str_detect(basename(.x), "^Module"))

# Read and compile lm_results_mean_VO2.csv across modules
compiled_df <- module_dirs %>%
  map_dfr(function(dir) {
    
    file_path <- file.path(dir, "lm_results_Surf_Temp.csv")
    
    if (!file.exists(file_path)) {
      return(NULL)
    }
    
    read_csv(file_path, show_col_types = FALSE) %>%
      mutate(
        Module = basename(dir)
      )
  })

# Sanity check
stopifnot(nrow(compiled_df) > 0)

# Filter significant results
sig_df <- compiled_df %>%
  filter(p.adj < 0.05)

# Optional: arrange for readability
sig_df <- sig_df %>%
  arrange(p.adj)

sig_df

write.csv(sig_df,here("Outputs/WGCNA/Gut3/genes_in_metabolic_pathways_GenesCorrelatedWith_Surf_Temp_FDR_0.05.csv"))

# Suppose norm.counts is your data frame:
# rows = samples, cols = genes
# rownames(norm.counts) = sample IDs like AS5, AS6, etc.

library(tidyverse)

# --------------------------
# 1. Example: your normalized counts
# rows = samples, cols = genes
#rownames(norm.counts) = sample IDs like "AS5", "AS6", etc.

# --------------------------
# 2. Your sample metadata (replace with your actual mapping!)
metadata <- data.frame(
  SampleID = rownames(norm.counts),
  Metabolic_State = c("N","T","D")[sample(1:3, nrow(norm.counts), replace=TRUE)]
)

# --------------------------
# 1. Reshape counts + join metadata
norm_long <- norm.counts %>%
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to="Gene", values_to="Expression") %>%
  left_join(metadata, by="SampleID")

# --------------------------
# 2. Genes of interest
metabolic_genes <- read_csv(here("Outputs/WGCNA/Gut3/genes_in_metabolic_pathways_GenesCorrelatedWith_Surf_Temp_FDR_0.05.csv"))
genes_of_interest <- metabolic_genes$Gene

# --------------------------
# 3. Relabel and reorder states
norm_long <- norm_long %>%
  mutate(Metabolic_State = recode(Metabolic_State,
                                  "N" = "Normothermy",
                                  "T" = "Transition",
                                  "D" = "Deep Torpor"),
         Metabolic_State = factor(Metabolic_State,
                                  levels = c("Normothermy", "Transition", "Deep Torpor")))


# --------------------------
# 5. Compute ANOVA p-values per gene
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

# --------------------------
# 6. Break into chunks of 20 genes
n_chunks <- ceiling(length(genes_of_interest) / 20)

for (i in seq_len(n_chunks)) {
  subset_genes <- genes_of_interest[((i - 1) * 20 + 1):min(i * 20, length(genes_of_interest))]
  
  df_plot <- norm_long %>%
    filter(Gene %in% subset_genes) %>%
    left_join(anova_results, by="Gene")
  
  p <- ggplot(df_plot,
              aes(x = Metabolic_State, y = Expression, fill = Metabolic_State)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    facet_wrap(~Gene + p_label, scales = "free_y") +
    scale_fill_manual(values = c("Normothermy" = "goldenrod",
                                 "Transition" = "midnightblue",
                                 "Deep Torpor" = "purple")) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  print(p)  
  ggsave(
    filename = here::here("Outputs/WGCNA/Gut3", paste0("metabolic_genes_chunk_Surf_Temp_", i, ".pdf")),
    plot = p,
    width = 14,
    height = 10
  )
}


##### Module membership (How does the eigengene expression correlate with gene expression in the module?)

# ==============================================================
# 0. Load libraries and options
# ==============================================================
library(WGCNA)
options(stringsAsFactors = FALSE)

# ==============================================================
# 1. Prepare your data
datExpr<-norm.counts.filtered
MEs_full <- module_eigengenes
moduleColors
# ==============================================================
# Assuming you already have:
# datExpr               : normalized counts matrix (Samples x Genes)
# MEs_full              : module eigengenes data.frame or matrix
# moduleColors          : vector of module colors for each gene

# Remove unwanted columns from MEs (grey + any non-ME like SampleID)
MEs <- MEs_full[, !(colnames(MEs_full) %in% c("MEgrey", "SampleID")), drop = FALSE]

# ==============================================================
# 2. Calculate full module membership (kME) and p-values
# ==============================================================
geneModuleMembership <- cor(datExpr, MEs, use = "p")
MMPvalue <- corPvalueStudent(geneModuleMembership, nrow(datExpr))

# Convert matrices to data.frames and explicitly name columns
geneModuleMembership <- as.data.frame(geneModuleMembership)
MMPvalue <- as.data.frame(MMPvalue)

colnames(geneModuleMembership) <- paste0("kME.", colnames(MEs))
colnames(MMPvalue) <- paste0("p.kME.", colnames(MEs))

# ==============================================================
# 3. Compile all results into a single data.frame
# ==============================================================
geneInfo <- data.frame(
  GeneID = colnames(datExpr),
  Assigned_Module = moduleColors,
  geneModuleMembership,
  MMPvalue,
  row.names = colnames(datExpr),
  check.names = FALSE
)

# Verify kME columns exist
cat("Columns in geneInfo starting with 'kME':\n")
print(grep("^kME", names(geneInfo), value = TRUE))

# ==============================================================
# 4. Extract intra-modular correlation (kME_Assigned)
# ==============================================================
geneInfo$kME_Assigned <- NA
geneInfo$p_kME_Assigned <- NA

# Prepare module color list
ME_colors_available <- gsub("^ME", "", colnames(MEs))
ME_colors_available <- tolower(ME_colors_available)
geneInfo$Assigned_Module <- tolower(geneInfo$Assigned_Module)

for (module_color in ME_colors_available) {
  
  kME_col_name <- paste0("kME.ME", module_color)
  p_kME_col_name <- paste0("p.kME.ME", module_color)
  
  # Check if columns exist
  if (!(kME_col_name %in% names(geneInfo))) {
    warning(paste("Column", kME_col_name, "not found! Skipping module:", module_color))
    next
  }
  
  # Identify genes in this module
  is_gene_in_module <- geneInfo$Assigned_Module == module_color
  if (sum(is_gene_in_module) == 0) next
  
  # Fill intra-modular kME and p-values
  geneInfo$kME_Assigned[is_gene_in_module] <- geneInfo[[kME_col_name]][is_gene_in_module]
  geneInfo$p_kME_Assigned[is_gene_in_module] <- geneInfo[[p_kME_col_name]][is_gene_in_module]
}

# ==============================================================
# 5. Generate final report for target modules
# ==============================================================
TARGET_MODULES <- c("magenta")

final_report <- geneInfo[
  geneInfo$Assigned_Module %in% TARGET_MODULES,
  c("GeneID", "Assigned_Module", "kME_Assigned", "p_kME_Assigned")
]

# Sort by module and descending kME
final_report_sorted <- final_report[order(
  final_report$Assigned_Module,
  -final_report$kME_Assigned
), ]

# Display top 10 as sanity check
cat("\n--- Top 10 Genes by Intra-Modular kME in Target Modules ---\n")
print(head(final_report_sorted, 10))

write.csv(final_report_sorted,here("Outputs/WGCNA/Gut3/IntraModuleMembership_Gut3.csv"))

######### Sample heatmap for all genes in one module across categorical trait states #######

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(here)

# --------------------------
# Step 1: Pick module
# --------------------------
TARGET_MODULE <- "magenta"

# --------------------------
# Step 2: Get genes in module
# --------------------------
genes_in_module <- names(bwnet$colors)[
  bwnet$colors == tolower(TARGET_MODULE)
]

stopifnot(length(genes_in_module) > 0)

# --------------------------
# Step 3: Subset normalized expression (genes x samples)
# --------------------------
expr_module <- t(
  norm.counts.filtered[, genes_in_module, drop = FALSE]
)

# --------------------------
# Step 4: Z-score per gene
# --------------------------
expr_module_scaled <- t(scale(t(expr_module)))

# --------------------------
# Step 5: Prepare metadata
# --------------------------
metadata_plot <- metadata %>%
  dplyr::select(SampleID, Metabolic_State) %>%
  mutate(
    Metabolic_State = recode(
      Metabolic_State,
      "N" = "Normothermy",
      "T" = "Transition",
      "D" = "Deep Torpor"
    ),
    Metabolic_State = factor(
      Metabolic_State,
      levels = c("Normothermy", "Transition", "Deep Torpor")
    )
  ) %>%
  column_to_rownames("SampleID")

# --------------------------
# Step 6: Align samples
# --------------------------
common_samples <- intersect(
  colnames(expr_module_scaled),
  rownames(metadata_plot)
)

expr_module_scaled <- expr_module_scaled[, common_samples, drop = FALSE]
metadata_plot      <- metadata_plot[common_samples, , drop = FALSE]

stopifnot(
  identical(colnames(expr_module_scaled), rownames(metadata_plot))
)

# --------------------------
# Step 7: Order columns by metabolic state
# --------------------------
sample_order <- order(metadata_plot$Metabolic_State)

expr_module_scaled <- expr_module_scaled[, sample_order, drop = FALSE]
metadata_plot      <- metadata_plot[sample_order, , drop = FALSE]

# --------------------------
# Diverging color scale (z-scored data)
# --------------------------
heatmap_colors <- colorRampPalette(
  c("#2166ac", "white", "#b2182b")
)(101)


# --------------------------
# Step 9: Annotation colors
# --------------------------
annotation_colors <- list(
  Metabolic_State = c(
    "Normothermy" = "#1b9e77",
    "Transition"  = "#d95f02",
    "Deep Torpor" = "#7570b3"
  )
)

# --------------------------
# Step 10: Plot heatmap
# --------------------------
pheatmap(
  expr_module_scaled,
  color             = heatmap_colors,
  annotation_col    = metadata_plot,
  annotation_colors = annotation_colors,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  fontsize          = 12,
  main              = paste("Heatmap of", TARGET_MODULE, "module genes"),
  filename = here::here(
    "Outputs/WGCNA/Gut3",
    paste0("Heatmap_", TARGET_MODULE, ".pdf")
  ),
  width  = 10,
  height = 12
)
