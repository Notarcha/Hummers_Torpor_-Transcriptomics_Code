# ==========================================================================
# SCREEN ONE OR MORE KEGG PATHWAYS (hsa identifiers) against YOUR DATASET and YOUR CUSTOM PANEL
#
#  1. pulls the KEGG gene list(s) (HGNC symbols from the hsa Entrez IDs; no alias forcing);
#     with several pathways the genes are combined and each gene records which pathway(s) it is in
#  2. marks which KEGG genes are rows of your rlog matrix
#  3. compares with your custom groups:
#       KEGG genes you did not list      -> candidate additions (those in your dataset matter)
#       custom genes not in KEGG         -> would be labelled "non-KEGG" additions
#       custom genes not in the dataset  -> cannot be analysed anyway
#  4. flags genes that also belong to "context" pathways (e.g. mTOR, OXPHOS), to see where
#     your custom genes sit when they are not in the main pathway(s)
#
# To reuse for another pathway: edit ONLY the block "EDIT THIS BLOCK".
# All CSVs are written to main_fig_dir.
# ==========================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(KEGGREST)
  library(org.Hs.eg.db); library(AnnotationDbi)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
summarise <- dplyr::summarise; arrange <- dplyr::arrange; count <- dplyr::count

# ==========================================================================
# EDIT THIS BLOCK
# ==========================================================================
FILE_PREFIX  <- "PPAR"
# Main pathway: hsa03320 PPAR signaling pathway (check the name printed when the script runs)
KEGG_IDS     <- c("hsa03320")
# Context pathways: NOT part of the set, only used to show where other genes belong
#   Fatty acid metabolism, Fatty acid degradation, Biosynthesis of unsaturated FA, Fatty acid biosynthesis,
#   Synthesis/degradation of ketone bodies, Glycolysis/gluconeogenesis, AMPK, FoxO, Adipocytokine,
#   Thermogenesis, Cholesterol metabolism, Peroxisome
overlap_pathways <- c("hsa01212", "hsa00071", "hsa01040", "hsa00061", "hsa00072", "hsa00010",
                      "hsa04152", "hsa04068", "hsa04920", "hsa04714", "hsa04979", "hsa04146")
main_fig_dir <- "/Users/notarcha/Desktop/PNAS_Code/NewOutputs"

custom_groups <- list(
  "PPAR Receptors & Co-regulators" = c("PPARA", "PPARD", "PPARG", "RXRA", "NR1H3"),
  "Fatty Acid Transport & Binding" = c("CD36", "FABP1", "FABP3", "SLC27A4", "DBI", "SCP2"),
  "Fatty Acid Activation & Oxidation" = c("ACSL1", "CPT1A", "CPT1B", "CPT2", "ACADL", "ACOX1", "EHHADH", "ACAA1"),
  "Lipid Storage & Adipogenesis" = c("PLIN1", "PLIN2", "LPL", "ANGPTL4", "ADIPOQ"),
  "Lipoprotein & Cholesterol Metabolism" = c("APOA1", "APOA5", "APOC3", "PLTP"),
  "Glucose, Ketone & Unsaturated FA Metabolism" = c("PCK1", "GK", "ME1", "HMGCS2", "SCD", "FADS2")
)
# ==========================================================================

dir.create(main_fig_dir, showWarnings = FALSE, recursive = TRUE)
tag <- if (length(KEGG_IDS) == 1) KEGG_IDS else FILE_PREFIX
out <- function(suffix) file.path(main_fig_dir, paste0("KEGG_", tag, "_", suffix))

# ---- 1. KEGG gene list(s) ----------------------------------------------------------------
pull_one <- function(id) {
  pw <- keggGet(id)[[1]]
  cat("Pathway:", id, "-", pw$NAME, "|", length(pw$GENE) / 2, "genes\n")
  g <- pw$GENE; entrez <- g[c(TRUE, FALSE)]; desc <- g[c(FALSE, TRUE)]
  tibble(Entrez = entrez, KEGG_symbol = sub(";.*", "", desc),
         Description = sub("^[^;]*; ?", "", desc) %>% str_remove_all("\\[.*?\\]") %>% str_trim(),
         EC = str_extract(desc, "(?<=\\[EC:)[^\\]]+"), KO = str_extract(desc, "(?<=\\[KO:)[^\\]]+"),
         source_pathway = id)
}
kegg_raw <- map_dfr(KEGG_IDS, pull_one)
ue   <- unique(kegg_raw$Entrez)
hs_u <- mapIds(org.Hs.eg.db, keys = ue, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
kegg_raw$Symbol <- unname(hs_u[kegg_raw$Entrez])
kegg_raw$Symbol <- ifelse(is.na(kegg_raw$Symbol), kegg_raw$KEGG_symbol, kegg_raw$Symbol)

kegg_tbl <- kegg_raw %>% group_by(Symbol) %>%
  summarise(Entrez = dplyr::first(Entrez), KEGG_symbol = dplyr::first(KEGG_symbol),
            Description = dplyr::first(Description), EC = dplyr::first(EC), KO = dplyr::first(KO),
            source_pathway = paste(sort(unique(source_pathway)), collapse = "; "), .groups = "drop")

# ---- 2. Dataset presence (exact hsa symbol; aliases are reported but NOT used) ---------------
data_genes <- rownames(read.csv(here("Results", "rlog_normalized_counts_STAR.csv"),
                                row.names = 1, check.names = FALSE))
kegg_tbl <- kegg_tbl %>% mutate(in_dataset = Symbol %in% data_genes)

miss <- kegg_tbl %>% filter(!in_dataset)
if (nrow(miss) > 0) {
  al <- AnnotationDbi::select(org.Hs.eg.db, keys = miss$Entrez, columns = "ALIAS", keytype = "ENTREZID") %>%
    filter(ALIAS %in% data_genes) %>% group_by(ENTREZID) %>%
    summarise(alias_in_dataset = paste(unique(ALIAS), collapse = "; "), .groups = "drop")
  kegg_tbl <- kegg_tbl %>% left_join(al, by = c("Entrez" = "ENTREZID"))
} else kegg_tbl$alias_in_dataset <- NA_character_

# ---- 3. Context pathways ------------------------------------------------------------------------
for (p in overlap_pathways) cat("Context pathway:", p, "-", tryCatch(keggGet(p)[[1]]$NAME, error = function(e) "?"), "\n")
other_ids <- map(setNames(overlap_pathways, overlap_pathways), function(p)
  tryCatch(sub("hsa:", "", unique(keggLink("hsa", paste0("path:", p)))), error = function(e) character(0)))
ctx <- function(entrez) paste(names(other_ids)[map_lgl(other_ids, ~ entrez %in% .x)], collapse = "; ")
kegg_tbl$also_in_pathways <- map_chr(kegg_tbl$Entrez, ctx)

# ---- 4. Compare with the custom panel -------------------------------------------------------------
custom_long <- imap_dfr(custom_groups, ~ tibble(custom_group = .y, Gene = .x))
custom_syms <- unique(custom_long$Gene)
# where do custom genes that are NOT in the main pathways sit? (look up their Entrez, flag context pathways)
cust_entrez <- suppressWarnings(AnnotationDbi::select(org.Hs.eg.db, keys = custom_syms, columns = "ENTREZID", keytype = "SYMBOL")) %>%
  distinct(SYMBOL, .keep_all = TRUE)
cust_ctx <- tibble(Gene = custom_syms, EntrezC = cust_entrez$ENTREZID[match(custom_syms, cust_entrez$SYMBOL)]) %>%
  mutate(context_pathways = map_chr(EntrezC, ~ if (is.na(.x)) "" else ctx(.x))) %>% select(Gene, context_pathways)

cmp <- tibble(Gene = union(kegg_tbl$Symbol, custom_long$Gene)) %>%
  left_join(kegg_tbl %>% select(Gene = Symbol, source_pathway, Description, EC, KO, alias_in_dataset, also_in_pathways, Entrez),
            by = "Gene") %>%
  left_join(custom_long %>% group_by(Gene) %>% summarise(custom_group = paste(custom_group, collapse = "; "), .groups = "drop"),
            by = "Gene") %>%
  left_join(cust_ctx, by = "Gene") %>%
  mutate(in_KEGG = !is.na(Entrez), in_custom = !is.na(custom_group), in_dataset = Gene %in% data_genes,
         status = case_when(in_KEGG & in_custom ~ "in both", in_KEGG ~ "KEGG only", TRUE ~ "custom only"),
         context_pathways = ifelse(is.na(context_pathways) | context_pathways == "", also_in_pathways, context_pathways)) %>%
  select(Gene, status, in_KEGG, in_custom, in_dataset, custom_group, source_pathway, Description, EC, KO,
         context_pathways, alias_in_dataset, Entrez) %>%
  arrange(status, custom_group, Gene)

grp_sum <- custom_long %>% left_join(cmp %>% select(Gene, in_KEGG, in_dataset), by = "Gene") %>%
  group_by(custom_group) %>%
  summarise(n_custom = n(), n_in_KEGG = sum(in_KEGG), n_in_dataset = sum(in_dataset),
            n_in_KEGG_and_dataset = sum(in_KEGG & in_dataset),
            custom_not_in_KEGG = paste(Gene[!in_KEGG], collapse = "; "),
            custom_not_in_dataset = paste(Gene[!in_dataset], collapse = "; "), .groups = "drop")

candidates <- cmp %>% filter(status == "KEGG only", in_dataset) %>%
  select(Gene, source_pathway, Description, context_pathways)

pathway_cov <- kegg_raw %>% mutate(in_dataset = Symbol %in% data_genes) %>%
  group_by(source_pathway) %>% summarise(KEGG_genes = n_distinct(Symbol), in_dataset = n_distinct(Symbol[in_dataset]), .groups = "drop")

# ---- 5. Report -----------------------------------------------------------------------------------------
cat("\n=== COVERAGE per pathway (exact hsa symbol) ===\n"); print(as.data.frame(pathway_cov))
cat("\nUnion of KEGG genes:", nrow(kegg_tbl), "| in your dataset:", sum(kegg_tbl$in_dataset), "\n")
cat("Custom genes:", nrow(custom_long), "| in KEGG:", sum(custom_long$Gene %in% kegg_tbl$Symbol),
    "| in dataset:", sum(custom_long$Gene %in% data_genes), "\n")
cat("\n=== CUSTOM GROUPS vs KEGG and DATASET ===\n"); print(as.data.frame(grp_sum), right = FALSE)
cat("\n=== CUSTOM genes NOT in the main KEGG pathway(s) (context pathways shown) ===\n")
print(cmp %>% filter(status == "custom only") %>% select(Gene, custom_group, in_dataset, context_pathways), n = Inf)
cat("\n=== KEGG genes you did NOT list that ARE in your dataset (", nrow(candidates), ") ===\n", sep = "")
print(candidates, n = 60, width = Inf)

write.csv(kegg_tbl,    out("gene_list_check.csv"), row.names = FALSE)
write.csv(cmp,         out("custom_vs_KEGG_comparison.csv"), row.names = FALSE)
write.csv(grp_sum,     out("custom_group_summary.csv"), row.names = FALSE)
write.csv(candidates,  out("candidate_additions_in_dataset.csv"), row.names = FALSE)
write.csv(pathway_cov, out("pathway_coverage.csv"), row.names = FALSE)
cat("\nSaved 5 CSVs to:", main_fig_dir, "\n")



# ==========================================================================
# PPAR SIGNALLING (v6): pre-specified gene-set tests, variance-aware, honestly labelled
#
# WHAT CHANGED vs v1
#  1. Two families of set tests, each labelled for what it asks:
#       COMPETITIVE     (camera)         "is the set shifted MORE than the other genes?"
#                         - per tissue
#                         - pooled across tissues: one limma-style model with a tissue term,
#                           bird random intercept (duplicateCorrelation) and unequal variances
#                           (arrayWeights by tissue x state), implemented by whitening the data
#                           so that camera()/mroast() can be run on the pooled torpor effect
#       SELF-CONTAINED  (set score)      "does the set score shift AT ALL?"
#                         - per tissue: Welch t-test
#                         - pooled: nlme::lme, bird random intercept, tissue x state variances
#     Two heatmaps are drawn, one per family, so no column mixes the two.
#  2. Heatmap colour = a VARIANCE-AWARE statistic (sets: Hedges g; genes: moderated t).
#     Printed in each cell: the mean difference D - N (in Z units) + significance symbols.
#  3. Master audit CSV: for every set x analysis, which genes were defined, found, used by
#     camera, dropped (and why), plus the p/q values. Plus a gene-membership CSV.
#  4. Unequal variances in the pooled models (see 1).
#
# To reuse for another pathway: edit ONLY the block "PATHWAY DEFINITION" (and main_fig_dir).
# ==========================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(limma)
  library(pheatmap); library(KEGGREST); library(org.Hs.eg.db); library(AnnotationDbi)
})
# nlme is used via nlme:: (ships with R); lme4 is no longer needed.
options(pillar.sigfig = 3)
select <- dplyr::select;       rename   <- dplyr::rename;   filter    <- dplyr::filter
mutate <- dplyr::mutate;       transmute <- dplyr::transmute; summarise <- dplyr::summarise
count  <- dplyr::count;        arrange  <- dplyr::arrange;  recode    <- dplyr::recode

# ---- styling ----------------------------------------------------------------
my_colors <- colorRampPalette(c("#b35806", "#f1a340", "#fee0b6", "#f7f7f7",
                                "#d8daeb", "#998ec3", "#542788"))(100)
tissue_order <- c("Proximal Gut", "Medial Gut", "Distal Gut", "Lungs", "Heart", "Liver", "Pectoral")

# ---- settings ---------------------------------------------------------------
min_n_per_group <- 3
fdr_cut         <- 0.05
run_pathview    <- FALSE   # KEGG maps are not informative here; set TRUE to draw hsa03320
set_g_lim       <- 2      # colour scale limit for set-level Hedges g (values beyond are clamped in colour only)
gene_t_lim      <- 5      # colour scale limit for gene-level moderated t
# Which groups get their own residual variance in the POOLED models.
# c("Tissue", "Metabolic_State") = one variance per tissue x state (recommended);
# "Tissue" = one variance per tissue.
var_group_by    <- c("Tissue", "Metabolic_State")
pooled_name     <- "ALL tissues (pooled)"
set.seed(1)

# ==========================================================================
# PATHWAY DEFINITION  (the only block that changes between pathways)
# Gene sets = (1) the whole KEGG pathway, (2) your curated functional groups, and
# (3) optionally the KEGG modules of the pathway. Every set is tested (camera etc.) and reported.
# ==========================================================================
FILE_PREFIX   <- "PPAR"
KEGG_IDS      <- c("hsa03320")   # PPAR signaling pathway
KEGG_ID       <- KEGG_IDS[1]
KEGG_SET_NAME <- "KEGG PPAR signaling pathway (hsa03320)"   # the whole-pathway KEGG set

# (2) Curated functional groups: defined from biology / gene nomenclature BEFORE looking at any
# results. They are all tested and reported (nothing is dropped after seeing results).
custom_groups <- list(
  # Curated panel (all 34 genes are in KEGG hsa03320) + additions from the audit. All members are KEGG genes.
  # Genes not in the count matrix are dropped automatically at the test stage.
  "PPAR Receptors & Co-regulators" = c("PPARA", "PPARD", "PPARG", "RXRA", "RXRG", "NR1H3"),
  "Fatty Acid Transport & Binding" = c("CD36", "FABP1", "FABP3", "FABP6", "FABP7", "SLC27A2", "SLC27A4", "SLC27A6", "DBI", "SCP2"),
  "Fatty Acid Activation & Oxidation" = c("ACSL1", "ACSL3", "ACSL4", "ACSL5", "ACSL6", "ACSBG1", "ACSBG2", "CPT1A", "CPT1B", "CPT2", "ACADL", "ACOX1", "ACOX2", "ACOX3", "EHHADH", "ACAA1"),
  "Lipid Storage & Adipogenesis" = c("PLIN1", "PLIN2", "LPL", "ANGPTL4", "ADIPOQ"),
  "Lipoprotein & Cholesterol Metabolism" = c("APOA1", "APOA5", "APOC3", "PLTP", "HMGCS1"),
  "Glucose, Ketone & Unsaturated FA Metabolism" = c("PCK1", "GK", "ME1", "ME3", "HMGCS2", "SCD", "SCD5", "FADS2")
)
# TRUE = keep only group genes whose symbol is in the KEGG pathway (exact hsa symbols; nothing is forced
# by alias or ortholog). Groups named in nonkegg_groups are exempt (deliberately added, labelled non-KEGG).
restrict_custom_to_kegg <- TRUE
nonkegg_groups          <- character(0)

# (2b) Groups taken straight from the KEGG pathway gene list by symbol pattern (regex on hsa symbols)
kegg_pattern_groups     <- list()

# (2c) Extra whole-pathway sets = this KEGG pathway MINUS the genes of other KEGG pathways
# (rule-based: no hand-picking). Used to avoid re-testing genes analysed in another pathway.
kegg_exclusion_sets     <- list()

# (3) OPTIONAL extra sets = KEGG MODULES linked to the pathway (names and genes come from KEGG):
#   NULL              -> none (default)
#   "auto"            -> every KEGG module linked to the pathway
#   c("M00142", ...)  -> only these module IDs
# Note: KEGG modules can be nested (overlapping genes) and add to the number of tests.
KEGG_MODULES  <- NULL

# Genes you want FLAGGED in tables/heatmap (e.g. nominal hits from a Results table). They are
# corrected as part of the full family, never privileged. Use character(0) if there are none.
flag_label    <- "Flagged genes"
posthoc_genes <- character(0)

# umbrella set: tested, but not used to assign genes to a row-block in the gene heatmap
aggregate_sets  <- c(KEGG_SET_NAME, names(kegg_exclusion_sets))
kegg_only_label <- "KEGG only (not in a curated group)"

define_sets <- function(all_genes) {
  found <- function(g) intersect(g, all_genes)
  entrez2sym <- function(ids) {
    s <- mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    unique(na.omit(unname(s)))
  }
  pw_ids   <- character(0)
  kegg_all <- tryCatch({
    pw_ids <- unique(unlist(lapply(KEGG_IDS, function(p) sub("hsa:", "", unique(keggLink("hsa", paste0("path:", p)))))))
    entrez2sym(pw_ids)
  }, error = function(e) {
    message("KEGG pathway lookup failed (", conditionMessage(e), "); the whole-pathway set is skipped.")
    character(0)
  })
  
  # KEGG modules belonging to the pathway (names and genes both come from KEGG)
  mod_groups <- list()
  if (!is.null(KEGG_MODULES) && length(pw_ids)) {
    mod_groups <- tryCatch({
      mods <- if (identical(KEGG_MODULES, "auto")) {
        ml <- keggLink("module", paste0("path:", KEGG_ID))
        unique(sub("^md:(hsa_)?", "", ml))
      } else sub("^md:(hsa_)?", "", KEGG_MODULES)
      mod_names <- tryCatch(keggList("module"), error = function(e) character(0))
      out <- list()
      for (m in mods) {
        ids <- character(0)
        for (src in c(paste0("md:hsa_", m), paste0("md:", m))) {     # organism-specific id first
          ids <- tryCatch(keggLink("hsa", src), error = function(e) character(0))
          if (length(ids)) break
        }
        ids <- intersect(sub("hsa:", "", unique(ids)), pw_ids)        # keep genes that are in the pathway
        if (!length(ids)) next
        nm <- unname(mod_names[paste0("md:", m)]); if (length(nm) == 0 || is.na(nm)) nm <- "module"
        out[[paste0(m, ": ", nm)]] <- entrez2sym(ids)
      }
      out
    }, error = function(e) {
      message("KEGG module lookup failed (", conditionMessage(e), "); only the whole-pathway set is used.")
      list()
    })
    cat("\nKEGG modules used as sub-groups:", if (length(mod_groups)) paste(names(mod_groups), collapse = " | ") else "none", "\n")
  }
  
  cg <- custom_groups
  if (restrict_custom_to_kegg) {
    rs      <- setdiff(names(cg), nonkegg_groups)
    dropped <- setdiff(unlist(cg[rs]), kegg_all)
    cg[rs]  <- lapply(cg[rs], function(g) intersect(g, kegg_all))
    cat("
Curated genes dropped because their symbol is not in KEGG ", KEGG_ID, " (", length(dropped), "): ",
        paste(dropped, collapse = ", "), "
", sep = "")
  }
  pat_groups <- lapply(kegg_pattern_groups, function(p) grep(p, kegg_all, value = TRUE))
  groups  <- c(cg, pat_groups, mod_groups)
  excl_sets <- lapply(kegg_exclusion_sets, function(p) {
    other <- tryCatch(entrez2sym(sub("hsa:", "", unique(keggLink("hsa", paste0("path:", p))))),
                      error = function(e) character(0))
    if (!length(other)) message("Could not pull ", p, "; the exclusion set will equal the whole pathway.")
    setdiff(kegg_all, other)
  })
  defined <- c(setNames(list(kegg_all), KEGG_SET_NAME), excl_sets, groups)
  sets    <- lapply(defined, found)
  attr(sets, "defined")   <- defined                 # genes the set was MEANT to contain
  attr(sets, "n_defined") <- lengths(defined)
  sets
}
# ==========================================================================

main_fig_dir <- "/Users/notarcha/Desktop/PNAS_Code/NewOutputs"
dir.create(main_fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. LOAD + CLEAN (N vs D only) -------------------------------------------
RNASeq_metadata  <- read.csv(here("Results", "RNASeq_metadata.csv"), row.names = 1)
Normalizedcounts <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"),
                             row.names = 1, check.names = FALSE)

metadata <- RNASeq_metadata %>% rownames_to_column("SampleID")
metadata$SampleID <- metadata$SampleID %>% str_trim() %>% str_replace_all("-", "_")
colnames(Normalizedcounts) <- colnames(Normalizedcounts) %>% str_trim() %>% str_replace_all("-", "_")
Normalizedcounts <- as.matrix(Normalizedcounts)
storage.mode(Normalizedcounts) <- "double"
stopifnot("BirdID" %in% colnames(metadata))

common_samples   <- intersect(metadata$SampleID, colnames(Normalizedcounts))
metadata         <- metadata %>% filter(SampleID %in% common_samples, Metabolic_State != "T")
Normalizedcounts <- Normalizedcounts[, metadata$SampleID]

metadata <- metadata %>%
  mutate(Tissue = recode(Tissue, "Gut1" = "Proximal Gut", "Gut2" = "Medial Gut",
                         "Gut3" = "Distal Gut", "Pect" = "Pectoral"),
         Tissue = factor(Tissue, levels = tissue_order),
         Metabolic_State = factor(Metabolic_State, levels = c("N", "D")),
         BirdID = factor(BirdID)) %>%
  arrange(Metabolic_State, Tissue)
Normalizedcounts <- Normalizedcounts[, metadata$SampleID]

cat("\nSamples per tissue x state (N vs D only):\n")
print(table(metadata$Tissue, metadata$Metabolic_State))

tissues_tab <- metadata %>% count(Tissue, Metabolic_State) %>%
  pivot_wider(names_from = Metabolic_State, values_from = n, values_fill = 0) %>%
  mutate(tested = N >= min_n_per_group & D >= min_n_per_group)
tested_tissues <- as.character(tissues_tab$Tissue[tissues_tab$tested])
untested_tissues <- setdiff(tissue_order, tested_tissues)
cat("Tissues tested:", paste(tested_tissues, collapse = ", "), "\n")
if (length(untested_tissues))
  cat("NOT tested (n <", min_n_per_group, "in a state):", paste(untested_tissues, collapse = ", "), "\n")

# ---- 2. GENE SETS --------------------------------------------------------------
sets_raw     <- define_sets(rownames(Normalizedcounts))
defined_list <- attr(sets_raw, "defined")
sets         <- Filter(function(g) length(g) >= 3, sets_raw)
attr(sets, "defined") <- NULL; attr(sets, "n_defined") <- NULL
cat("\nPre-specified gene sets (genes found in the data):\n")
print(tibble(Set = names(sets_raw), n_defined = unname(lengths(defined_list)[names(sets_raw)]),
             n_found_in_data = lengths(sets_raw), kept = names(sets_raw) %in% names(sets)), width = Inf)
set_genes <- unique(unlist(sets_raw))
cat(flag_label, "outside every set (not in the tested family):",
    paste(setdiff(posthoc_genes, set_genes), collapse = ", "), "\n")

# ---- 3. PER-TISSUE Z-SCORES + CHECK ----------------------------------------------
z_scale <- function(m) {
  m  <- as.matrix(m)
  mu <- rowMeans(m); s <- apply(m, 1, sd)
  z  <- (m - mu) / s
  z[!is.finite(s) | s == 0, ] <- 0
  z
}
if (length(set_genes) == 0)
  stop("No genes were found for any set. First rownames: ", paste(head(rownames(Normalizedcounts)), collapse = ", "))
zmat <- matrix(NA_real_, nrow = length(set_genes), ncol = nrow(metadata),
               dimnames = list(set_genes, metadata$SampleID))
for (tis in levels(metadata$Tissue)) {
  ids <- metadata$SampleID[metadata$Tissue == tis]
  if (length(ids) < 2) next
  zmat[, ids] <- z_scale(Normalizedcounts[set_genes, ids, drop = FALSE])
}
z_check <- map_dfr(levels(metadata$Tissue), function(tis) {
  ids <- metadata$SampleID[metadata$Tissue == tis]
  m <- zmat[, ids, drop = FALSE]; s <- apply(m, 1, sd)
  tibble(Tissue = tis, n_samples = length(ids), max_abs_gene_mean = max(abs(rowMeans(m))))
})
stopifnot(all(z_check$max_abs_gene_mean < 1e-8))

# ---- helpers -------------------------------------------------------------------
fmt_genes <- function(x) if (length(x)) paste(sort(x), collapse = "; ") else ""

# Hedges g (D - N) with approximate 95% CI, and pooled within-group SD
hedges <- function(d_vals, n_vals) {
  nD <- length(d_vals); nN <- length(n_vals); df <- nD + nN - 2
  sp <- sqrt(((nD - 1) * var(d_vals) + (nN - 1) * var(n_vals)) / df)
  g  <- (1 - 3 / (4 * df - 1)) * (mean(d_vals) - mean(n_vals)) / sp
  se <- sqrt((nD + nN) / (nD * nN) + g^2 / (2 * (nD + nN)))
  tibble(g = g, g_lo = g - 1.96 * se, g_hi = g + 1.96 * se, sd_pooled = sp)
}

# Audit rows (static part): one row per set for a given analysis
audit_static <- function(analysis, Yrows, tested) {
  map_dfr(names(sets_raw), function(s) {
    dfn <- defined_list[[s]]; fnd <- sets_raw[[s]]
    used <- if (tested) intersect(fnd, Yrows) else character(0)
    reason <- if (!tested) paste0("tissue has < ", min_n_per_group, " birds in N or D")
    else if (length(fnd) < 3) "< 3 set genes found in data"
    else if (length(used) < 3) "< 3 set genes left after zero-variance filter"
    else ""
    tibble(Set = s, Analysis = analysis,
           n_defined = length(dfn), n_found_in_data = length(fnd), n_used_in_camera = length(used),
           n_background_genes = if (tested) length(Yrows) - length(used) else NA_integer_,
           camera_run = tested && reason == "", reason_not_run = reason,
           genes_defined = fmt_genes(dfn),
           genes_used_in_camera = fmt_genes(used),
           genes_not_in_data = fmt_genes(setdiff(dfn, fnd)),
           genes_dropped_zero_variance = if (tested) fmt_genes(setdiff(fnd, Yrows)) else "")
  })
}

# ---- 4. PER-TISSUE TESTS -----------------------------------------------------------
res_list <- lapply(tested_tissues, function(tis) {
  md  <- metadata %>% filter(Tissue == tis)
  idN <- md$SampleID[md$Metabolic_State == "N"]
  idD <- md$SampleID[md$Metabolic_State == "D"]
  Y   <- as.matrix(Normalizedcounts[, md$SampleID])
  Y   <- Y[apply(Y, 1, var) > 0, , drop = FALSE]
  design <- model.matrix(~ Metabolic_State, data = md)
  
  # gene level: moderated t (variance shrunk toward a prior; uses all genes)
  fit <- eBayes(lmFit(Y, design))
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none") %>% rownames_to_column("Gene")
  gene_res <- tt %>% filter(Gene %in% set_genes) %>%
    transmute(Tissue = tis, Gene, logFC_rlog = logFC, t, p = P.Value)
  gene_res$shift_Z <- rowMeans(zmat[gene_res$Gene, idD, drop = FALSE]) -
    rowMeans(zmat[gene_res$Gene, idN, drop = FALSE])
  
  # COMPETITIVE: camera (inter-gene correlation estimated per set) + mroast 'mixed'
  idx <- lapply(sets, function(g) which(rownames(Y) %in% g))
  idx <- idx[lengths(idx) >= 3]
  cam_raw <- camera(Y, idx, design, contrast = 2, inter.gene.cor = NA)
  cam_res <- data.frame(Set = rownames(cam_raw), Tissue = tis,
                        cam_raw[, setdiff(names(cam_raw), c("PValue", "FDR")), drop = FALSE],
                        p_camera = cam_raw$PValue, stringsAsFactors = FALSE, row.names = NULL)
  mr <- tryCatch(mroast(Y, idx, design, contrast = 2, nrot = 999), error = function(e) NULL)
  cam_res$p_mixed <- if (!is.null(mr) && "PValue.Mixed" %in% names(mr)) mr[cam_res$Set, "PValue.Mixed"] else NA_real_
  
  # SELF-CONTAINED: Welch t on the per-sample set score + Hedges g
  sc_res <- map_dfr(names(sets), function(s) {
    score <- colMeans(zmat[sets[[s]], md$SampleID, drop = FALSE])
    bind_cols(tibble(Tissue = tis, Set = s, nN = length(idN), nD = length(idD),
                     mean_N = mean(score[idN]), mean_D = mean(score[idD]),
                     shift_set_Z = mean(score[idD]) - mean(score[idN]),
                     p_selfcontained = t.test(score[idD], score[idN])$p.value),
              hedges(score[idD], score[idN]))
  })
  list(gene = gene_res, cam = cam_res, sc = sc_res,
       audit = audit_static(tis, rownames(Y), TRUE), yrows = rownames(Y))
})
names(res_list) <- tested_tissues

gene_set_map <- tibble(Gene = set_genes) %>%
  mutate(Set = map_chr(Gene, function(g) {
    base <- sets_raw[setdiff(names(sets_raw), aggregate_sets)]
    hit  <- names(base)[map_lgl(base, ~ g %in% .x)]
    if (length(hit)) hit[1] else kegg_only_label
  }))

# genes that are only in the whole-pathway KEGG list (in no curated group) are left out of the
# gene-level tables/figure; they still count in the whole-pathway set test
gene_tbl <- bind_rows(map(res_list, "gene")) %>%
  left_join(gene_set_map, by = "Gene") %>%
  filter(Set != kegg_only_label) %>%
  mutate(q = p.adjust(p, "BH"), in_posthoc_list = Gene %in% posthoc_genes)
cam_tbl <- bind_rows(map(res_list, "cam")) %>%
  mutate(q_camera = p.adjust(p_camera, "BH"), q_mixed = p.adjust(p_mixed, "BH"))
sc_tbl  <- bind_rows(map(res_list, "sc"))  %>% mutate(q_selfcontained = p.adjust(p_selfcontained, "BH"))

# ---- 4b. POOLED GENE-LEVEL TEST (one test per gene across all tested tissues) -----------
# limma on the full transcriptome: expression ~ Tissue + State, bird random intercept
# (duplicateCorrelation), unequal variances (arrayWeights by tissue x state), moderated t.
# BH is applied across the genes in the sets only (one family of n_set_genes tests).
pooled_gene <- function() {
  md <- metadata %>% filter(Tissue %in% tested_tissues) %>% mutate(Tissue = droplevels(Tissue))
  Y  <- as.matrix(Normalizedcounts[, md$SampleID]); Y <- Y[apply(Y, 1, var) > 0, , drop = FALSE]
  design  <- model.matrix(~ Tissue + Metabolic_State, data = md)
  var_grp <- interaction(md[var_group_by], drop = TRUE)
  aw  <- arrayWeights(Y, design = design, var.group = var_grp)
  awm <- matrix(aw, nrow(Y), length(aw), byrow = TRUE)
  dc  <- duplicateCorrelation(Y, design, block = md$BirdID, weights = awm)
  rho <- dc$consensus.correlation
  if (!is.finite(rho) || rho < 0) rho <- 0
  rho <- min(rho, 0.95)
  fit <- eBayes(lmFit(Y, design, block = md$BirdID, correlation = rho, weights = awm))
  topTable(fit, coef = "Metabolic_StateD", number = Inf, sort.by = "none") %>%
    rownames_to_column("Gene") %>% filter(Gene %in% set_genes) %>%
    transmute(Gene, t_pooled = t, p_pooled = P.Value, logFC_pooled = logFC)
}
gene_pooled_tbl <- pooled_gene() %>% filter(Gene %in% gene_tbl$Gene) %>%
  left_join(gene_tbl %>% group_by(Gene) %>%
              summarise(mean_shift_Z = mean(shift_Z), n_tissues_up = sum(shift_Z > 0),
                        n_tissues = n(), .groups = "drop"), by = "Gene") %>%
  mutate(q_pooled = p.adjust(p_pooled, "BH")) %>%
  left_join(gene_set_map, by = "Gene")
cat("\n=== GENE LEVEL, pooled across tissues (BH over", nrow(gene_pooled_tbl), "genes) ===\n",
    "nominal p < 0.05:", sum(gene_pooled_tbl$p_pooled < 0.05), "| BH q <", fdr_cut, ":",
    sum(gene_pooled_tbl$q_pooled < fdr_cut), "\n")
print(as_tibble(gene_pooled_tbl) %>% arrange(p_pooled) %>%
        select(Gene, Set, t_pooled, p_pooled, q_pooled, mean_shift_Z, n_tissues_up), n = 15, width = Inf)

# ---- 5a. POOLED COMPETITIVE TEST (camera on the pooled torpor effect) ------------------
# Model: expression ~ Tissue + State, with
#   * unequal variances: arrayWeights() with one variance per group in `var_group_by`
#   * bird random intercept: duplicateCorrelation() -> consensus rho
# camera()/mroast() cannot take a block/correlation argument, so the data and design are
# WHITENED (GLS transform: scale by sqrt(weights), then multiply each bird's block by the
# inverse Cholesky factor of its correlation matrix). Ordinary camera on the whitened
# data then tests the GLS torpor effect.
pooled_competitive <- function() {
  md <- metadata %>% filter(Tissue %in% tested_tissues) %>% mutate(Tissue = droplevels(Tissue))
  Y  <- as.matrix(Normalizedcounts[, md$SampleID]); Y <- Y[apply(Y, 1, var) > 0, , drop = FALSE]
  design   <- model.matrix(~ Tissue + Metabolic_State, data = md)
  coef_idx <- which(colnames(design) == "Metabolic_StateD")
  var_grp  <- interaction(md[var_group_by], drop = TRUE)
  
  aw  <- arrayWeights(Y, design = design, var.group = var_grp)
  awm <- matrix(aw, nrow(Y), length(aw), byrow = TRUE)
  dc  <- duplicateCorrelation(Y, design, block = md$BirdID, weights = awm)
  rho <- dc$consensus.correlation
  if (!is.finite(rho) || rho < 0) { message("Pooled bird correlation <= 0; set to 0."); rho <- 0 }
  rho <- min(rho, 0.95)
  cat("\nPooled model: bird consensus correlation =", round(rho, 3), "\n")
  
  Linv <- diag(nrow(md))
  for (b in unique(as.character(md$BirdID))) {
    i <- which(as.character(md$BirdID) == b); k <- length(i)
    if (k > 1) {
      R <- (1 - rho) * diag(k) + rho * matrix(1, k, k)
      Linv[i, i] <- solve(t(chol(R)))
    }
  }
  Ys <- sweep(Y, 2, sqrt(aw), "*") %*% t(Linv)
  Xs <- Linv %*% (design * sqrt(aw))
  colnames(Ys) <- md$SampleID; colnames(Xs) <- colnames(design)
  
  idx <- lapply(sets, function(g) which(rownames(Ys) %in% g))
  idx <- idx[lengths(idx) >= 3]
  cam_raw <- camera(Ys, idx, Xs, contrast = coef_idx, inter.gene.cor = NA)
  cam_res <- data.frame(Set = rownames(cam_raw), Tissue = pooled_name,
                        cam_raw[, setdiff(names(cam_raw), c("PValue", "FDR")), drop = FALSE],
                        p_camera = cam_raw$PValue, stringsAsFactors = FALSE, row.names = NULL)
  mr <- tryCatch(mroast(Ys, idx, Xs, contrast = coef_idx, nrot = 999), error = function(e) NULL)
  cam_res$p_mixed <- if (!is.null(mr) && "PValue.Mixed" %in% names(mr)) mr[cam_res$Set, "PValue.Mixed"] else NA_real_
  
  wt <- tibble(group = levels(var_grp), n_samples = as.integer(table(var_grp)),
               array_weight = as.numeric(tapply(aw, var_grp, mean)),
               relative_variance = 1 / array_weight, bird_rho = rho)
  list(cam = cam_res, weights = wt, yrows = rownames(Ys))
}
pc <- pooled_competitive()
cam_pooled <- pc$cam %>% mutate(q_camera = p.adjust(p_camera, "BH"), q_mixed = p.adjust(p_mixed, "BH"))
cat("\nPooled-model variance weights (relative variance > 1 = noisier group, down-weighted):\n")
print(pc$weights, width = Inf)

# ---- 5b. POOLED SELF-CONTAINED TEST (set score, lme, bird intercept, unequal variances) ---
md_t <- metadata %>% filter(Tissue %in% tested_tissues) %>% mutate(Tissue = droplevels(Tissue))

fit_lme <- function(fml, d, het) {
  args <- list(fixed = fml, random = ~ 1 | BirdID, data = d, method = "ML",
               control = nlme::lmeControl(maxIter = 200, msMaxIter = 200, returnObject = TRUE))
  if (het) args$weights <- nlme::varIdent(form = ~ 1 | grp)
  tryCatch(suppressWarnings(do.call(nlme::lme, args)), error = function(e) NULL)
}
fit_pair <- function(f1, f0, d) {
  for (het in c(TRUE, FALSE)) {
    m1 <- fit_lme(f1, d, het); m0 <- fit_lme(f0, d, het)
    if (!is.null(m1) && !is.null(m0))
      return(list(m1 = m1, m0 = m0, var_model = if (het) "unequal variances (varIdent)" else "equal variances (fallback)"))
  }
  list(m1 = NULL, m0 = NULL, var_model = "failed")
}
lrt_p <- function(m1, m0) {
  if (is.null(m1) || is.null(m0)) return(NA_real_)
  stat <- as.numeric(2 * (logLik(m1) - logLik(m0)))
  df   <- attr(logLik(m1), "df") - attr(logLik(m0), "df")
  pchisq(max(stat, 0), df, lower.tail = FALSE)
}

sd_by_set <- sc_tbl %>% group_by(Set) %>% summarise(sd_rms = sqrt(mean(sd_pooled^2)), .groups = "drop")
pooled_sc <- map_dfr(names(sets), function(s) {
  d <- md_t; d$score <- colMeans(zmat[sets[[s]], d$SampleID, drop = FALSE])
  d$grp <- interaction(d[var_group_by], drop = TRUE)
  st <- fit_pair(score ~ Tissue + Metabolic_State, score ~ Tissue, d)
  it <- fit_pair(score ~ Metabolic_State * Tissue, score ~ Metabolic_State + Tissue, d)
  est <- if (!is.null(st$m1)) unname(nlme::fixef(st$m1)["Metabolic_StateD"]) else NA_real_
  tis_shift <- d %>% group_by(Tissue, Metabolic_State) %>% summarise(m = mean(score), .groups = "drop") %>%
    pivot_wider(names_from = Metabolic_State, values_from = m) %>% mutate(s = D - N)
  tibble(Set = s, n_genes = length(sets[[s]]), variance_model = st$var_model,
         pooled_shift_Z = est,
         g_pooled = est / sd_by_set$sd_rms[match(s, sd_by_set$Set)],
         p_pooled_sc = lrt_p(st$m1, st$m0),
         n_tissues_down = sum(tis_shift$s < 0), n_tissues = nrow(tis_shift),
         p_interaction_EXPLORATORY = lrt_p(it$m1, it$m0))
}) %>% mutate(q_pooled_sc = p.adjust(p_pooled_sc, "BH"),
              q_interaction_EXPLORATORY = p.adjust(p_interaction_EXPLORATORY, "BH"))

# ---- 6. SUMMARIES + TABLES -----------------------------------------------------------
set_tbl <- cam_tbl %>% left_join(sc_tbl, by = c("Tissue", "Set")) %>% arrange(q_camera)
cat("\n=== COMPETITIVE, per tissue (camera; BH over", nrow(cam_tbl), "set x tissue tests) ===\n")
print(as_tibble(set_tbl) %>%
        select(any_of(c("Tissue", "Set", "NGenes", "Correlation", "Direction", "p_camera", "q_camera",
                        "p_mixed", "q_mixed", "shift_set_Z", "g", "q_selfcontained"))), n = 20, width = Inf)
cat("\n=== COMPETITIVE, pooled across tissues (camera on whitened data; BH over", nrow(cam_pooled), "sets) ===\n")
print(as_tibble(cam_pooled) %>% arrange(q_camera), n = Inf, width = Inf)
cat("\n=== SELF-CONTAINED, pooled across tissues (lme; BH over", nrow(pooled_sc), "sets) ===\n")
print(as_tibble(pooled_sc) %>% arrange(q_pooled_sc), n = Inf, width = Inf)

ph <- gene_tbl %>% filter(in_posthoc_list)
cat("\n=== GENE LEVEL ===\n")
cat("Tests in the family:", nrow(gene_tbl), "| nominal p < 0.05:", sum(gene_tbl$p < 0.05),
    "| expected by chance:", round(0.05 * nrow(gene_tbl)), "| BH q <", fdr_cut, ":", sum(gene_tbl$q < fdr_cut), "\n")
cat(flag_label, ":", nrow(ph), "tests | nominal p < 0.05:", sum(ph$p < 0.05),
    "(chance ~", round(0.05 * nrow(ph), 1), ") | BH q <", fdr_cut, ":", sum(ph$q < fdr_cut), "\n")

# ---- 6b. BACKGROUND CHECK: is the rest of the transcriptome shifted? -------------------
# (mean limma t of non-set genes ~ 0 and ~50% up => no transcriptome-wide shift)
bg_check <- map_dfr(tested_tissues, function(tis) {
  md <- metadata %>% filter(Tissue == tis)
  Y  <- Normalizedcounts[, md$SampleID]; Y <- Y[apply(Y, 1, var) > 0, , drop = FALSE]
  tt <- eBayes(lmFit(Y, model.matrix(~ Metabolic_State, data = md)))$t[, 2]
  in_set <- names(tt) %in% set_genes
  tibble(Tissue = tis, n_background = sum(!in_set),
         mean_t_background = mean(tt[!in_set]), frac_background_up = mean(tt[!in_set] > 0),
         n_set_genes = sum(in_set), mean_t_set_genes = mean(tt[in_set]), frac_set_up = mean(tt[in_set] > 0))
})
cat("\n=== BACKGROUND CHECK (background should be ~0 / ~50% up if there is no global shift) ===\n")
print(bg_check, width = Inf)

out <- function(suffix) file.path(main_fig_dir, paste0(FILE_PREFIX, suffix))
write.csv(gene_tbl,    out("_gene_level_FDR.csv"), row.names = FALSE)
write.csv(gene_pooled_tbl, out("_gene_level_pooled_FDR.csv"), row.names = FALSE)
write.csv(set_tbl,     out("_set_level_per_tissue_FDR.csv"), row.names = FALSE)
write.csv(cam_pooled,  out("_set_level_pooled_COMPETITIVE_camera.csv"), row.names = FALSE)
write.csv(pooled_sc,   out("_set_level_pooled_SELFCONTAINED_lme.csv"), row.names = FALSE)
write.csv(pc$weights,  out("_pooled_model_variance_weights.csv"), row.names = FALSE)
write.csv(bg_check,    out("_background_check.csv"), row.names = FALSE)

# ---- 7. MASTER AUDIT CSV --------------------------------------------------------------
# One row per set x analysis (each tested tissue, each untested tissue, and the pooled model):
# which genes were defined / found / used by camera / dropped (and why), plus the camera results.
audit_tbl <- bind_rows(
  unname(map(res_list, "audit")),
  map(untested_tissues, ~ audit_static(.x, character(0), FALSE)),
  list(audit_static(pooled_name, pc$yrows, TRUE)))
cam_all <- bind_rows(
  cam_tbl    %>% mutate(Analysis = as.character(Tissue), FDR_family = "per-tissue (all set x tissue tests)"),
  cam_pooled %>% mutate(Analysis = pooled_name,          FDR_family = "pooled (all sets)")) %>%
  select(any_of(c("Analysis", "Set", "NGenes", "Correlation", "Direction", "p_camera", "q_camera",
                  "p_mixed", "q_mixed", "FDR_family"))) %>%
  rename(camera_NGenes = NGenes, camera_inter_gene_cor = Correlation, camera_direction = Direction)
master_audit <- audit_tbl %>%
  left_join(cam_all, by = c("Analysis", "Set")) %>%
  mutate(Analysis = factor(Analysis, levels = c(tissue_order, pooled_name))) %>%
  arrange(Analysis, Set)
write.csv(master_audit, out("_MASTER_AUDIT.csv"), row.names = FALSE)

# gene membership: every gene x set, with the tissues in which that gene went into a camera run
yrows <- map(res_list, "yrows")
cam_ran <- audit_tbl %>% filter(camera_run, Analysis %in% tested_tissues) %>% select(Set, Analysis)
membership <- imap_dfr(defined_list, ~ tibble(Set = .y, Gene = .x)) %>%
  mutate(in_data = Gene %in% rownames(Normalizedcounts),
         set_was_tested = Set %in% names(sets),
         flagged_Table1 = Gene %in% posthoc_genes)
membership$tissues_in_camera <- map2_chr(membership$Set, membership$Gene, function(s, g) {
  tt <- cam_ran$Analysis[cam_ran$Set == s]
  paste(tt[map_lgl(tt, ~ g %in% yrows[[.x]])], collapse = "; ")
})
write.csv(membership, out("_gene_membership_audit.csv"), row.names = FALSE)

# ---- 8. FIGURES: SET x TISSUE (two honestly-labelled heatmaps) ------------------------
to_mat <- function(tbl, row_col, col_col, val_col, rows, cols, fill = NA_real_) {
  m <- matrix(fill, length(rows), length(cols), dimnames = list(rows, cols))
  m[cbind(match(tbl[[row_col]], rows), match(tbl[[col_col]], cols))] <- tbl[[val_col]]
  m
}
star_fun <- function(q) ifelse(is.na(q), "", ifelse(q < 0.01, "**", ifelse(q < fdr_cut, "*", "")))

col_names <- c(tissue_order, "All tissues")
g_mat     <- to_mat(sc_tbl, "Set", "Tissue", "g",           names(sets), col_names)
shift_mat <- to_mat(sc_tbl, "Set", "Tissue", "shift_set_Z", names(sets), col_names)
g_mat[pooled_sc$Set, "All tissues"]     <- pooled_sc$g_pooled
shift_mat[pooled_sc$Set, "All tissues"] <- pooled_sc$pooled_shift_Z

draw_set_heatmap <- function(star_mat, file, title_txt, pooled_label) {
  lab <- matrix("", nrow(shift_mat), ncol(shift_mat), dimnames = dimnames(shift_mat))
  ok  <- !is.na(shift_mat)
  lab[ok] <- paste0(sprintf("%+.2f", shift_mat[ok]), "\n", star_mat[ok])
  gm  <- pmax(pmin(g_mat, set_g_lim), -set_g_lim)
  tiff(file, width = 14, height = 7.5, units = "in", res = 300, compression = "lzw")
  pheatmap(gm, display_numbers = lab, fontsize_number = 11, number_color = "black",
           cluster_rows = FALSE, cluster_cols = FALSE, gaps_col = length(tissue_order),
           labels_col = c(tissue_order, pooled_label),
           color = my_colors, breaks = seq(-set_g_lim, set_g_lim, length.out = 101), na_col = "grey85",
           border_color = "white", fontsize = 11, fontsize_row = 12, fontsize_col = 12, main = title_txt)
  dev.off()
}

# 8a. COMPETITIVE (camera only; mroast is self-contained, so its '+' lives in figure 8b)
star_comp <- to_mat(mutate(cam_tbl, star = star_fun(q_camera)),
                    "Set", "Tissue", "star", names(sets), col_names, fill = "")
star_comp[cam_pooled$Set, "All tissues"] <- star_fun(cam_pooled$q_camera)
draw_set_heatmap(
  star_comp, out("_Set_level_COMPETITIVE_camera_Heatmap.tiff"),
  paste0(FILE_PREFIX, " - COMPETITIVE (camera): shifted MORE than other genes?\n",
         "Colour = Hedges g (D - N).  Text = mean set shift D - N (Z), then symbols.\n",
         "* camera BH q < ", fdr_cut, " (** < 0.01);  grey = not tested\n",
         "Last column: pooled camera (bird random effect, unequal variances)"),
  "All tissues\n(pooled camera)")

# 8b. SELF-CONTAINED (Welch / lme stars; '+' = any-direction mroast 'mixed', also self-contained)
sc_plot <- sc_tbl %>% left_join(cam_tbl %>% select(Tissue, Set, q_mixed), by = c("Tissue", "Set"))
star_sc <- to_mat(mutate(sc_plot, star = paste0(star_fun(q_selfcontained),
                                                ifelse(!is.na(q_mixed) & q_mixed < fdr_cut, "+", ""))),
                  "Set", "Tissue", "star", names(sets), col_names, fill = "")
qm_pooled <- cam_pooled$q_mixed[match(pooled_sc$Set, cam_pooled$Set)]
star_sc[pooled_sc$Set, "All tissues"] <- paste0(star_fun(pooled_sc$q_pooled_sc),
                                                ifelse(!is.na(qm_pooled) & qm_pooled < fdr_cut, "+", ""))
draw_set_heatmap(
  star_sc, out("_Set_level_SELFCONTAINED_Heatmap.tiff"),
  paste0(FILE_PREFIX, " - SELF-CONTAINED (set score): does the set shift AT ALL?\n",
         "Colour = Hedges g (D - N).  Text = mean set shift D - N (Z), then symbols.\n",
         "* Welch t BH q < ", fdr_cut, " (** < 0.01);  + any-direction (mroast mixed) q < ", fdr_cut, ";  grey = not tested\n",
         "Last column: lme, bird random intercept + unequal variances (LRT)"),
  "All tissues\n(pooled lme)")

# ---- 9. FIGURE: ALL GENES IN THE SETS ---------------------------------------------------
gene_levels <- setdiff(unique(gene_tbl$Set), kegg_only_label)   # only groups that appear in the gene tables
clean_lab <- function(x) trimws(sub("\\s*\\((non-)?KEGG[^)]*\\)\\s*$", "", x))   # legend text without (KEGG)/(non-KEGG)
gene_order  <- gene_tbl %>% group_by(Gene, Set) %>% summarise(mean_shift = mean(shift_Z), .groups = "drop") %>%
  mutate(Set = factor(Set, levels = gene_levels)) %>% arrange(Set, mean_shift)

gcols <- c(tissue_order, "All tissues")
gt   <- to_mat(gene_tbl, "Gene", "Tissue", "t",       gene_order$Gene, gcols)
gsh  <- to_mat(gene_tbl, "Gene", "Tissue", "shift_Z", gene_order$Gene, gcols)
gst  <- to_mat(mutate(gene_tbl, star = ifelse(q < fdr_cut, "*", "")), "Gene", "Tissue", "star",
               gene_order$Gene, gcols, fill = "")
pi_  <- match(gene_order$Gene, gene_pooled_tbl$Gene)
gt[, "All tissues"]  <- gene_pooled_tbl$t_pooled[pi_]
gsh[, "All tissues"] <- gene_pooled_tbl$mean_shift_Z[pi_]
gst[, "All tissues"] <- star_fun(gene_pooled_tbl$q_pooled[pi_])
glab <- matrix("", nrow(gsh), ncol(gsh), dimnames = dimnames(gsh))
ok   <- !is.na(gsh)
glab[ok] <- paste0(sprintf("%+.1f", gsh[ok]), gst[ok])
gt_c <- pmax(pmin(gt, gene_t_lim), -gene_t_lim)

ann_row <- data.frame(Set = factor(clean_lab(as.character(gene_order$Set)), levels = clean_lab(gene_levels)),
                      row.names = gene_order$Gene)
ann_col <- list(Set = setNames(scales::hue_pal()(length(gene_levels)), clean_lab(gene_levels)))
gene_gaps <- cumsum(table(gene_order$Set)); gene_gaps <- unique(gene_gaps[gene_gaps > 0 & gene_gaps < nrow(gt)])
tiff(out("_All_Set_Genes_FDR_Heatmap.tiff"),
     width = 11, height = max(8.5, 0.22 * nrow(gt) + 2), units = "in", res = 300, compression = "lzw")
pheatmap(gt_c, display_numbers = glab, fontsize_number = 6.5, number_color = "black",
         cluster_rows = FALSE, cluster_cols = FALSE, gaps_col = length(tissue_order),
         annotation_row = ann_row, annotation_colors = ann_col, gaps_row = gene_gaps,
         labels_col = c("Proximal\nGut", "Medial\nGut", "Distal\nGut", "Lungs", "Heart", "Liver", "Pectoral", "All tissues\n(pooled)"),
         angle_col = "0",
         color = my_colors, breaks = seq(-gene_t_lim, gene_t_lim, length.out = 101), na_col = "grey85",
         border_color = NA, fontsize = 9, fontsize_row = 7, fontsize_col = 10,
         main = paste0(FILE_PREFIX, ": every gene in the pre-specified sets\n",
                       "Colour = limma moderated t (D - N)\n",
                       "Text = mean shift D - N (Z); * BH q < ", fdr_cut, ", ** < 0.01\n",
                       "Tissue columns: BH over ", nrow(gene_tbl), " gene x tissue tests\n",
                       "Last column: pooled limma (bird effect, unequal variances);\n",
                       "text = mean of tissue shifts; BH over ", nrow(gene_pooled_tbl), " genes"))
dev.off()

cat("\nDone. Tables, audit CSVs and figures are in:", main_fig_dir, "\n")

# ==========================================================================
# 10. OPTIONAL: PATHVIEW MAPS (DESCRIPTIVE ONLY: node colour = mean Z shift D - N)
# ==========================================================================
if (run_pathview) {
  suppressPackageStartupMessages(library(pathview))
  pv_dir <- file.path(main_fig_dir, "pathview")
  dir.create(pv_dir, showWarnings = FALSE, recursive = TRUE)
  
  get_symbols_from_ids <- function(id_string) {
    if (is.na(id_string) || id_string == "") return(NA)
    ids  <- unlist(strsplit(as.character(id_string), ","))
    syms <- tryCatch(mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
                     error = function(e) ids)
    syms[is.na(syms)] <- ids[is.na(syms)]
    paste(syms, collapse = ", ")
  }
  
  original_wd <- getwd()
  setwd(pv_dir)
  for (tis in tested_tissues) {
    md <- metadata %>% filter(Tissue == tis)
    idN <- md$SampleID[md$Metabolic_State == "N"]; idD <- md$SampleID[md$Metabolic_State == "D"]
    diff_vector <- rowMeans(zmat[, idD, drop = FALSE]) - rowMeans(zmat[, idN, drop = FALSE])
    
    entrez_map <- AnnotationDbi::select(org.Hs.eg.db, keys = names(diff_vector), columns = "ENTREZID", keytype = "SYMBOL")
    pv_df <- merge(data.frame(SYMBOL = names(diff_vector), value = diff_vector), entrez_map, by = "SYMBOL")
    pv_df <- pv_df[!is.na(pv_df$ENTREZID), ]
    pv_df <- aggregate(value ~ ENTREZID, data = pv_df, mean)
    gene_data_vector <- setNames(pv_df$value, pv_df$ENTREZID)
    
    pv_out <- tryCatch(
      pathview(gene.data = gene_data_vector, pathway.id = KEGG_ID, species = "hsa",
               out.suffix = paste0(FILE_PREFIX, "_", gsub(" ", "", tis)), limit = list(gene = 2, cpd = 1),
               low = "goldenrod2", mid = "gray", high = "purple", kegg.native = TRUE),
      error = function(e) { message("pathview failed for ", tis, ": ", conditionMessage(e)); NULL })
    
    if (!is.null(pv_out) && !is.null(pv_out$plot.data.gene)) {
      nd <- pv_out$plot.data.gene
      audit <- data.frame(Node_Label = nd$labels, My_Mapped_Entrez = nd$all.mapped, stringsAsFactors = FALSE)
      audit$Tissue <- tis
      audit <- audit[!is.na(audit$My_Mapped_Entrez) & audit$My_Mapped_Entrez != "", ]
      audit$My_Genes_In_Bin <- sapply(audit$My_Mapped_Entrez, get_symbols_from_ids)
      write.csv(audit, file.path(main_fig_dir, paste0(FILE_PREFIX, "_Node_Audit_", gsub(" ", "", tis), ".csv")), row.names = FALSE)
    }
  }
  setwd(original_wd)
  cat("Pathview maps written to:", pv_dir, "\n")
}