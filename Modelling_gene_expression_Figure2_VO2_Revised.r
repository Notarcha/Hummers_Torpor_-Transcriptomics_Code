# ==========================================================================
# Hummingbird torpor: temperature-responsive genes, linear vs threshold-like
# FULL SCRIPT = revised classification (sections 0-7) + within-state
# permutation null (section 8). Run top to bottom.
#
# Changes vs. the ORIGINAL analysis:
#  1. Z-scores computed PER TISSUE (original grouped by gene only).
#  2. Bird random intercept (original: tissue intercept only).
#  3. Screen = FDR-corrected LRT of a NEUTRAL model (temp + state) vs null,
#     replacing the uncorrected pooled ANOVA that favoured the state model.
#  4. Fair model set: linear (1 df) vs early step (1 df) vs late step (1 df)
#     vs 3-level categorical (2 df).
#  5. Ambiguous genes reported as ambiguous (not folded into Early/Late).
#  6. Direct step-vs-continuous test: within-state temperature slope.
#  7. Tukey tests on pre-selected genes removed (circular + pseudoreplicated).
#  8. Permutation null (shuffle temperature among birds WITHIN state, rerun
#     the whole pipeline incl. screening) to calibrate #6 and the class split.
# ==========================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(lme4)
  library(ComplexHeatmap); library(circlize); library(here)
})
set.seed(1)

# ---- 0. SETTINGS ---------------------------------------------------------
predictor       <- "log_VO2"   # or "Surf_Temp" for the main analysis
fdr_screen      <- 0.05          # FDR for "is the gene responsive at all"
aic_tol         <- 2             # models within this of the best are "supported"
n_cores         <- max(1, parallel::detectCores() - 1)  # set to 1 on Windows

run_permutation <- TRUE          # section 8; set FALSE to skip
n_perm          <- 100           # try 5 first to check timing; 100-200 for the paper
cand_p          <- 0.05          # permutations re-screen only genes with observed
# screen p < cand_p (speed-up; set 1 for exact, ~10x slower)

output_dir <- "/Users/notarcha/Desktop/PNAS_Code/NewOutputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# predictor prefix so a log_VO2 run does not overwrite the Surf_Temp run
out <- function(f) file.path(output_dir, paste0(predictor, "_", f))

# ---- 1. LOAD + ALIGN -----------------------------------------------------
seq_metadata <- read.csv(here("Results", "RNASeq_metadata.csv"), stringsAsFactors = FALSE)
rlog_counts  <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"),
                         row.names = 1, check.names = FALSE)
metabolic    <- read.csv(here("WGCNA_data", "metadata_surfTemp_mlO2_min_euthTime.csv"),
                         stringsAsFactors = FALSE)

# NOTE: the original rowSums(rlog) > 10 filter does nothing meaningful on
# log-scale values. Apply your expression filter to RAW counts upstream
# (e.g. >=10 counts in >=3 samples of at least one tissue) and re-export.

sample_info <- seq_metadata %>%
  left_join(metabolic %>% select(BirdID, Surf_Temp, mean_VO2), by = "BirdID") %>%
  rename(Mean_VO2 = mean_VO2)

common      <- intersect(colnames(rlog_counts), sample_info$X)
sample_info <- sample_info[match(common, sample_info$X), ] %>%
  mutate(
    State  = factor(Metabolic_State, levels = c("N", "T", "D"),
                    labels = c("Normothermy", "Transition", "Deep Torpor")),
    # full tissue names (same labels as the other figures); level order = old alphabetical order, so colours are unchanged
    Tissue = forcats::fct_recode(factor(Tissue), `Proximal Gut` = "Gut1", `Medial Gut` = "Gut2",
                                 `Distal Gut` = "Gut3", Pectoral = "Pect"),
    BirdID = factor(BirdID),
    X_raw  = if (predictor == "Surf_Temp") Surf_Temp else log10(Mean_VO2)
  ) %>%
  filter(!is.na(X_raw), !is.na(State)) %>%
  arrange(X_raw) %>%
  mutate(
    X_c   = as.numeric(scale(X_raw)),
    Early = as.numeric(State != "Normothermy"),  # step between N and T
    Late  = as.numeric(State == "Deep Torpor")   # step between T and D
  )

# bird-level sanity checks (each bird: one state, one predictor value)
stopifnot(
  all(tapply(as.character(sample_info$State), sample_info$BirdID, function(x) length(unique(x))) == 1),
  all(tapply(sample_info$X_raw, sample_info$BirdID, function(x) length(unique(x))) == 1)
)

# which samples were lost, and why (reconcile with the Methods: "3 outliers removed")
message("Samples in count matrix: ", ncol(rlog_counts), " | retained: ", nrow(sample_info))
dropped <- setdiff(colnames(rlog_counts), sample_info$X)
if (length(dropped) > 0) message("Not retained: ", paste(dropped, collapse = ", "))

counts <- as.matrix(rlog_counts[, sample_info$X])

# ---- 2. PER-TISSUE Z-SCORES (matches Methods; original was pooled) -------
zmat <- counts
for (tt in levels(sample_info$Tissue)) {
  idx <- which(sample_info$Tissue == tt)
  x   <- counts[, idx, drop = FALSE]
  mu  <- rowMeans(x)
  zmat[, idx] <- (x - mu) / sqrt(rowSums((x - mu)^2) / (length(idx) - 1))
}
zmat <- zmat[apply(zmat, 1, function(r) all(is.finite(r))), ]
message("Genes retained: ", nrow(zmat), " | samples: ", ncol(zmat))

# sample sizes per state x tissue (report these in the paper)
print(with(sample_info, table(Tissue, State)))

# ---- 3. HELPERS ----------------------------------------------------------
fit_ml <- function(formula, dat) {
  tryCatch(
    suppressMessages(suppressWarnings(
      lmer(formula, data = dat, REML = FALSE,
           control = lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE))
    )),
    error = function(e) NULL
  )
}
ll  <- function(m) as.numeric(logLik(m))
dfm <- function(m) attr(logLik(m), "df")

# parallel/serial lapply; extra args (e.g. si = ...) are passed through to f
run_par <- function(x, f, ...) {
  r <- if (n_cores > 1) parallel::mclapply(x, f, ..., mc.cores = n_cores) else lapply(x, f, ...)
  Filter(function(z) !inherits(z, "try-error"), r)
}

# ---- 4. STAGE 1: IS THE GENE RESPONSIVE AT ALL? --------------------------
# Neutral screen: model containing BOTH continuous temp and discrete state
# (nests the linear and step models) vs intercept-only. LRT, BH-corrected.
# `si` = sample table (so permutations can pass a shuffled copy).
stage1 <- function(g, si) {
  dat <- si; dat$y <- zmat[g, ]
  m0 <- fit_ml(y ~ 1 + (1 | BirdID), dat)
  mu <- fit_ml(y ~ X_c + State + (1 | BirdID), dat)
  if (is.null(m0) || is.null(mu)) return(tibble(Gene = g, p_union = NA_real_))
  lrt <- max(2 * (ll(mu) - ll(m0)), 0)
  tibble(Gene = g, p_union = pchisq(lrt, dfm(mu) - dfm(m0), lower.tail = FALSE))
}

res1 <- bind_rows(run_par(rownames(zmat), stage1, si = sample_info)) %>%
  mutate(q_union = p.adjust(p_union, "BH"))
sig_genes <- res1 %>% filter(q_union < fdr_screen) %>% pull(Gene)
message("Responsive genes (FDR < ", fdr_screen, "): ", length(sig_genes),
        " of ", nrow(res1), " tested")

# ---- 5. STAGE 2: FAIR MODEL COMPARISON ON RESPONSIVE GENES ---------------
stage2 <- function(g, si) {
  dat <- si; dat$y <- zmat[g, ]
  fits <- lapply(list(
    lin   = y ~ X_c   + (1 | BirdID),                  # 1 df, continuous
    early = y ~ Early + (1 | BirdID),                  # 1 df, N | T+D
    late  = y ~ Late  + (1 | BirdID),                  # 1 df, N+T | D
    cat   = y ~ State + (1 | BirdID),                  # 2 df, any pattern
    union = y ~ X_c + State + (1 | BirdID)             # 3 df
  ), fit_ml, dat = dat)
  if (any(vapply(fits, is.null, logical(1)))) return(NULL)
  aic <- vapply(fits, AIC, numeric(1))
  # Does temperature explain anything beyond state? (within-state slope)
  p_within <- pchisq(max(2 * (ll(fits$union) - ll(fits$cat)), 0), 1, lower.tail = FALSE)
  tibble(Gene = g, AIC_lin = aic[["lin"]], AIC_early = aic[["early"]],
         AIC_late = aic[["late"]], AIC_cat = aic[["cat"]], p_within = p_within)
}

classify_gene <- function(aic, tol = 2) {
  supported <- names(aic)[aic - min(aic) <= tol]
  # categorical nests the steps; its being "supported" alongside a step
  # model adds no information, so prefer the simpler step model
  if ("cat" %in% supported && any(c("early", "late") %in% supported))
    supported <- setdiff(supported, "cat")
  has_lin <- "lin" %in% supported
  nonlin  <- setdiff(supported, "lin")
  broad <- if (has_lin && length(nonlin) == 0) "Linear"
  else if (!has_lin) "Non-linear"
  else "Indistinguishable"
  fine <- broad
  if (broad == "Non-linear") {
    fine <- if (identical(nonlin, "early")) "Early step"
    else if (identical(nonlin, "late")) "Late step"
    else if (identical(nonlin, "cat")) "Graded / non-monotonic"
    else "Early/late step unresolved"
  }
  list(broad = broad, fine = fine)
}

res2 <- bind_rows(run_par(sig_genes, stage2, si = sample_info))

# descriptive state means (per-tissue Z) for direction and magnitudes
state_means <- sapply(levels(sample_info$State), function(s)
  rowMeans(zmat[, sample_info$State == s, drop = FALSE]))
desc <- as_tibble(state_means, rownames = "Gene") %>%
  mutate(Early_Change = abs(Transition - Normothermy),
         Late_Change  = abs(`Deep Torpor` - Transition),
         State_Var    = apply(cbind(Normothermy, Transition, `Deep Torpor`), 1, var),
         Direction    = ifelse(`Deep Torpor` > Normothermy, "Up in torpor", "Down in torpor"))

results <- res2 %>%
  rowwise() %>%
  mutate(cl = list(classify_gene(c(lin = AIC_lin, early = AIC_early,
                                   late = AIC_late, cat = AIC_cat), aic_tol))) %>%
  ungroup() %>%
  mutate(Class_Broad = map_chr(cl, "broad"), Class_Fine = map_chr(cl, "fine"),
         Delta_AIC_fair = AIC_lin - pmin(AIC_early, AIC_late)) %>%   # equal-df comparison
  select(-cl) %>%
  left_join(desc, by = "Gene") %>%
  left_join(res1, by = "Gene")

# ---- 6. SUMMARY + OUTPUT -------------------------------------------------
cat("\n--- Broad classes (N genes; % of responsive genes) ---\n")
print(results %>% count(Class_Broad) %>% mutate(pct = round(100 * n / sum(n), 1)))
cat("\n--- Fine classes by direction ---\n")
print(table(results$Class_Fine, results$Direction))
cat("\n--- Within-state temperature slope: fraction p < 0.05 (null expectation ~5%) ---\n")
print(mean(results$p_within < 0.05))

write.csv(res1,    out("screen_all_genes.csv"), row.names = FALSE)
write.csv(results, out("classified_responsive_genes.csv"), row.names = FALSE)
for (k in unique(results$Class_Fine)) {
  write.csv(filter(results, Class_Fine == k),
            out(paste0("genes_", gsub("[^A-Za-z]+", "_", k), ".csv")),
            row.names = FALSE)
}

# ---- 7. FIGURES ----------------------------------------------------------
state_colors <- c("Normothermy" = "#7570B3", "Transition" = "#1B9E77", "Deep Torpor" = "#D95F02")

# 7a. Equal-df comparison: linear vs best 1-df step model
p_aic <- ggplot(results, aes(Delta_AIC_fair)) +
  annotate("rect", xmin = -aic_tol, xmax = aic_tol, ymin = 0, ymax = Inf, fill = "grey85", alpha = 0.4) +
  geom_density(fill = "#5ab4ac", alpha = 0.55) +
  geom_vline(xintercept = c(-aic_tol, aic_tol), linetype = "dashed", colour = "firebrick") +
  coord_cartesian(xlim = c(-15, 15)) +
  labs(x = "AIC(linear) - AIC(best step model)   [>0 favours step]", y = "Density",
       title = "Linear vs best step model, equal parameters",
       subtitle = paste0("Responsive genes, N = ", nrow(results))) +
  theme_minimal(base_size = 12)
ggsave(out("delta_aic_equal_df.png"), p_aic, width = 8, height = 5, dpi = 300)
print(table(cut(results$Delta_AIC_fair, c(-Inf, -aic_tol, aic_tol, Inf),
                labels = c("linear better", "tie", "step better"))))

# 7b. Direct step-vs-continuous diagnostic (observed; calibrated in section 8)
p_within_obs <- ggplot(results, aes(p_within)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey60", colour = "white") +
  labs(x = "p: temperature effect beyond state (within-state slope)", y = "Genes",
       title = "Is there continuous temperature dependence within states?",
       subtitle = "Uniform = pure step; excess near 0 = continuous component") +
  theme_minimal(base_size = 12)
ggsave(out("within_state_slope_pvalues.png"), p_within_obs, width = 7, height = 4.5, dpi = 300)

# 7c. Class trajectories at bird level (descriptive; no tests on selected genes)
cls <- results %>% filter(Class_Fine %in% c("Linear", "Early step", "Late step")) %>%
  select(Gene, Class_Fine, Direction)
traj <- as_tibble(zmat[cls$Gene, , drop = FALSE], rownames = "Gene") %>%
  pivot_longer(-Gene, names_to = "X", values_to = "Z") %>%
  inner_join(cls, by = "Gene") %>%
  inner_join(select(sample_info, X, BirdID, State, X_raw), by = "X") %>%
  group_by(Class_Fine, Direction, BirdID, State, X_raw) %>%
  summarise(Z = mean(Z), .groups = "drop")   # one point per bird

# state-level mean +/- SE computed explicitly (stat_summary would group by the
# continuous x, giving n = 1 per group and NA error bars)
traj_sum <- traj %>%
  group_by(Class_Fine, Direction, State) %>%
  summarise(x = mean(X_raw), m = mean(Z), se = sd(Z) / sqrt(n()), .groups = "drop")

p_traj <- ggplot(traj, aes(X_raw, Z, colour = State)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_pointrange(data = traj_sum, aes(x = x, y = m, ymin = m - se, ymax = m + se),
                  colour = "black", shape = 18, size = 0.5, inherit.aes = FALSE) +
  facet_grid(Direction ~ Class_Fine, scales = "free_y") +
  scale_colour_manual(values = state_colors) +
  labs(x = predictor, y = "Mean per-tissue Z (one point per bird)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave(out("class_trajectories_bird_level.png"), p_traj, width = 11, height = 7, dpi = 300)

# 7d. Heatmap (descriptive; top 100 per class by state variance).
#     Columns ordered warm -> cold (normothermy -> transition -> torpor);
#     tissue track added so per-bird stripes can be read as such.
top <- results %>% filter(Class_Fine %in% c("Linear", "Early step", "Late step")) %>%
  group_by(Class_Fine) %>% slice_max(State_Var, n = 100, with_ties = FALSE) %>% ungroup()
si_h <- sample_info[order(sample_info$X_raw, decreasing = TRUE), ]
mat  <- zmat[top$Gene, si_h$X]
tissue_cols <- setNames(scales::hue_pal()(nlevels(si_h$Tissue)), levels(si_h$Tissue))
ht_opt$message <- FALSE
top_annot <- HeatmapAnnotation(
  `Metabolic Rate`  = si_h$Mean_VO2,
  `Surf Temp`       = si_h$Surf_Temp,
  `Metabolic State` = si_h$State,
  Tissue            = si_h$Tissue,
  col = list(`Metabolic State` = state_colors, Tissue = tissue_cols))
ht <- Heatmap(mat, name = "Per-tissue Z",
              col = colorRamp2(c(-2, 0, 2), c("goldenrod3", "white", "purple")),
              top_annotation = top_annot,
              row_split = factor(top$Class_Fine, levels = c("Linear", "Early step", "Late step")),
              cluster_row_slices = FALSE, cluster_rows = TRUE, cluster_columns = FALSE,
              show_row_names = TRUE, row_names_side = "right",
              row_names_gp = gpar(fontsize = min(9, max(2.5, 700 / nrow(mat)))),  # scales with number of genes
              show_column_names = FALSE,
              use_raster = FALSE)   # raster needs cairo/XQuartz, which is missing on this Mac
png(out("heatmap_top100_per_class.png"), width = 12, height = 14, units = "in", res = 300)
draw(ht, merge_legend = TRUE)
dev.off()

# ==========================================================================
# 8. PERMUTATION NULL: shuffle temperature among BIRDS WITHIN each state
#    Keeps state structure and each bird's tissues together; breaks any real
#    within-state expression-temperature relationship. Each permutation reruns
#    screen -> model comparison -> classification, so screening-induced
#    inflation is reproduced in the null.
# ==========================================================================
if (run_permutation) {
  m_tested <- sum(!is.na(res1$p_union))
  cand <- res1 %>% filter(!is.na(p_union), p_union < cand_p) %>% pull(Gene)
  # BH uses n = all genes tested (others treated as p = 1), matching the
  # observed analysis; genes with observed p >= cand_p essentially never pass
  # once temperature is randomised.
  message("Permutation: ", length(cand), " candidate genes re-screened per permutation (of ", m_tested, ")")
  
  bins <- seq(0, 1, 0.05)
  obs <- tibble(
    n_sig = nrow(results),
    n_within = sum(results$p_within < 0.05),
    frac_within = mean(results$p_within < 0.05),
    frac_linear = mean(results$Class_Broad == "Linear"),
    frac_nonlin = mean(results$Class_Broad == "Non-linear"),
    frac_indist = mean(results$Class_Broad == "Indistinguishable")
  )
  obs_hist <- hist(results$p_within, breaks = bins, plot = FALSE)$counts / nrow(results)
  
  permute_temp <- function(si) {
    birds <- si %>% distinct(BirdID, State, X_raw)
    stopifnot(nrow(birds) == n_distinct(si$BirdID))
    birds <- birds %>% group_by(State) %>%
      mutate(X_perm = X_raw[sample.int(n())]) %>%   # sample(x) misbehaves for length-1 x
      ungroup()
    si %>% select(-X_raw) %>%
      left_join(select(birds, BirdID, X_raw = X_perm), by = "BirdID") %>%
      mutate(X_c = as.numeric(scale(X_raw)))        # row order preserved -> still matches zmat
  }
  
  empty_perm <- tibble(n_sig = 0L, n_within = 0L, frac_within = NA_real_,
                       n_linear = 0L, n_nonlin = 0L, n_indist = 0L,
                       hist = list(rep(NA_real_, length(bins) - 1)))
  
  run_perm <- function(si) {
    r1 <- bind_rows(run_par(cand, stage1, si = si)) %>% filter(!is.na(p_union))
    r1$q <- p.adjust(r1$p_union, "BH", n = m_tested)
    sg <- r1$Gene[r1$q < fdr_screen]
    if (length(sg) == 0) return(empty_perm)
    r2 <- bind_rows(run_par(sg, stage2, si = si))
    if (nrow(r2) == 0) return(empty_perm)
    broad <- vapply(seq_len(nrow(r2)), function(i)
      classify_gene(c(lin = r2$AIC_lin[i], early = r2$AIC_early[i],
                      late = r2$AIC_late[i], cat = r2$AIC_cat[i]), aic_tol)$broad,
      character(1))
    tibble(n_sig = length(sg),
           n_within = sum(r2$p_within < 0.05),
           frac_within = mean(r2$p_within < 0.05),
           n_linear = sum(broad == "Linear"),
           n_nonlin = sum(broad == "Non-linear"),
           n_indist = sum(broad == "Indistinguishable"),
           hist = list(hist(r2$p_within, breaks = bins, plot = FALSE)$counts / nrow(r2)))
  }
  
  perm_res <- vector("list", n_perm)
  t0 <- Sys.time()
  for (b in seq_len(n_perm)) {
    perm_res[[b]] <- run_perm(permute_temp(sample_info))
    if (b %% 5 == 0 || b == n_perm || b == 1) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message(sprintf("perm %d/%d | %.1f min elapsed | ~%.1f min remaining",
                      b, n_perm, el, el / b * (n_perm - b)))
      saveRDS(perm_res, out("perm_partial.rds"))
    }
  }
  
  perm <- bind_rows(perm_res, .id = "perm") %>%
    mutate(frac_linear = n_linear / n_sig, frac_nonlin = n_nonlin / n_sig,
           frac_indist = n_indist / n_sig)
  write.csv(select(perm, -hist), out("permutation_null_summary.csv"), row.names = FALSE)
  
  # empirical one-sided p-values (+1 correction)
  emp_p <- function(null, observed) {
    null <- null[is.finite(null)]
    (1 + sum(null >= observed)) / (1 + length(null))
  }
  mets <- list(
    "N responsive genes"          = list(obs$n_sig,       perm$n_sig),
    "Frac. within-state p < 0.05" = list(obs$frac_within, perm$frac_within),
    "Frac. Linear"                = list(obs$frac_linear, perm$frac_linear),
    "Frac. Non-linear"            = list(obs$frac_nonlin, perm$frac_nonlin),
    "Frac. Indistinguishable"     = list(obs$frac_indist, perm$frac_indist)
  )
  cmp <- imap_dfr(mets, function(m, nm) tibble(
    metric = nm, observed = m[[1]],
    null_mean = mean(m[[2]], na.rm = TRUE),
    null_p95  = unname(quantile(m[[2]], 0.95, na.rm = TRUE)),
    p_emp     = emp_p(m[[2]], m[[1]])))
  cat("\n--- Observed vs permutation null (", n_perm, " permutations) ---\n", sep = "")
  print(as.data.frame(cmp), digits = 3)
  excess <- obs$frac_within - mean(perm$frac_within, na.rm = TRUE)
  cat("\nContinuous-component estimate = observed frac(p<0.05) - null mean = ",
      round(excess, 3), " (~", round(excess * obs$n_sig), " genes)\n", sep = "")
  write.csv(cmp, out("permutation_null_comparison.csv"), row.names = FALSE)
  
  # A. within-state p-value histogram: observed vs null band
  H <- do.call(rbind, perm$hist)
  H <- H[complete.cases(H), , drop = FALSE]
  band <- tibble(mid = bins[-length(bins)] + 0.025, obs = obs_hist,
                 lo = apply(H, 2, quantile, 0.025), mean = colMeans(H),
                 hi = apply(H, 2, quantile, 0.975))
  pA <- ggplot(band, aes(mid)) +
    geom_col(aes(y = obs), width = 0.047, fill = "grey60") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "firebrick", alpha = 0.25) +
    geom_line(aes(y = mean), colour = "firebrick", linewidth = 0.8) +
    labs(x = "p: temperature effect beyond state (within-state slope)",
         y = "Proportion of responsive genes",
         title = "Within-state temperature effect: observed vs permutation null",
         subtitle = paste0("Grey = observed; red = null mean and 95% band (",
                           nrow(H), " permutations incl. screening)")) +
    theme_minimal(base_size = 12)
  ggsave(out("permutation_within_state_hist.png"), pA, width = 8, height = 5, dpi = 300)
  
  # B. class proportions: observed vs null distribution
  long_null <- perm %>%
    select(perm, Linear = frac_linear, `Non-linear` = frac_nonlin,
           Indistinguishable = frac_indist) %>%
    pivot_longer(-perm, names_to = "Class", values_to = "frac") %>% filter(is.finite(frac))
  long_obs <- tibble(Class = c("Linear", "Non-linear", "Indistinguishable"),
                     frac = c(obs$frac_linear, obs$frac_nonlin, obs$frac_indist))
  pB <- ggplot(long_null, aes(Class, frac)) +
    geom_violin(fill = "grey85", colour = "grey40") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    geom_point(data = long_obs, colour = "firebrick", size = 4, shape = 18) +
    labs(x = NULL, y = "Proportion of responsive genes",
         title = "Class proportions: observed (red) vs permutation null",
         subtitle = "Null = all temperature information at state level only") +
    theme_minimal(base_size = 12)
  ggsave(out("permutation_class_proportions.png"), pB, width = 7, height = 5, dpi = 300)
}

# ---- 9. RECOMMENDED SENSITIVITY RUN --------------------------------------
# Re-run the whole script with predictor <- "log_VO2": Transition birds were
# defined by surface temperature, and VO2 had already collapsed by then
# (Fig 1A). Outputs are prefixed by predictor, so nothing is overwritten.



# ==========================================================================
# Heatmap only: reads the saved classification CSV, so NO model refitting.
# Needs: <predictor>_classified_responsive_genes.csv from the main script.
# ==========================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(ComplexHeatmap); library(circlize); library(here)
})

predictor  <- "log_VO2"   # must match the run whose CSV you want to plot
output_dir <- "/Users/notarcha/Desktop/PNAS_Code/NewOutputs"
out <- function(f) file.path(output_dir, paste0(predictor, "_", f))
top_n_per_class <- 100

# ---- load data (same alignment as the main script; takes seconds) ----------
seq_metadata <- read.csv(here("Results", "RNASeq_metadata.csv"), stringsAsFactors = FALSE)
rlog_counts  <- read.csv(here("Results", "rlog_normalized_counts_STAR.csv"),
                         row.names = 1, check.names = FALSE)
metabolic    <- read.csv(here("WGCNA_data", "metadata_surfTemp_mlO2_min_euthTime.csv"),
                         stringsAsFactors = FALSE)

sample_info <- seq_metadata %>%
  left_join(metabolic %>% select(BirdID, Surf_Temp, mean_VO2), by = "BirdID") %>%
  rename(Mean_VO2 = mean_VO2)
common      <- intersect(colnames(rlog_counts), sample_info$X)
sample_info <- sample_info[match(common, sample_info$X), ] %>%
  mutate(
    State  = factor(Metabolic_State, levels = c("N", "T", "D"),
                    labels = c("Normothermy", "Transition", "Deep Torpor")),
    # full tissue names (same labels as the other figures); level order = old alphabetical order, so colours are unchanged
    Tissue = forcats::fct_recode(factor(Tissue), `Proximal Gut` = "Gut1", `Medial Gut` = "Gut2",
                                 `Distal Gut` = "Gut3", Pectoral = "Pect"),
    X_raw  = if (predictor == "Surf_Temp") Surf_Temp else log10(Mean_VO2)
  ) %>%
  filter(!is.na(X_raw), !is.na(State)) %>%
  arrange(X_raw)

counts <- as.matrix(rlog_counts[, sample_info$X])
zmat <- counts                                   # per-tissue Z-scores
for (tt in levels(sample_info$Tissue)) {
  idx <- which(sample_info$Tissue == tt)
  x   <- counts[, idx, drop = FALSE]
  mu  <- rowMeans(x)
  zmat[, idx] <- (x - mu) / sqrt(rowSums((x - mu)^2) / (length(idx) - 1))
}
zmat <- zmat[apply(zmat, 1, function(r) all(is.finite(r))), ]

# ---- saved results ---------------------------------------------------------
results <- read.csv(out("classified_responsive_genes.csv"), stringsAsFactors = FALSE,
                    check.names = FALSE)

# ---- heatmap ---------------------------------------------------------------
state_colors <- c("Normothermy" = "#7570B3", "Transition" = "#1B9E77", "Deep Torpor" = "#D95F02")
classes <- c("Linear", "Early step", "Late step", "Graded / non-monotonic")

top <- results %>% filter(Class_Fine %in% classes) %>%
  group_by(Class_Fine) %>% slice_max(State_Var, n = top_n_per_class, with_ties = FALSE) %>% ungroup()

si_h <- sample_info[order(sample_info$X_raw, decreasing = TRUE), ]   # warm -> cold
mat  <- zmat[top$Gene, si_h$X]
tissue_cols <- setNames(scales::hue_pal()(nlevels(si_h$Tissue)), levels(si_h$Tissue))

ht_opt$message <- FALSE
top_annot <- HeatmapAnnotation(
  `Metabolic Rate`  = si_h$Mean_VO2,
  `Surf Temp`       = si_h$Surf_Temp,
  `Metabolic State` = si_h$State,
  Tissue            = si_h$Tissue,
  col = list(`Metabolic State` = state_colors, Tissue = tissue_cols))

ht <- Heatmap(mat, name = "Per-tissue Z",
              col = colorRamp2(c(-2, 0, 2), c("goldenrod3", "white", "purple")),
              top_annotation = top_annot,
              row_split = factor(top$Class_Fine, levels = classes),
              cluster_row_slices = FALSE, cluster_rows = TRUE, cluster_columns = FALSE,
              show_row_names = TRUE, row_names_side = "right",
              row_names_gp = gpar(fontsize = min(9, max(2.5, 700 / nrow(mat)))),
              show_column_names = FALSE,
              use_raster = FALSE)

png(out("heatmap_top100_per_class_genenames.png"), width = 12, height = 16, units = "in", res = 300)   # taller: four classes now
draw(ht, merge_legend = TRUE)
dev.off()
message("Saved: ", out("heatmap_top100_per_class_genenames.png"))

# ==========================================================================
# Per-tissue Z-score expression vs VO2 for selected genes, faceted by
# gene, with the best-fitting model (lowest AIC) drawn over the data.
#
# Run in the same session as the main script (needs: zmat, sample_info, out(),
# predictor). Rebuilds the few model columns itself, so it also works after the
# heatmap-only block has overwritten sample_info.
# ==========================================================================
library(tidyverse); library(lme4)

genes <- c("NR1D1", "RSRP1", "SGK1", "LUC7L3", "BHLHE40",
           "RBM12B", "OTX1", "CLK4", "TADA2B", "VAMP1")

absent <- setdiff(genes, rownames(zmat))
if (length(absent)) message("Not found in zmat (check symbol vs LOC ID): ", paste(absent, collapse = ", "))
genes <- intersect(genes, rownames(zmat))

# Class_Broad from the main analysis (Linear / Non-linear / Indistinguishable);
# genes that did not pass the pooled screen are labelled as such
class_lookup <- setNames(results$Class_Broad, results$Gene)

state_colors <- c("Normothermy" = "#7570B3", "Transition" = "#1B9E77", "Deep Torpor" = "#D95F02")

si <- sample_info %>%
  mutate(BirdID = factor(BirdID),
         X_c    = as.numeric(scale(X_raw)),
         Early  = as.numeric(State != "Normothermy"),
         Late   = as.numeric(State == "Deep Torpor"))

# same candidate models as the main analysis (bird random intercept, ML for AIC)
forms <- list(
  "Linear"      = y ~ X_c   + (1 | BirdID),
  "Early step"  = y ~ Early + (1 | BirdID),
  "Late step"   = y ~ Late  + (1 | BirdID),
  "State means" = y ~ State + (1 | BirdID)      # graded / non-monotonic
)

fit_gene <- function(g) {
  d <- si; d$y <- zmat[g, si$X]
  fits <- lapply(forms, function(f)
    suppressMessages(suppressWarnings(
      lmer(f, data = d, REML = FALSE, control = lmerControl(check.conv.singular = "ignore")))))
  aic  <- vapply(fits, AIC, numeric(1))
  best <- names(which.min(aic))
  gap  <- sort(aic)[2] - sort(aic)[1]            # AIC margin over the runner-up
  d$fit   <- predict(fits[[best]], re.form = NA) # fixed-effects curve (no bird effect)
  d$Gene  <- g
  cls <- unname(class_lookup[g])
  if (is.na(cls)) cls <- "Not responsive (pooled FDR >= 0.05)"
  d$label <- sprintf("%s | %s\nbest fit: %s (\u0394AIC to next = %.1f)", g, cls, best, gap)
  d
}

plot_df <- map_dfr(genes, fit_gene) %>%
  mutate(label = factor(label, levels = unique(label)))   # keep the order given in `genes`

# Trend line.
#  - Surface temperature: the states occupy separate x-ranges, so the fitted values
#    sorted by x already form a clean line (as in the earlier plot).
#  - VO2: the states overlap on x, so a line through the sorted fitted values zigzags
#    between state levels. Instead, smooth the model's fitted values over x: a local
#    average of what the best-fit model predicts at each VO2.
fit_layer <- if (predictor == "Surf_Temp") {
  geom_line(aes(y = fit), colour = "black", linewidth = 1)
} else {
  geom_smooth(aes(y = fit), method = "loess", formula = y ~ x, span = 0.4,
              method.args = list(degree = 1), se = FALSE,
              colour = "black", linewidth = 1)
}

# x-axis = the predictor used in the models (log10 VO2 when predictor = "log_VO2").
# To plot raw VO2 instead, use Mean_VO2 in place of X_raw below.
p <- ggplot(plot_df, aes(X_raw, y)) +
  geom_point(aes(colour = State, shape = Tissue), size = 2, alpha = 0.8) +
  fit_layer +
  facet_wrap(~ label, ncol = 5, scales = "free_y") +
  scale_colour_manual(values = state_colors) +
  scale_shape_manual(values = rep(c(16, 17, 15, 3, 7, 8, 18), length.out = nlevels(si$Tissue))) +
  labs(x = ifelse(predictor == "Surf_Temp", "Surface temperature (\u00B0C)", "log10 VO2 (ml/min)"),
       y = "Per-tissue Z-score",
       title = ifelse(predictor == "Surf_Temp",
                      "Expression vs temperature, with best-fit model (lowest AIC)",
                      "Expression vs VO2, with best-fit model (lowest AIC)"),
       subtitle = ifelse(predictor == "Surf_Temp",
                         "Points = samples (bird x tissue), coloured by state; black = fixed-effects fit; header shows Class_Broad from the main analysis",
                         "Points = samples (bird x tissue), coloured by state; black = best-fit model's fitted values, smoothed over VO2; header shows Class_Broad")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(out("selected_genes_scatter_bestfit.png"), p, width = 17, height = 8, dpi = 300)