# Load libraries
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(lubridate)
library(here)
library(stringr)
library(zoo)
library(ggrepel)
library(patchwork)
library(effectsize)

# Set global ggplot theme with explicit BOLD settings for titles AND tick text
theme_set(
  theme_minimal(base_size = 16) +
    theme(
      axis.title        = element_text(face = "bold", size = 18, color = "black"),
      axis.title.x      = element_text(face = "bold", size = 18, color = "black"),
      axis.title.y      = element_text(face = "bold", size = 18, color = "black"),
      axis.text         = element_text(face = "bold", size = 14, color = "black"),
      axis.text.x       = element_text(face = "bold", size = 14, color = "black"),
      axis.text.y       = element_text(face = "bold", size = 14, color = "black"),
      plot.title        = element_text(face = "bold", size = 20, hjust = 0.5, color = "black"),
      strip.text        = element_text(face = "bold", size = 16, color = "black"),
      legend.title      = element_text(face = "bold", size = 16, color = "black"),
      legend.text       = element_text(face = "bold", size = 14, color = "black")
    )
)

# Set up main output directories relative to project root
main_fig_dir <- here("Outputs")
dir.create(main_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(main_fig_dir, "Individual_Plots"), showWarnings = FALSE, recursive = TRUE)

# List all analyzed files
all_files <- list.files(
  here("Hummingbird_MetabolicData"),
  pattern = "_analyzed\\.csv$",
  full.names = TRUE
)

# Define the blacklist based on lab notes (equipment unplugged, WVP/O2 freakouts)
bad_files <- c(
  "CAAN02_0623_2208", 
  "CAAN03_0624_2150", 
  "CAAN04_0628_2306", 
  "CAAN04_0628_2318", 
  "CAAN05_0629_2138", 
  "CAAN05_0629_2148", 
  "CAAN05_0629_2223", 
  "CAAN05_0629_2338"  
)

# Extract BirdID from file names
extract_bird_id <- function(file_path) {
  bird_id <- str_extract(basename(file_path), "^[A-Z0-9]+")
  tibble(file = file_path, BirdID = bird_id)
}

metadata_df <- map_dfr(all_files, extract_bird_id) %>%
  filter(!str_detect(file, paste(bad_files, collapse = "|")))

# Read and process individual CSV files
read_and_rebuild_datetime <- function(file, bird_id) {
  df <- read_csv(file, show_col_types = FALSE)
  
  if (!all(c("Year", "Month", "Day", "Time_hours", "VO2_ml_min") %in% names(df))) {
    return(NULL)
  }
  
  df %>%
    mutate(
      BirdID = bird_id,
      Time_hours = as.numeric(Time_hours),
      datetime = make_datetime(Year, Month, Day) + seconds(Time_hours * 3600)
    ) %>%
    dplyr::select(BirdID, datetime, VO2_ml_min)
}
all_data <- map2_dfr(metadata_df$file, metadata_df$BirdID, read_and_rebuild_datetime)

# Join metabolic state metadata
metadata <- read_csv(here("Hummingbird_MetabolicData", "RNASeq_metadata.csv"), show_col_types = FALSE) %>%
  dplyr::select(BirdID, Metabolic_State) %>%
  distinct()

all_data <- all_data %>%
  left_join(metadata, by = "BirdID") %>%
  mutate(Metabolic_State = recode(Metabolic_State,
                                  "N" = "Normothermy",
                                  "D" = "Deep torpor",
                                  "T" = "Transition"))

# Directly define euthanasia timepoints
euth_events <- tibble(
  BirdID = c("CAAN01", "CAAN02", "CAAN03", "CAAN04", "CAAN05", "CAAN06",
             "CAAN07", "CAAN08", "CAAN09", "CAAN10", "CAAN11", "CAAN12",
             "CAAN13", "CAAN14", "CAAN15", "CAAN16", "CAAN17", "CAAN18"),
  euth_datetime = as.POSIXct(c(
    "2021-06-22 01:30:00", "2021-06-24 01:02:00", "2021-06-25 02:11:00",
    "2021-06-29 00:46:00", "2021-06-30 00:50:00", "2021-07-03 00:38:45",
    "2021-07-06 22:57:00", "2021-07-08 02:40:00", "2021-07-09 00:48:00",
    "2021-07-10 01:59:00", "2021-07-13 00:54:00", "2021-07-14 01:39:00",
    "2021-07-19 23:06:00", "2021-07-20 23:37:00", "2021-07-21 22:35:00",
    "2021-07-22 23:13:00", "2021-07-26 23:33:00", "2021-07-27 23:42:00"
  ), tz = "UTC")
)

summary_times <- all_data %>%
  group_by(BirdID) %>%
  summarise(
    start_time = min(datetime, na.rm = TRUE),
    end_time = max(datetime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(euth_events, by = "BirdID")

# Plot with euthanasia markers
ggplot(all_data, aes(x = datetime, y = VO2_ml_min, color = Metabolic_State)) +
  geom_line(alpha = 0.8) +
  geom_vline(
    data = euth_events,
    mapping = aes(xintercept = as.numeric(euth_datetime)),
    color = "red", linetype = "dashed", size = 0.8
  ) +
  facet_wrap(~ BirdID, scales = "free_x") +
  labs(
    x = "Absolute Time",
    y = expression(bold(VO[2]~"(ml/min)")),
    title = "Metabolic Rate with Euthanasia Timepoints (Red Dashed Lines)",
    color = "Metabolic State"
  ) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", size = 18, color = "black"),
    axis.text = element_text(face = "bold", size = 14, color = "black")
  )

ggsave(file.path(main_fig_dir, "Metabolic rate vs Time (Euthenasia markers).tiff"),
       dpi = 300, width = 16, height = 12, units = "in", compression = "lzw")

# Extract last 10 minutes of data
last10_df <- all_data %>%
  group_by(BirdID) %>%
  arrange(datetime) %>%
  slice_tail(n = 2400) %>%
  ungroup()

summary_stats <- last10_df %>%
  group_by(BirdID, Metabolic_State) %>%
  summarise(
    mean_VO2 = mean(VO2_ml_min, na.rm = TRUE),
    median_VO2 = median(VO2_ml_min, na.rm = TRUE),
    sd_VO2 = sd(VO2_ml_min, na.rm = TRUE),
    .groups = "drop")

# Boxplot colored by metabolic state
ggplot(last10_df, aes(x = BirdID, y = VO2_ml_min, fill = Metabolic_State)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    x = "Bird ID",
    y = expression(bold(VO[2]~"(ml/min)")),
    title = "VO2 (ml/min) in Last 10 Minutes Before Euthanasia",
    fill = "Metabolic State"
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 18, color = "black"),
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14, color = "black")
  )

ggsave(file.path(main_fig_dir, "VO2_mlPerMin vs BirdID.tiff"),
       dpi = 300, width = 8, height = 6, units = "in", compression = "lzw")

last10_df <- last10_df %>%
  mutate(Metabolic_State = factor(Metabolic_State, levels = c("Normothermy", "Transition", "Deep torpor")))

last10_summary <- last10_df %>%
  group_by(BirdID, Metabolic_State) %>%
  summarise(mean_VO2 = mean(VO2_ml_min, na.rm = TRUE), .groups = "drop")

surface_temp <- read_csv(here("Hummingbird_MetabolicData", "Euthenasia.csv"), show_col_types = FALSE) %>%
  select(BirdID, Surf_Temp) %>%
  distinct()

last10_summary <- merge(last10_summary, surface_temp, by="BirdID")
last10_summary <- merge(last10_summary, summary_times[c("BirdID", "end_time")], by ="BirdID")

write_csv(last10_summary, here("Outputs", "metadata_surfTemp_mlO2_min_euthTime.csv"))

surface_temp_state_BirdID <- merge(metadata, surface_temp, by="BirdID") %>%
  mutate(Metabolic_State = recode(Metabolic_State,
                                  "D" = "Deep Torpor",
                                  "N" = "Normothermy",
                                  "T" = "Transition")) %>%
  mutate(Metabolic_State = factor(Metabolic_State, levels = c("Normothermy", "Transition", "Deep Torpor")))

# Linear models and ANOVAs
lm_vo2_temp <- lm(mean_VO2 ~ Surf_Temp, data = last10_summary)
lm_vo2_temp_state <- lm(mean_VO2 ~ Surf_Temp + Metabolic_State, data = last10_summary)

anova_vo2 <- aov(VO2_ml_min ~ Metabolic_State, data = last10_df)
last10_df <- last10_df %>% left_join(surface_temp, by = "BirdID")
anova_temp <- aov(Surf_Temp ~ Metabolic_State, data = last10_df)


# ==============================================================================
# FIGURE 1 GENERATION (2x2 GRID VIA PATCHWORK)
# ==============================================================================

# Panel A: CAAN16 Trace
bird_id <- "CAAN16"
bird_data <- all_data %>% filter(BirdID == bird_id)
bird_euth <- euth_events %>% filter(BirdID == bird_id)

p1 <- ggplot(bird_data, aes(x = datetime, y = VO2_ml_min)) +
  geom_line(color = "purple", alpha = 0.6, size = 1) +
  geom_vline(
    data = bird_euth,
    aes(xintercept = as.numeric(euth_datetime)),
    color = "goldenrod3", linetype = "dashed", size = 0.8
  ) +
  labs(
    x = "Absolute Time",
    y = expression(bold(VO[2]~"(ml/min)")),
    title = "Sample metabolic trace of a bird entering Torpor"
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 18, color = "black"),
    axis.text  = element_text(face = "bold", size = 14, color = "black")
  )

# Panel B: Surface Temp Boxplot
st_filtered <- surface_temp_state_BirdID %>% 
  filter(Metabolic_State %in% c("Normothermy", "Deep Torpor"))

p2 <- ggplot(st_filtered, aes(x = Metabolic_State, y = Surf_Temp, fill = Metabolic_State)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, height = 0, size = 2.5, alpha = 0.6, color = "black") +
  labs(
    title = "Surface Temperature by Physiological State",
    x = "Physiological State",
    y = "Surface Temperature (°C)"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold", size = 18, color = "black"),
    axis.text  = element_text(face = "bold", size = 14, color = "black")
  ) +
  scale_fill_manual(values = c("Normothermy" = "#1b9e77", "Deep Torpor" = "#7570b3"))

# Panel C: Regression
p3 <- ggplot(last10_summary, aes(x = Surf_Temp, y = mean_VO2)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Metabolic Rate vs Surface Temperature",
    x = "Surface Temperature (°C)",
    y = expression(bold(Mean~VO[2]~"(ml/min)"))
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 18, color = "black"),
    axis.text  = element_text(face = "bold", size = 14, color = "black")
  )

# Panel D: VO2 Boxplot
vo2_filtered_df <- last10_df %>% filter(Metabolic_State %in% c("Normothermy", "Deep torpor"))
vo2_filtered_summary <- last10_summary %>% filter(Metabolic_State %in% c("Normothermy", "Deep torpor"))

p4 <- ggplot(vo2_filtered_df, aes(x = Metabolic_State, y = VO2_ml_min, fill = Metabolic_State)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(
    data = vo2_filtered_summary,
    aes(x = Metabolic_State, y = mean_VO2),
    color = "black", size = 2
  ) + 
  labs(
    x = "Physiological State",
    y = expression(bold(Mean~VO[2]~"(ml/min)")),
    title = "Metabolic rate by Physiological State",
    fill = "Metabolic State"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold", size = 18, color = "black"),
    axis.text  = element_text(face = "bold", size = 14, color = "black")
  ) + 
  scale_fill_manual(values = c("Normothermy" = "#1b9e77", "Deep torpor" = "#7570b3"))

# Combine into 2x2 Grid using Patchwork and apply theme globally across subplots
fig1_grid <- ((p1 | p2) / (p3 | p4)) & 
  theme(
    axis.title   = element_text(face = "bold", size = 18, color = "black"),
    axis.title.x = element_text(face = "bold", size = 18, color = "black"),
    axis.title.y = element_text(face = "bold", size = 18, color = "black"),
    axis.text    = element_text(face = "bold", size = 14, color = "black"),
    axis.text.x  = element_text(face = "bold", size = 14, color = "black"),
    axis.text.y  = element_text(face = "bold", size = 14, color = "black"),
    plot.title   = element_text(face = "bold", size = 20, hjust = 0.5, color = "black"),
    plot.tag     = element_text(face = "bold", size = 22, color = "black")
  )

fig1_grid <- fig1_grid + plot_annotation(tag_levels = 'A')

ggsave(file.path(main_fig_dir, "Figure1_26thAugust2026.tiff"),
       plot = fig1_grid,
       dpi = 300,
       width = 16,
       height = 12,
       units = "in",
       compression = "lzw"
)


######## Individual plots for each bird #####
bird_list <- unique(all_data$BirdID)

walk(bird_list, function(id) {
  bird_df <- all_data %>% filter(BirdID == id)
  bird_euth <- euth_events %>% filter(BirdID == id)
  
  p <- ggplot(bird_df, aes(x = datetime, y = VO2_ml_min)) +
    geom_line(color = "black", alpha = 0.8, size = 0.7) +
    geom_vline(data = bird_euth, 
               aes(xintercept = as.numeric(euth_datetime)), 
               color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Metabolic Rate Trace -", id),
      subtitle = paste("State:", unique(bird_df$Metabolic_State)),
      x = "Time",
      y = expression(bold(VO[2]~"(ml/min)"))
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 18, color = "black"),
      axis.text  = element_text(face = "bold", size = 14, color = "black")
    )
  
  file_name <- file.path(main_fig_dir, "Individual_Plots", paste0("Trace_", id, ".tiff"))
  ggsave(file_name, plot = p, width = 10, height = 6, dpi = 300, compression = "lzw")
})


# ==============================================================================
# DUAL-AXIS TORPOR KINETICS
# ==============================================================================
deep_torpor_birds <- c("CAAN02", "CAAN03", "CAAN04", "CAAN05", "CAAN14", "CAAN16")

all_data_smooth <- all_data %>%
  filter(BirdID %in% deep_torpor_birds) %>%
  group_by(BirdID) %>%
  arrange(datetime) %>%
  mutate(
    time_diff_mins = as.numeric(difftime(datetime, lag(datetime), units = "mins")),
    gap_flag = if_else(is.na(time_diff_mins) | time_diff_mins > 5, 1, 0),
    segment_id = cumsum(gap_flag)
  ) %>%
  group_by(BirdID, segment_id) %>%
  mutate(
    VO2_smooth = zoo::rollmean(VO2_ml_min, k = 1200, fill = NA, align = "center")
  ) %>%
  ungroup()

temp_data <- read_csv(here("Hummingbird_MetabolicData", "Torpor_data_2021.xlsx - IR_notes.csv"), show_col_types = FALSE)

temp_clean <- temp_data %>%
  filter(BirdID %in% deep_torpor_birds) %>%
  mutate(Ts_max = as.numeric(Ts_max), Tamb = as.numeric(Tamb)) %>% 
  filter(!is.na(Ts_max)) %>%
  mutate(
    hour_raw = Time %/% 100,
    min_raw  = Time %% 100,
    is_next_day = (hour_raw >= 24) | (hour_raw < 12),
    hour_clean = ifelse(hour_raw >= 24, hour_raw - 24, hour_raw),
    datetime_base = make_datetime(Year, Month, Day, hour = hour_clean, min = min_raw),
    datetime = if_else(is_next_day, datetime_base + days(1), datetime_base)
  ) %>%
  dplyr::select(BirdID, datetime, Ts_max, Tamb) %>%
  arrange(BirdID, datetime)

analyze_torpor_phases <- function(df, target_bird) {
  bird_df <- df %>% filter(BirdID == target_bird) %>% arrange(datetime)
  if(nrow(bird_df) < 10) return(NULL)
  
  start_time <- min(bird_df$datetime)
  bird_df$t_min <- as.numeric(difftime(bird_df$datetime, start_time, units = "mins"))
  
  fit <- smooth.spline(bird_df$t_min, bird_df$Ts_max, spar = 0.5)
  t_seq <- seq(min(bird_df$t_min), max(bird_df$t_min), length.out = 300)
  pred <- predict(fit, t_seq, deriv = 1) 
  
  slope_df <- tibble(time_min = pred$x, slope = pred$y, datetime = start_time + (pred$x * 60))
  
  min_slope_idx <- which.min(slope_df$slope)
  min_slope <- slope_df$slope[min_slope_idx]
  threshold <- min_slope * 0.20
  
  pre_peak_slopes <- slope_df$slope[1:min_slope_idx]
  start_idx <- max(which(pre_peak_slopes > threshold))
  if(is.infinite(start_idx)) start_idx <- 1
  entry_time <- slope_df$datetime[start_idx]
  
  completion_df <- bird_df %>% filter(datetime > entry_time) %>% mutate(temp_diff = Ts_max - Tamb) %>% filter(temp_diff <= 3) 
  
  if(nrow(completion_df) > 0) {
    completion_time <- min(completion_df$datetime)
    is_complete <- TRUE
  } else {
    completion_time <- NA
    is_complete <- FALSE
  }
  
  tibble(BirdID = target_bird, Entry_Time = entry_time, Completion_Time = completion_time, Is_Complete = is_complete)
}

torpor_phases <- map_dfr(deep_torpor_birds, ~analyze_torpor_phases(temp_clean, .x))

coeff <- 20 
cols <- c("Raw VO2"="gray70", "Smoothed VO2"="darkgreen", "Surface Temp"="goldenrod3", "Ambient Temp"="black", "Entry (Start)"="purple", "Completion (Diff <= 3°C)"="purple")

ggplot() +
  geom_line(data = all_data_smooth, aes(x = datetime, y = VO2_ml_min, color = "Raw VO2", group = interaction(BirdID, segment_id)), size = 0.3, alpha = 0.5) +
  geom_line(data = all_data_smooth, aes(x = datetime, y = VO2_smooth, color = "Smoothed VO2", group = interaction(BirdID, segment_id)), size = 1) +
  geom_line(data = temp_clean, aes(x = datetime, y = Tamb / coeff, color = "Ambient Temp"), linetype = "dotted", size = 0.8) +
  geom_line(data = temp_clean, aes(x = datetime, y = Ts_max / coeff, color = "Surface Temp"), size = 0.8) +
  geom_point(data = temp_clean, aes(x = datetime, y = Ts_max / coeff, color = "Surface Temp"), size = 1.5) +
  geom_vline(data = torpor_phases, aes(xintercept = as.numeric(Entry_Time), color = "Entry (Start)"), linetype = "dashed", size = 1) +
  geom_vline(data = torpor_phases, aes(xintercept = as.numeric(Completion_Time), color = "Completion (Diff <= 3°C)"), linetype = "solid", size = 1, na.rm = TRUE) +
  facet_wrap(~ BirdID, scales = "free_x") +
  scale_y_continuous(name = expression(bold(VO[2]~"(ml/min)")), sec.axis = sec_axis(~ . * coeff, name = "Temperature (°C)")) +
  scale_color_manual(name = "Legend", values = cols) +
  labs(title = "Torpor Kinetics: Improved Entry Detection", x = "Time") +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(face = "bold", size = 18, color = "black"),
    axis.text.x  = element_text(face = "bold", size = 14, color = "black"),
    axis.title.y.left  = element_text(color = "darkgreen", face = "bold", size = 18),
    axis.text.y.left   = element_text(color = "darkgreen", face = "bold", size = 14),
    axis.title.y.right = element_text(color = "goldenrod3", face = "bold", size = 18),
    axis.text.y.right  = element_text(color = "goldenrod3", face = "bold", size = 14),
    panel.grid.minor   = element_blank()
  )

ggsave(file.path(main_fig_dir, "Torpor_Kinetics_ImprovedEntry_Clean.tiff"), width = 14, height = 9, bg = "white", compression = "lzw")

# Calculate Durations
cooling_duration_df <- torpor_phases %>%
  mutate(Cooling_Duration_Mins = round(as.numeric(difftime(Completion_Time, Entry_Time, units = "mins")), 1)) %>%
  dplyr::select(BirdID, Entry_Time, Completion_Time, Cooling_Duration_Mins)

write_csv(cooling_duration_df, here("Outputs", "Torpor_Cooling_Durations.csv"))

ggplot(cooling_duration_df, aes(x = BirdID, y = Cooling_Duration_Mins, fill = "gray50")) +
  geom_col(alpha = 0.8, color = "black") +
  geom_text(aes(label = paste0(Cooling_Duration_Mins, " min")), vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Time to Torpor Entry", subtitle = "Duration of Cooling Phase", y = "Duration (Minutes)", x = "Bird ID") +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold", size = 18, color = "black"),
    axis.text  = element_text(face = "bold", size = 14, color = "black")
  )

ggsave(file.path(main_fig_dir, "Torpor_Cooling_Duration_Barplot.tiff"), width = 8, height = 6, bg = "white", compression = "lzw")


# ==============================================================================
# STATISTICAL COMPARISONS: NORMOTHERMY VS. DEEP TORPOR
# ==============================================================================

# Filter dataset to include strictly Normothermy and Deep Torpor groups
stat_df <- last10_summary %>%
  filter(Metabolic_State %in% c("Normothermy", "Deep torpor", "Deep Torpor")) %>%
  mutate(Metabolic_State = factor(Metabolic_State, levels = c("Normothermy", "Deep torpor")))

# ------------------------------------------------------------------------------
# 1. SUMMARY STATISTICS
# ------------------------------------------------------------------------------
summary_stats_group <- stat_df %>%
  group_by(Metabolic_State) %>%
  summarise(
    N           = n(),
    Mean_VO2    = mean(mean_VO2, na.rm = TRUE),
    SD_VO2      = sd(mean_VO2, na.rm = TRUE),
    Median_VO2  = median(mean_VO2, na.rm = TRUE),
    IQR_VO2     = IQR(mean_VO2, na.rm = TRUE),
    Mean_Temp   = mean(Surf_Temp, na.rm = TRUE),
    SD_Temp     = sd(Surf_Temp, na.rm = TRUE),
    Median_Temp = median(Surf_Temp, na.rm = TRUE),
    IQR_Temp    = IQR(Surf_Temp, na.rm = TRUE),
    .groups     = "drop"
  )

cat("\n================ GROUP SUMMARY STATISTICS ================\n")
print(summary_stats_group)

# ------------------------------------------------------------------------------
# 2. METABOLIC RATE (VO2) COMPARISONS
# ------------------------------------------------------------------------------
# Welch's two-sample t-test (Unpaired, unequal variances)
t_test_vo2 <- t.test(mean_VO2 ~ Metabolic_State, data = stat_df, var.equal = FALSE)

# Non-parametric Wilcoxon Rank Sum (Mann-Whitney U) test
wilcox_vo2 <- wilcox.test(mean_VO2 ~ Metabolic_State, data = stat_df, exact = FALSE)

# Cohen's d (Effect Size using effectsize package)
cohen_vo2  <- cohens_d(mean_VO2 ~ Metabolic_State, data = stat_df)

# ------------------------------------------------------------------------------
# 3. SURFACE TEMPERATURE COMPARISONS
# ------------------------------------------------------------------------------
# Welch's two-sample t-test (Unpaired, unequal variances)
t_test_temp <- t.test(Surf_Temp ~ Metabolic_State, data = stat_df, var.equal = FALSE)

# Non-parametric Wilcoxon Rank Sum (Mann-Whitney U) test
wilcox_temp <- wilcox.test(Surf_Temp ~ Metabolic_State, data = stat_df, exact = FALSE)

# Cohen's d (Effect Size using effectsize package)
cohen_temp  <- cohens_d(Surf_Temp ~ Metabolic_State, data = stat_df)

# ------------------------------------------------------------------------------
# 4. PRINT FORMATTED STATISTICAL RESULTS TO CONSOLE
# ------------------------------------------------------------------------------
cat("\n================ METABOLIC RATE (VO2) STATISTICAL RESULTS ================\n")
cat(sprintf("Welch's t-test    : t = %.3f, df = %.2f, p-value = %.5e\n", 
            t_test_vo2$statistic, t_test_vo2$parameter, t_test_vo2$p.value))
cat(sprintf("Mann-Whitney U    : W = %.1f, p-value = %.5e\n", 
            wilcox_vo2$statistic, wilcox_vo2$p.value))
cat(sprintf("Cohen's d         : d = %.3f [95%% CI: %.3f, %.3f]\n", 
            cohen_vo2$Cohens_d, cohen_vo2$CI_low, cohen_vo2$CI_high))

cat("\n================ SURFACE TEMPERATURE STATISTICAL RESULTS ================\n")
cat(sprintf("Welch's t-test    : t = %.3f, df = %.2f, p-value = %.5e\n", 
            t_test_temp$statistic, t_test_temp$parameter, t_test_temp$p.value))
cat(sprintf("Mann-Whitney U    : W = %.1f, p-value = %.5e\n", 
            wilcox_temp$statistic, wilcox_temp$p.value))
cat(sprintf("Cohen's d         : d = %.3f [95%% CI: %.3f, %.3f]\n", 
            cohen_temp$Cohens_d, cohen_temp$CI_low, cohen_temp$CI_high))

# ------------------------------------------------------------------------------
# 5. COMPILE AND EXPORT STATISTICAL SUMMARY TABLE
# ------------------------------------------------------------------------------
stat_summary_table <- tibble(
  Variable           = c("Metabolic Rate (VO2)", "Surface Temperature"),
  Normothermy_Mean   = c(summary_stats_group$Mean_VO2[summary_stats_group$Metabolic_State == "Normothermy"],
                         summary_stats_group$Mean_Temp[summary_stats_group$Metabolic_State == "Normothermy"]),
  Normothermy_SD     = c(summary_stats_group$SD_VO2[summary_stats_group$Metabolic_State == "Normothermy"],
                         summary_stats_group$SD_Temp[summary_stats_group$Metabolic_State == "Normothermy"]),
  Torpor_Mean        = c(summary_stats_group$Mean_VO2[summary_stats_group$Metabolic_State == "Deep torpor"],
                         summary_stats_group$Mean_Temp[summary_stats_group$Metabolic_State == "Deep torpor"]),
  Torpor_SD          = c(summary_stats_group$SD_VO2[summary_stats_group$Metabolic_State == "Deep torpor"],
                         summary_stats_group$SD_Temp[summary_stats_group$Metabolic_State == "Deep torpor"]),
  Welch_t_stat       = c(t_test_vo2$statistic, t_test_temp$statistic),
  Welch_df           = c(t_test_vo2$parameter, t_test_temp$parameter),
  Welch_p_value      = c(t_test_vo2$p.value, t_test_temp$p.value),
  Wilcoxon_W         = c(wilcox_vo2$statistic, wilcox_temp$statistic),
  Wilcoxon_p_value   = c(wilcox_vo2$p.value, wilcox_temp$p.value),
  Cohens_d           = c(cohen_vo2$Cohens_d, cohen_temp$Cohens_d)
)

write_csv(stat_summary_table, file.path(main_fig_dir, "Normothermy_vs_Torpor_Stats_Summary.csv"))
cat("\nStatistical summary exported to: Normothermy_vs_Torpor_Stats_Summary.csv\n")