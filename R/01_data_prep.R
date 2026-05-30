# =============================================================================
# 01_data_prep.R
# Build site-week and district-week tables from raw trap surveillance data
#
# INPUT:  SBD_17_23_xy   — raw trap data frame (loaded from data/raw/)
#                          Required columns:
#                            site_code, collection_date, disease_week,
#                            Total_Aedes, longitude, latitude
#
# OUTPUT: data/processed/site_week.rds
#         data/processed/district_week.rds
# =============================================================================

source(here::here("R", "00_config.R"))

# -----------------------------------------------------------------------------
# 1) Load raw data
# NOTE: Update the filename below to match your actual raw data file.
# -----------------------------------------------------------------------------
SBD_17_23_xy <- readRDS(file.path(data_raw_path, "SBD_17_23_xy.rds"))
# Alternatively: SBD_17_23_xy <- read.csv(file.path(data_raw_path, "SBD_17_23_xy.csv"))

# -----------------------------------------------------------------------------
# 2) Add CDC epidemiological week and year
# Uses lubridate::epiyear() and epiweek() which follow CDC/MMWR convention.
# disease_week is verified to match epiweek (sanity check below).
# -----------------------------------------------------------------------------
dat <- SBD_17_23_xy %>%
  mutate(
    wk_epiyear = epiyear(collection_date),
    wk_epiweek = epiweek(collection_date)
  )

# Sanity check: disease_week should match epiweek
tab <- dat %>%
  mutate(match = (as.integer(disease_week) == wk_epiweek)) %>%
  count(match)
print(tab)  # expect all TRUE

# Inspect mismatches if any
mism <- dat %>%
  filter(as.integer(disease_week) != wk_epiweek) %>%
  dplyr::select(collection_date, disease_week, wk_epiyear, wk_epiweek) %>%
  arrange(collection_date)
if (nrow(mism) > 0) warning("Found ", nrow(mism), " epiweek mismatches — inspect `mism`.")

# -----------------------------------------------------------------------------
# 3) Build SITE-WEEK table
# One row per unique (site_code × epi-year × epi-week).
# Columns: trap effort (trap_nights), total count, catch rate, week start date.
# -----------------------------------------------------------------------------
site_week <- dat %>%
  mutate(
    week_year = wk_epiyear,
    week_num  = as.integer(disease_week)
  ) %>%
  group_by(site_code, week_year, week_num) %>%
  summarise(
    Total_Aedes_week = sum(Total_Aedes, na.rm = TRUE),
    trap_nights      = n(),                              # 1 row = 1 trap-night
    longitude        = median(longitude, na.rm = TRUE),
    latitude         = median(latitude,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    week_start = MMWRweek2Date(MMWRyear = week_year, MMWRweek = week_num),
    rate_week  = Total_Aedes_week / pmax(trap_nights, 1L)
  ) %>%
  arrange(site_code, week_start)

# Add UTM coordinates and cyclic seasonality terms
site_week <- site_week %>%
  to_xy(lon = "longitude", lat = "latitude") %>%
  mutate(
    week_cc  = pmin(pmax(as.integer(week_num), 1L), 52L),  # clipped for cyclic spline
    week_sin = sin(2 * pi * week_num / 52),
    week_cos = cos(2 * pi * week_num / 52)
  )

# -----------------------------------------------------------------------------
# 4) Build DISTRICT-WEEK table
# Aggregates all sites; used for time-series visualization and district-level
# covariate summaries.
# -----------------------------------------------------------------------------
district_week <- dat %>%
  mutate(
    week_year = wk_epiyear,
    week_num  = as.integer(disease_week)
  ) %>%
  group_by(week_year, week_num) %>%
  summarise(
    Total_Aedes_week = sum(Total_Aedes, na.rm = TRUE),
    trap_nights      = n(),
    .groups = "drop"
  ) %>%
  mutate(
    yearweek   = sprintf("%04d-E%02d", week_year, week_num),
    offset_log = log(pmax(trap_nights, 1L)),
    rate_week  = Total_Aedes_week / trap_nights,
    week_start = MMWRweek2Date(MMWRyear = week_year, MMWRweek = week_num)
  ) %>%
  arrange(week_start)

# -----------------------------------------------------------------------------
# 5) Train / test split (train: 2017–2022; test: 2023)
# -----------------------------------------------------------------------------
train <- filter(site_week, week_year <= 2022)
test  <- filter(site_week, week_year == 2023)

cat("Training rows (2017–2022):", nrow(train), "\n")
cat("Test rows    (2023)      :", nrow(test),  "\n")

# -----------------------------------------------------------------------------
# 6) Quick diagnostic summaries
# -----------------------------------------------------------------------------
cat("\n--- site_week summary ---\n")
print(summary(site_week[, c("Total_Aedes_week", "trap_nights", "rate_week")]))

cat("\n--- district_week summary ---\n")
print(summary(district_week[, c("Total_Aedes_week", "trap_nights", "rate_week")]))

# -----------------------------------------------------------------------------
# 7) Save processed data
# -----------------------------------------------------------------------------
saveRDS(site_week,     file.path(data_proc_path, "site_week.rds"))
saveRDS(district_week, file.path(data_proc_path, "district_week.rds"))

message("✓ 01_data_prep.R complete — site_week and district_week saved.")
