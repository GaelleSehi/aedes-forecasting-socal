# =============================================================================
# 02_covariates.R
# Attach environmental covariates and autoregressive lags to site_week
#
# INPUT:  data/processed/site_week.rds
#         data/processed/district_week.rds
#         dat  — raw trap data with covariate columns (from 01_data_prep.R
#                or reloaded from raw)
#
# OUTPUT: data/processed/site_week_full.rds
#         data/processed/district_week_sc.rds
#
# Covariate columns expected in raw data:
#   NDWI14_scale  — 14-day lagged NDWI (normalized difference water index)
#   Psum14_scale  — 14-day cumulative precipitation (scaled)
#   Tlag14_scale  — 14-day lagged mean temperature (scaled)
#   DEM_scale     — digital elevation model (scaled)
#   Built_scale   — built environment index (scaled)
# =============================================================================

source(here::here("R", "00_config.R"))

site_week     <- readRDS(file.path(data_proc_path, "site_week.rds"))
district_week <- readRDS(file.path(data_proc_path, "district_week.rds"))

# Raw data must still be available for covariate aggregation
# (re-load here if not already in environment from 01_data_prep.R)
if (!exists("dat")) {
  SBD_17_23_xy <- readRDS(file.path(data_raw_path, "SBD_17_23_xy.rds"))
  dat <- SBD_17_23_xy %>%
    mutate(
      wk_epiyear = epiyear(collection_date),
      wk_epiweek = epiweek(collection_date)
    )
}

# -----------------------------------------------------------------------------
# 1) Site-week covariates
# Time-varying (weekly means) and static-ish (median across trap-nights)
# covariates aggregated to the site × week level.
# -----------------------------------------------------------------------------
site_week_cov <- dat %>%
  mutate(
    week_year = wk_epiyear,
    week_num  = as.integer(disease_week)
  ) %>%
  group_by(site_code, week_year, week_num) %>%
  summarise(
    NDWI14_week = mean(NDWI14_scale, na.rm = TRUE),   # surface water (lagged)
    Psum14_week = mean(Psum14_scale, na.rm = TRUE),   # precipitation (lagged)
    Tlag14_week = mean(Tlag14_scale, na.rm = TRUE),   # temperature (lagged)
    DEM_week    = median(DEM_scale,   na.rm = TRUE),  # elevation (static)
    Built_week  = median(Built_scale, na.rm = TRUE),  # built environment
    .groups = "drop"
  )

site_week <- site_week %>%
  left_join(site_week_cov, by = c("site_code", "week_year", "week_num"))

# -----------------------------------------------------------------------------
# 2) District-week covariates (effort-weighted mean across sites)
# -----------------------------------------------------------------------------
district_week_cov <- site_week %>%
  group_by(week_year, week_num) %>%
  summarise(
    NDWI14_week = w_mean(NDWI14_week, trap_nights),
    Psum14_week = w_mean(Psum14_week, trap_nights),
    Tlag14_week = w_mean(Tlag14_week, trap_nights),
    DEM_week    = w_mean(DEM_week,    trap_nights),
    Built_week  = w_mean(Built_week,  trap_nights),
    .groups = "drop"
  )

district_week <- district_week %>%
  left_join(district_week_cov, by = c("week_year", "week_num")) %>%
  mutate(
    week_sin = sin(2 * pi * week_num / 52),
    week_cos = cos(2 * pi * week_num / 52)
  )

# Same seasonality terms for site_week
site_week <- site_week %>%
  mutate(
    week_sin = sin(2 * pi * week_num / 52),
    week_cos = cos(2 * pi * week_num / 52)
  )

# -----------------------------------------------------------------------------
# 3) Re-scale weekly covariates using TRAINING years only (≤ 2022)
# This prevents data leakage: test-year scaling is based on train statistics.
# -----------------------------------------------------------------------------
dw_train <- dplyr::filter(district_week, week_year <= 2022)
dw_test  <- dplyr::filter(district_week, week_year == 2023)

# Compute mean and SD from training data only
sc_stats <- dw_train %>%
  summarise(
    NDWI_mu = mean(NDWI14_week, na.rm = TRUE), NDWI_sd = sd(NDWI14_week, na.rm = TRUE),
    Psum_mu = mean(Psum14_week, na.rm = TRUE), Psum_sd = sd(Psum14_week, na.rm = TRUE),
    Tlag_mu = mean(Tlag14_week, na.rm = TRUE), Tlag_sd = sd(Tlag14_week, na.rm = TRUE),
    DEM_mu  = mean(DEM_week,    na.rm = TRUE), DEM_sd  = sd(DEM_week,    na.rm = TRUE),
    Built_mu= mean(Built_week,  na.rm = TRUE), Built_sd= sd(Built_week,  na.rm = TRUE)
  )

scale_cols <- function(df, st) {
  df %>%
    mutate(
      NDWI_sc = (NDWI14_week - st$NDWI_mu) / pmax(st$NDWI_sd, .Machine$double.eps),
      Psum_sc = (Psum14_week - st$Psum_mu) / pmax(st$Psum_sd, .Machine$double.eps),
      Tlag_sc = (Tlag14_week - st$Tlag_mu) / pmax(st$Tlag_sd, .Machine$double.eps),
      DEM_sc  = (DEM_week    - st$DEM_mu ) / pmax(st$DEM_sd,  .Machine$double.eps),
      Built_sc= (Built_week  - st$Built_mu) / pmax(st$Built_sd,.Machine$double.eps)
    )
}

dw_train <- scale_cols(dw_train, sc_stats)
dw_test  <- scale_cols(dw_test,  sc_stats)

district_week_sc <- dplyr::bind_rows(dw_train, dw_test)

# -----------------------------------------------------------------------------
# 4) Autoregressive lags on log-rate scale (within-site, in time order)
# log1p(rate) lags at 1, 2, and 4 weeks.
#
# NOTE: lag4 causes ~46% missingness in 2023 evaluation rows and is
# EXCLUDED from the primary models. It is computed here for exploratory use.
# Primary models use log1p_lag1 and log1p_lag2 only.
# -----------------------------------------------------------------------------
site_week <- site_week %>%
  arrange(site_code, week_year, week_num) %>%
  group_by(site_code) %>%
  mutate(
    rate_t      = Total_Aedes_week / pmax(trap_nights, 1L),
    log1p_rate  = log1p(rate_t),
    log1p_lag1  = dplyr::lag(log1p_rate, 1),
    log1p_lag2  = dplyr::lag(log1p_rate, 2),
    log1p_lag4  = dplyr::lag(log1p_rate, 4)   # high missingness — use with caution
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 5) Covariate diagnostic checks
# -----------------------------------------------------------------------------
cat("\n--- Covariate NA profile (site_week) ---\n")
cov_cols <- c("NDWI14_week", "Psum14_week", "Tlag14_week", "DEM_week", "Built_week")
print(colSums(is.na(site_week[cov_cols])))

cat("\n--- AR lag NA proportions ---\n")
lag_cols <- c("log1p_lag1", "log1p_lag2", "log1p_lag4")
print(sapply(site_week[lag_cols], function(x) round(mean(is.na(x)), 3)))

# Verify lags are NA only at start of each site's time series
site_week %>%
  arrange(site_code, week_start) %>%
  group_by(site_code) %>%
  summarise(
    n       = n(),
    na_lag1 = sum(is.na(log1p_lag1)),
    na_lag2 = sum(is.na(log1p_lag2))
  ) %>%
  slice_head(n = 10) %>%
  print()

# -----------------------------------------------------------------------------
# 6) Save
# -----------------------------------------------------------------------------
saveRDS(site_week,        file.path(data_proc_path, "site_week_full.rds"))
saveRDS(district_week_sc, file.path(data_proc_path, "district_week_sc.rds"))
saveRDS(sc_stats,         file.path(data_proc_path, "scale_stats.rds"))

message("✓ 02_covariates.R complete — covariates and AR lags added.")
