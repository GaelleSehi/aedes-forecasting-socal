# =============================================================================
# 05_rolling_forecast.R
# Rolling 2-week-ahead forecasting — PRIMARY ANALYSIS
#
# For each week in 2023, models are re-trained on all data up to that week
# and used to predict mosquito abundance 2 weeks forward (horizon h=2).
# This is a strictly prospective evaluation with no data leakage.
#
# Baseline: persistence — carry forward the observed rate at time t
# Models:   GAM, RF, RFsp (same specifications as static analysis)
#
# OUTPUT:
#   data/processed/roll_results.rds      — site-level predictions all 2023 weeks
#   data/processed/metrics_rolling.rds   — overall performance metrics
#   output/metrics_rolling.csv
# =============================================================================

source(here::here("R", "00_config.R"))

site_week <- readRDS(file.path(data_proc_path, "site_week_full.rds"))
library(purrr)

# -----------------------------------------------------------------------------
# 1) Align for 2-week-ahead prediction
# Adds horizon-2 targets (count_lead2, rate_lead2) and current log-rate
# (y_lograte_t) used for the persistence baseline.
# -----------------------------------------------------------------------------
align_h2 <- function(df) {
  df %>%
    arrange(site_code, week_year, week_num) %>%
    group_by(site_code) %>%
    mutate(
      rate_t        = Total_Aedes_week / pmax(trap_nights, 1L),
      y_lograte_t   = log1p(rate_t),
      count_lead2   = dplyr::lead(Total_Aedes_week, 2),
      trapn_lead2   = dplyr::lead(trap_nights, 2),
      rate_lead2    = count_lead2 / pmax(trapn_lead2, 1L),
      lograte_lead2 = log1p(rate_lead2)
    ) %>%
    ungroup()
}

sw <- align_h2(site_week)

# -----------------------------------------------------------------------------
# 2) Model formulas and predictor sets
# lag4 excluded from all models due to ~46% missingness in evaluation rows.
# -----------------------------------------------------------------------------
gam_formula_h2 <- count_lead2 ~
  offset(log(pmax(trapn_lead2, 1L))) +
  s(X, Y, bs = "tp", k = 50) +
  DEM_week + Built_week +
  s(NDWI14_week, bs = "ts") +
  s(Psum14_week, bs = "ts") +
  s(Tlag14_week, bs = "ts", k = 10) +
  s(week_cc,     bs = "cc", k = 10) +
  s(week_year,   bs = "ts", k = 4) +
  ti(Tlag14_week, week_cc, bs = c("ts", "cc"), k = c(5, 10)) +
  log1p_lag1 + log1p_lag2

rf_vars <- c(
  "X", "Y",
  "NDWI14_week", "Psum14_week", "Tlag14_week",
  "DEM_week", "Built_week",
  "week_cc", "week_year",
  "log1p_lag1", "log1p_lag2"
)
p <- length(rf_vars)

# -----------------------------------------------------------------------------
# 3) Fit-and-predict function for a single rolling cut date
# -----------------------------------------------------------------------------
fit_and_predict_h2 <- function(train_df, eval_df) {

  # -- GAM --
  gam_train <- train_df %>%
    tidyr::drop_na(
      count_lead2, trapn_lead2,
      X, Y, DEM_week, Built_week,
      NDWI14_week, Psum14_week, Tlag14_week,
      week_cc, week_year, log1p_lag1, log1p_lag2
    )

  Gam_h2 <- mgcv::gam(
    gam_formula_h2,
    data      = gam_train,
    family    = mgcv::nb(link = "log"),
    method    = "REML",
    gamma     = 1.3,
    knots     = list(week_cc = c(0.5, 52.5)),
    select    = TRUE,
    na.action = na.exclude
  )

  # -- RF --
  rf_train <- train_df %>%
    select(all_of(c("rate_lead2", rf_vars))) %>%
    tidyr::drop_na()

  RF_h2 <- ranger::ranger(
    as.formula(paste("rate_lead2 ~", paste(rf_vars, collapse = " + "))),
    data          = rf_train,
    num.trees     = 110,
    mtry          = max(1L, floor(p / 3)),
    min.node.size = 5,
    importance    = "permutation"
  )

  # -- RFsp --
  rfsp_out   <- build_rfsp_features(
    train_df, eval_df,
    bandwidth    = 500,
    q            = 10,
    k_nn         = 8,
    moran_filter = TRUE
  )
  train_rfsp <- rfsp_out$train
  eval_rfsp  <- rfsp_out$test
  pcs        <- rfsp_out$kept_pcs
  rfsp_vars  <- c(pcs, rf_vars[!rf_vars %in% c("X", "Y")])
  p_sp       <- length(rfsp_vars)

  rfsp_train <- train_rfsp %>%
    select(all_of(c("rate_lead2", rfsp_vars))) %>%
    tidyr::drop_na()

  RFsp_h2 <- ranger::ranger(
    as.formula(paste("rate_lead2 ~", paste(rfsp_vars, collapse = " + "))),
    data          = rfsp_train,
    num.trees     = 120,
    mtry          = max(1L, floor(p_sp / 3)),
    min.node.size = 5,
    importance    = "permutation"
  )

  # -- Predictions --
  pred_gam_count <- predict(Gam_h2, newdata = eval_df, type = "response")
  pred_gam_rate  <- pred_gam_count / pmax(eval_df$trapn_lead2, 1L)
  pred_rf_rate   <- pmax(predict(RF_h2,  data = eval_df)$predictions, 0)

  keys <- c("site_code", "week_year", "week_num")
  eval_rfsp_pred <- eval_rfsp %>%
    dplyr::mutate(
      pred_rfsp_rate = pmax(predict(RFsp_h2, data = eval_rfsp)$predictions, 0)
    ) %>%
    dplyr::select(dplyr::all_of(keys), pred_rfsp_rate)

  out <- eval_df %>%
    dplyr::mutate(
      pred_gam_rate = pred_gam_rate,
      pred_rf_rate  = pred_rf_rate
    ) %>%
    dplyr::left_join(eval_rfsp_pred, by = keys)

  out
}

# -----------------------------------------------------------------------------
# 4) Rolling loop over all 2023 cut weeks
# At each cut_date t, train on all weeks ≤ t, predict rows at t
# whose t+2 truth is observed.
# -----------------------------------------------------------------------------
cut_weeks <- sw %>%
  filter(week_year == 2023) %>%
  distinct(week_start) %>%
  arrange(week_start) %>%
  pull(week_start)

cat("Running rolling forecast over", length(cut_weeks), "cut weeks...\n")

roll_results <- purrr::map_dfr(cut_weeks, function(cut_date) {

  train_df <- sw %>% filter(week_start <= cut_date)

  eval_df <- sw %>%
    filter(
      week_start == cut_date,
      is.finite(count_lead2), is.finite(trapn_lead2),
      is.finite(lograte_lead2), is.finite(y_lograte_t)
    ) %>%
    select(
      site_code, week_year, week_num, week_start,
      count_lead2, trapn_lead2, lograte_lead2, y_lograte_t,
      X, Y, DEM_week, Built_week,
      NDWI14_week, Psum14_week, Tlag14_week,
      week_cc, week_year, log1p_lag1, log1p_lag2
    )

  if (nrow(eval_df) == 0) return(tibble::tibble())

  preds <- fit_and_predict_h2(train_df, eval_df)
  preds %>% mutate(cut_date = cut_date)
})

cat("Rolling forecast complete.\n")
beepr::beep("mario")  # optional completion sound

# -----------------------------------------------------------------------------
# 5) Back-transform to rate scale and compute overall metrics
# -----------------------------------------------------------------------------
results <- roll_results %>%
  mutate(
    obs_rate_t2  = pmax(expm1(lograte_lead2), 0),   # observed at t+2
    base_rate_t2 = pmax(expm1(y_lograte_t),  0)     # persistence: rate at t
  )

mask <- with(results,
             is.finite(obs_rate_t2) & is.finite(base_rate_t2) &
               is.finite(pred_gam_rate) & is.finite(pred_rf_rate) &
               is.finite(pred_rfsp_rate))

results_aa <- results[mask, ]
mae_base   <- with(results_aa, mae(obs_rate_t2, base_rate_t2))

metrics_overall <- results_aa %>%
  transmute(
    obs  = obs_rate_t2,
    base = base_rate_t2,
    gam  = pred_gam_rate,
    rf   = pred_rf_rate,
    rfsp = pred_rfsp_rate
  ) %>%
  pivot_longer(-obs, names_to = "model", values_to = "pred") %>%
  group_by(model) %>%
  summarise(
    MAE  = mae(obs, pred),
    RMSE = rmse(obs, pred),
    R2   = r2(obs, pred),
    n    = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    Skill_vs_base = ifelse(model == "base", 0, 1 - (MAE / mae_base))
  )

cat("\n--- Rolling forecast metrics (2-week ahead, 2023) ---\n")
print(metrics_overall)

# -----------------------------------------------------------------------------
# 6) Save
# -----------------------------------------------------------------------------
saveRDS(results,         file.path(data_proc_path, "roll_results.rds"))
saveRDS(metrics_overall, file.path(data_proc_path, "metrics_rolling.rds"))
write.csv(metrics_overall, file.path(output_path, "metrics_rolling.csv"), row.names = FALSE)

message("✓ 05_rolling_forecast.R complete — results and metrics saved.")
