# =============================================================================
# 04_static_models.R
# Static forecasting: train on 2017–2022, evaluate on 2023
#
# This is a "static" (non-rolling) analysis where all training data are used
# at once and the model is applied once to held-out 2023 data.
# The primary rolling analysis is in 05_rolling_forecast.R.
#
# Models:
#   Baseline   — seasonal site × week-of-year climatology
#   GAM        — negative binomial GAM with counts + log(offset)
#   RF         — Random Forest predicting rate per trap-night
#   RFsp       — RF with kernel PCA spatial features (RFsp)
#
# OUTPUT:
#   data/processed/pred_2023_static.rds
#   figures/Figure2_static_timeseries.pdf
#   output/metrics_static_weekly.csv
#   output/metrics_static_siteweek.csv
# =============================================================================

source(here::here("R", "00_config.R"))

site_week <- readRDS(file.path(data_proc_path, "site_week_full.rds"))

# -----------------------------------------------------------------------------
# 1) Static train/test split
# -----------------------------------------------------------------------------
train <- site_week %>% filter(week_year <= 2022)
test  <- site_week %>% filter(week_year == 2023)

# -----------------------------------------------------------------------------
# 2) Baseline: site × week-of-year seasonal climatology
# Primary: mean rate for each (site_code, week_cc) in training years.
# Fallback: district mean for that week_cc if site has no history.
# -----------------------------------------------------------------------------
base_site_week   <- train %>%
  group_by(site_code, week_cc) %>%
  summarise(pred_base_rate = mean(rate_t, na.rm = TRUE), .groups = "drop")

base_global_week <- train %>%
  group_by(week_cc) %>%
  summarise(pred_base_rate_global = mean(rate_t, na.rm = TRUE), .groups = "drop")

test <- test %>%
  left_join(base_site_week,   by = c("site_code", "week_cc")) %>%
  left_join(base_global_week, by = "week_cc") %>%
  mutate(pred_base_rate = dplyr::coalesce(pred_base_rate, pred_base_rate_global))

# -----------------------------------------------------------------------------
# 3) GAM: negative binomial, counts + log(trap_nights) offset
# Smooths: 2D spatial (tp), seasonal (cyclic cc), inter-annual trend,
# environmental covariates, temperature × season interaction.
# Linear: parametric AR lags, DEM, built environment.
# -----------------------------------------------------------------------------
gam_formula_static <- Total_Aedes_week ~
  offset(log(pmax(trap_nights, 1L))) +
  s(X, Y, bs = "tp", k = 50) +
  DEM_week + Built_week +
  s(NDWI14_week, bs = "ts") +
  s(Psum14_week, bs = "ts") +
  s(Tlag14_week, bs = "ts", k = 10) +
  s(week_cc,     bs = "cc", k = 10) +
  ti(Tlag14_week, week_cc, bs = c("ts", "cc"), k = c(5, 10)) +
  log1p_lag1 + log1p_lag2

gam_train <- train %>%
  tidyr::drop_na(
    Total_Aedes_week, trap_nights,
    X, Y, DEM_week, Built_week,
    NDWI14_week, Psum14_week, Tlag14_week, week_cc,
    log1p_lag1, log1p_lag2
  )

cat("GAM training rows:", nrow(gam_train), "\n")

Gam_static <- mgcv::gam(
  gam_formula_static,
  data      = gam_train,
  family    = mgcv::nb(link = "log"),
  method    = "REML",
  gamma     = 1.3,
  knots     = list(week_cc = c(0.5, 52.5)),
  select    = TRUE,
  na.action = na.exclude
)

cat("\n--- GAM summary ---\n")
print(summary(Gam_static))

# Predict counts then convert to rate
test_pred_gam_count <- predict(Gam_static, newdata = test, type = "response")
test_pred_gam_rate  <- test_pred_gam_count / pmax(test$trap_nights, 1L)

# -----------------------------------------------------------------------------
# 4) RF: Random Forest predicting rate per trap-night
# Spatial position included as raw UTM coordinates (X, Y).
# -----------------------------------------------------------------------------
rf_vars <- c(
  "X", "Y",
  "NDWI14_week", "Psum14_week", "Tlag14_week",
  "DEM_week", "Built_week",
  "week_cc", "log1p_lag1", "log1p_lag2"
)

rf_train <- train %>%
  select(rate_t, all_of(rf_vars)) %>%
  tidyr::drop_na()

cat("RF training rows:", nrow(rf_train), "\n")

RF_static <- ranger::ranger(
  as.formula(paste("rate_t ~", paste(rf_vars, collapse = " + "))),
  data          = rf_train,
  num.trees     = 110,
  mtry          = max(1L, floor(length(rf_vars) / 3)),
  min.node.size = 5,
  importance    = "permutation"
)

test_pred_rf_rate <- pmax(predict(RF_static, data = test)$predictions, 0)

# -----------------------------------------------------------------------------
# 5) RFsp: RF with kernel PCA spatial features replacing raw X, Y
# -----------------------------------------------------------------------------
rfsp_out <- build_rfsp_features(
  train, test,
  bandwidth    = 500,
  q            = 10,
  k_nn         = 8,
  moran_filter = TRUE
)

train_rfsp <- rfsp_out$train
test_rfsp  <- rfsp_out$test
pcs        <- rfsp_out$kept_pcs

rfsp_vars  <- c(pcs, rf_vars[!rf_vars %in% c("X", "Y")])
p_sp       <- length(rfsp_vars)

rfsp_train <- train_rfsp %>%
  select(rate_t, all_of(rfsp_vars)) %>%
  tidyr::drop_na()

cat("RFsp training rows:", nrow(rfsp_train), "\n")

RFsp_static <- ranger::ranger(
  as.formula(paste("rate_t ~", paste(rfsp_vars, collapse = " + "))),
  data          = rfsp_train,
  num.trees     = 120,
  mtry          = max(1L, floor(p_sp / 3)),
  min.node.size = 5,
  importance    = "permutation"
)

test_pred_rfsp_rate <- pmax(predict(RFsp_static, data = test_rfsp)$predictions, 0)

# -----------------------------------------------------------------------------
# 6) Assemble predictions for 2023
# -----------------------------------------------------------------------------
pred_2023 <- test %>%
  mutate(
    obs_rate       = rate_t,
    pred_base_rate = pred_base_rate,
    pred_gam_rate  = test_pred_gam_rate,
    pred_rf_rate   = test_pred_rf_rate,
    pred_rfsp_rate = test_pred_rfsp_rate
  )

saveRDS(pred_2023, file.path(data_proc_path, "pred_2023_static.rds"))

# -----------------------------------------------------------------------------
# 7) Error metrics — (A) district weekly mean and (B) site-week level
# -----------------------------------------------------------------------------

## A) Weekly mean
weekly_wide <- pred_2023 %>%
  group_by(week_start) %>%
  summarise(
    obs_rate       = mean(obs_rate,       na.rm = TRUE),
    pred_base_rate = mean(pred_base_rate, na.rm = TRUE),
    pred_gam_rate  = mean(pred_gam_rate,  na.rm = TRUE),
    pred_rf_rate   = mean(pred_rf_rate,   na.rm = TRUE),
    pred_rfsp_rate = mean(pred_rfsp_rate, na.rm = TRUE),
    .groups = "drop"
  )

metrics_weeklymean <- weekly_wide %>%
  pivot_longer(cols = c(pred_base_rate, pred_gam_rate, pred_rf_rate, pred_rfsp_rate),
               names_to = "model", values_to = "pred") %>%
  mutate(model = recode(model,
    pred_base_rate = "Persistence", pred_gam_rate = "GAM",
    pred_rf_rate   = "RF",          pred_rfsp_rate = "RFsp"
  )) %>%
  group_by(model) %>%
  summarise(
    MAE  = mae(obs_rate,  pred),
    RMSE = rmse(obs_rate, pred),
    Bias = bias(obs_rate, pred),
    R2   = r2(obs_rate,   pred),
    n    = sum(is.finite(obs_rate) & is.finite(pred)),
    .groups = "drop"
  )

## B) Site-week
metrics_siteweek <- pred_2023 %>%
  pivot_longer(cols = c(pred_base_rate, pred_gam_rate, pred_rf_rate, pred_rfsp_rate),
               names_to = "model", values_to = "pred") %>%
  mutate(model = recode(model,
    pred_base_rate = "Persistence", pred_gam_rate = "GAM",
    pred_rf_rate   = "RF",          pred_rfsp_rate = "RFsp"
  )) %>%
  group_by(model) %>%
  summarise(
    MAE  = mae(obs_rate,  pred),
    RMSE = rmse(obs_rate, pred),
    Bias = bias(obs_rate, pred),
    R2   = r2(obs_rate,   pred),
    n    = sum(is.finite(obs_rate) & is.finite(pred)),
    .groups = "drop"
  )

cat("\n--- Static metrics (district weekly mean) ---\n"); print(metrics_weeklymean)
cat("\n--- Static metrics (site-week) ---\n");            print(metrics_siteweek)

write.csv(metrics_weeklymean, file.path(output_path, "metrics_static_weekly.csv"),   row.names = FALSE)
write.csv(metrics_siteweek,   file.path(output_path, "metrics_static_siteweek.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 8) Figure 2: Weekly mean time series (Observed + all models)
# -----------------------------------------------------------------------------
weekly_long_ts <- weekly_wide %>%
  pivot_longer(-week_start, names_to = "series", values_to = "rate") %>%
  mutate(
    series = recode(series,
      obs_rate       = "Observed",
      pred_base_rate = "Persistence",
      pred_gam_rate  = "GAM",
      pred_rf_rate   = "RF",
      pred_rfsp_rate = "RFsp"
    ),
    rate = pmax(rate, 0)
  )

# Peak week and magnitude per model (for annotation)
peak_long <- weekly_wide %>%
  pivot_longer(-week_start, names_to = "series", values_to = "rate") %>%
  mutate(series = recode(series,
    obs_rate = "Observed", pred_base_rate = "Persistence",
    pred_gam_rate = "GAM", pred_rf_rate = "RF", pred_rfsp_rate = "RFsp"
  )) %>%
  group_by(series) %>%
  slice_max(rate, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(peak_week = week_start, peak_mag = rate)

Figure2 <- ggplot(weekly_long_ts, aes(x = week_start, y = rate, color = series)) +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(data = peak_long,
             aes(x = peak_week, y = peak_mag, color = series),
             size = 3.6, show.legend = FALSE) +
  geom_text(data = peak_long,
            aes(x = peak_week, y = peak_mag, label = format(peak_week, "%b")),
            nudge_y = 5, size = 4.3, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = cb_pal,
                     breaks = c("Observed","Persistence","GAM","RF","RFsp")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
               expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "Weekly mean Aedes per trap-night", color = NULL) +
  theme_classic(base_size = 18) +
  theme(
    axis.title.y  = element_text(size = 17, face = "bold"),
    axis.text     = element_text(size = 17),
    legend.position = "right",
    legend.text   = element_text(size = 17)
  )

ggsave(
  file.path(fig_path, "Figure2_static_timeseries.pdf"),
  plot = Figure2, width = 10, height = 6, device = "pdf"
)

message("✓ 04_static_models.R complete — predictions, metrics, and Figure 2 saved.")
