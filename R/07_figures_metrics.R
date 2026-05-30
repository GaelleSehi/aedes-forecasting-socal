# =============================================================================
# 07_figures_metrics.R
# All manuscript figures from rolling forecast results + exported metric tables
#
# Requires: data/processed/roll_results.rds
#           data/processed/metrics_rolling.rds
#
# Produces:
#   Figure 3  — Weekly MAE + RMSE (smoothed), rolling re-estimation
#   Figure 4  — 2-week-ahead forecasts with 95% CI bands per model
#   Figure 5  — Predicted vs observed scatter plots (Pearson r)
#   Figure 6  — Trap-level mean error maps (July–October 2023)  [see 06_spatial_maps.R]
#   FigureS1  — Weekly relative MAE + RMSE (supplementary)
#   FigureS2  — Relative error maps (supplementary)             [see 06_spatial_maps.R]
#   output/metrics_rolling.csv  (already written by 05_rolling_forecast.R)
# =============================================================================

source(here::here("R", "00_config.R"))

results         <- readRDS(file.path(data_proc_path, "roll_results.rds"))
metrics_overall <- readRDS(file.path(data_proc_path, "metrics_rolling.rds"))

# Back-transform to rate scale if not already present
results <- results %>%
  mutate(
    obs_rate_t2  = pmax(expm1(lograte_lead2), 0),
    base_rate_t2 = pmax(expm1(y_lograte_t),   0)
  )

mask       <- with(results,
                   is.finite(obs_rate_t2) & is.finite(base_rate_t2) &
                     is.finite(pred_gam_rate) & is.finite(pred_rf_rate) &
                     is.finite(pred_rfsp_rate))
results_aa <- results[mask, ]

# =============================================================================
# FIGURE 3: Weekly MAE + RMSE (smoothed) — rolling re-estimation
# =============================================================================
week_var <- if ("cut_date" %in% names(results_aa)) "cut_date" else "week_start"

weekly_long <- results_aa %>%
  transmute(
    week = .data[[week_var]],
    obs  = obs_rate_t2,
    base = base_rate_t2,
    gam  = pred_gam_rate,
    rf   = pred_rf_rate,
    rfsp = pred_rfsp_rate
  ) %>%
  pivot_longer(c(base, gam, rf, rfsp), names_to = "model", values_to = "pred") %>%
  filter(is.finite(obs), is.finite(pred)) %>%
  mutate(model = factor(model,
                        levels = c("base","gam","rf","rfsp"),
                        labels = c("Persistence","GAM","RF","RFsp")))

metrics_by_week <- weekly_long %>%
  group_by(week, model) %>%
  summarise(
    MAE  = mae(obs, pred),
    RMSE = rmse(obs, pred),
    n    = dplyr::n(),
    .groups = "drop"
  )

Figure3 <- metrics_by_week %>%
  pivot_longer(c(MAE, RMSE), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = week, y = Value, color = model, group = model)) +
  geom_line(alpha = 0.35, linewidth = 0.8) +
  geom_smooth(se = FALSE, linewidth = 1.2) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = cb_pal[c("Persistence","GAM","RF","RFsp")]) +
  labs(x = NULL, y = "Error", color = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text    = element_text(size = 15, face = "bold"),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 14),
    legend.title  = element_text(size = 15, face = "bold"),
    legend.text   = element_text(size = 14)
  )

ggsave(file.path(fig_path, "Figure3_weekly_error.pdf"),
       plot = Figure3, width = 10, height = 8, device = "pdf")

# =============================================================================
# FIGURE 4: 2-week-ahead forecasts with 95% CI (mean ± 1.96 SE across sites)
# =============================================================================
long_h2 <- results_aa %>%
  transmute(
    week = as.Date(dplyr::coalesce(cut_date, week_start)),
    obs  = obs_rate_t2,
    base = base_rate_t2,
    gam  = pred_gam_rate,
    rf   = pred_rf_rate,
    rfsp = pred_rfsp_rate
  ) %>%
  pivot_longer(c(base, gam, rf, rfsp), names_to = "model", values_to = "pred") %>%
  filter(!is.na(week), is.finite(obs), is.finite(pred)) %>%
  mutate(model = factor(model,
                        levels = c("base","gam","rf","rfsp"),
                        labels = c("A. Persistence","B. GAM","C. RF","D. RFsp")))

summ_h2 <- long_h2 %>%
  group_by(model, week) %>%
  summarise(
    mean_pred = mean(pred, na.rm = TRUE),
    se_pred   = sd(pred, na.rm = TRUE) / pmax(sqrt(dplyr::n()), 1),
    lo        = pmax(mean_pred - 1.96 * se_pred, 0),
    hi        = mean_pred + 1.96 * se_pred,
    mean_obs  = mean(obs, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(week = as.Date(week))

pal_h2 <- c(
  "A. Persistence" = "#D55E00",
  "B. GAM"         = "#009E73",
  "C. RF"          = "#0072B2",
  "D. RFsp"        = "#CC79A7"
)

Figure4 <- ggplot(summ_h2, aes(x = week)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.16) +
  geom_line(aes(y = mean_obs), color = "black", linewidth = 1.0, lineend = "round") +
  geom_line(aes(y = mean_pred, color = model), linewidth = 1.25, lineend = "round") +
  facet_wrap(~model, ncol = 2, scales = "fixed") +
  scale_color_manual(values = pal_h2, guide = "none") +
  scale_fill_manual(values = pal_h2, guide = "none") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
               expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Weekly mean Aedes per trap-night",
       caption = "Black line: observed; shaded bands: mean ± 1.96 SE across sites.") +
  theme_minimal(base_size = 18) +
  theme(
    strip.text       = element_text(size = 16, face = "bold"),
    axis.title.y     = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_path, "Figure4_rolling_timeseries_CI.pdf"),
       plot = Figure4, width = 10, height = 7, device = "pdf")

# =============================================================================
# FIGURE 5: Predicted vs observed scatter plots (Pearson r)
# =============================================================================
plot_df <- results_aa %>%
  transmute(
    week = dplyr::coalesce(cut_date, week_start),
    obs  = obs_rate_t2,
    base = base_rate_t2,
    gam  = pred_gam_rate,
    rf   = pred_rf_rate,
    rfsp = pred_rfsp_rate
  ) %>%
  pivot_longer(c(base, gam, rf, rfsp), names_to = "model", values_to = "pred") %>%
  filter(is.finite(obs), is.finite(pred)) %>%
  mutate(model = factor(model,
                        levels = c("base","gam","rf","rfsp"),
                        labels = c("A. Baseline","B. GAM","C. RF","D. RFsp")))

r_df <- plot_df %>%
  group_by(model) %>%
  summarise(
    r = suppressWarnings(cor(obs, pred, method = "pearson", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f", r))

lim <- quantile(c(plot_df$obs, plot_df$pred), 0.99, na.rm = TRUE)
r_df <- r_df %>% mutate(x = 0.06 * lim, y = 0.92 * lim)

Figure5 <- ggplot(plot_df, aes(obs, pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.35, size = 1) +
  geom_text(data = r_df, aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE) +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
  facet_wrap(~model, ncol = 2) +
  scale_x_continuous("Observed rate (site-level, t+2)") +
  scale_y_continuous("Predicted rate (site-level, t+2)") +
  labs(title = "") +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(fig_path, "Figure5_scatter_obs_pred.pdf"),
       plot = Figure5, width = 10, height = 8, device = "pdf")

# =============================================================================
# SUPPLEMENTARY FIGURE S1: Weekly relative MAE + RMSE
# =============================================================================
weekly_long_rel <- results_aa %>%
  transmute(
    week = .data[[week_var]],
    obs  = obs_rate_t2,
    base = base_rate_t2,
    gam  = pred_gam_rate,
    rf   = pred_rf_rate,
    rfsp = pred_rfsp_rate
  ) %>%
  pivot_longer(c(base, gam, rf, rfsp), names_to = "model", values_to = "pred") %>%
  filter(is.finite(obs), is.finite(pred)) %>%
  mutate(model = factor(model,
                        levels = c("base","gam","rf","rfsp"),
                        labels = c("Persistence","GAM","RF","RFsp")))

rel_metrics_by_week <- weekly_long_rel %>%
  group_by(week, model) %>%
  summarise(
    Relative_MAE  = rel_mae(obs, pred, c = 1),
    Relative_RMSE = rel_rmse(obs, pred, c = 1),
    .groups = "drop"
  )

FigureS1 <- rel_metrics_by_week %>%
  pivot_longer(c(Relative_MAE, Relative_RMSE), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
    Relative_MAE  = "Relative MAE",
    Relative_RMSE = "Relative RMSE"
  )) %>%
  ggplot(aes(week, Value, color = model)) +
  geom_line(alpha = 0.35) +
  geom_smooth(se = FALSE, linewidth = 0.9) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = cb_pal[c("Persistence","GAM","RF","RFsp")]) +
  labs(x = NULL, y = "Relative error (%)", color = NULL) +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_path, "FigureS1_relative_error.pdf"),
       plot = FigureS1, width = 8, height = 6, device = "pdf")

# =============================================================================
# Print final metric table
# =============================================================================
cat("\n=== Rolling forecast performance (2-week ahead, 2023) ===\n")
print(metrics_overall)

message("✓ 07_figures_metrics.R complete — Figures 3–5 and S1 saved.")
