# =============================================================================
# 03_eda_figures.R
# Exploratory data analysis and Figure 1 (operational context)
#
# Produces:
#   Figure 1A — spatial map of mean site-level catch rate (2017–2022)
#   Figure 1B — inter-annual peak catch rate by year
#   Figure 1C — district-wide weekly time series with 2023 highlighted
#   figures/Figure1_operational_context.pdf
# =============================================================================

source(here::here("R", "00_config.R"))

site_week     <- readRDS(file.path(data_proc_path, "site_week_full.rds"))
district_week <- readRDS(file.path(data_proc_path, "district_week_sc.rds"))

# Load study-area boundary shapefile
# NOTE: place your boundary shapefile in inst/boundaries/district/
boundary_raw <- sf::st_read(boundary_path, quiet = TRUE)
boundary_ll  <- boundary_raw %>%
  sf::st_make_valid() %>%
  sf::st_transform(4326)

# -----------------------------------------------------------------------------
# 1) Site-level hotspot summary (training period only: 2017–2022)
# Cap at 99th percentile for visualization — does not affect analysis.
# -----------------------------------------------------------------------------
site_hotspot <- site_week %>%
  filter(week_year <= 2022) %>%
  group_by(site_code) %>%
  summarise(
    longitude   = median(longitude, na.rm = TRUE),
    latitude    = median(latitude,  na.rm = TRUE),
    total_cnt   = sum(Total_Aedes_week, na.rm = TRUE),
    total_eff   = sum(trap_nights,      na.rm = TRUE),
    mean_rate_w = total_cnt / pmax(total_eff, 1),
    .groups = "drop"
  )

cap_p99 <- quantile(site_hotspot$mean_rate_w, 0.99, na.rm = TRUE)
site_hotspot <- site_hotspot %>%
  mutate(
    rate_plot  = pmin(mean_rate_w, cap_p99),
    alpha_plot = pmax(rate_plot / cap_p99, 0.05)
  )

# -----------------------------------------------------------------------------
# 2) Panel A: Spatial map of mean catch rate
# -----------------------------------------------------------------------------
p_map <- ggplot() +
  geom_sf(data = boundary_ll, fill = NA, linewidth = 0.4, color = "grey40") +
  geom_point(
    data = site_hotspot,
    aes(x = longitude, y = latitude, color = rate_plot, alpha = alpha_plot),
    size = 2
  ) +
  coord_sf(crs = 4326) +
  scale_alpha(range = c(0.10, 1), guide = "none") +
  scale_color_gradient(
    low  = "grey90",
    high = "darkred",
    name = "Mean rate\n(2017–2022)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x        = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 8),
    legend.position    = c(1.1, 0.30),
    legend.justification = c(1, 0.5),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.margin      = margin(2, 2, 2, 2)
  )

# -----------------------------------------------------------------------------
# 3) Panel B: Inter-annual peak catch rate
# -----------------------------------------------------------------------------
peaks_by_year <- district_week %>%
  group_by(week_year) %>%
  summarise(
    peak_week = week_start[which.max(rate_week)],
    peak_rate = max(rate_week, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(peak_label = format(peak_week, "%b %d"))

p_peak <- ggplot(peaks_by_year, aes(x = week_year, y = peak_rate)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  geom_text(aes(label = peak_label), vjust = -0.7, size = 3) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(x = "Year", y = "Peak mosquitoes per trap-night") +
  theme_minimal(base_size = 12)

# -----------------------------------------------------------------------------
# 4) Panel C: District-wide weekly time series (2017–2023)
# 2023 is highlighted in a heavier line; peak week marked with dashed line.
# -----------------------------------------------------------------------------
peak_week_2023 <- district_week %>%
  filter(week_year == 2023) %>%
  slice_max(rate_week, n = 1, with_ties = FALSE) %>%
  pull(week_start)

p_ts <- ggplot(district_week, aes(week_start, rate_week)) +
  geom_line(alpha = 0.5) +
  geom_line(data = subset(district_week, week_year == 2023), linewidth = 0.9) +
  geom_vline(xintercept = peak_week_2023, linetype = 2) +
  labs(x = "Week start", y = "Mean Aedes per trap-night") +
  coord_cartesian(ylim = c(0, 65)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_minimal(base_size = 16)

# -----------------------------------------------------------------------------
# 5) Overlay: all years by week of year (grey gradient, 2023 in black)
# -----------------------------------------------------------------------------
gray_vals <- c(
  "2017" = "gray90", "2018" = "gray80", "2019" = "gray70",
  "2020" = "gray60", "2021" = "gray50", "2022" = "gray40", "2023" = "gray0"
)

district_week_plot <- district_week %>%
  filter(week_year >= 2017)

p_overlay_gray <- ggplot(
  district_week_plot,
  aes(x = week_num, y = rate_week,
      group = factor(week_year), color = factor(week_year))
) +
  geom_line(aes(linewidth = factor(week_year) == "2023"), alpha = 0.80) +
  scale_color_manual(values = gray_vals) +
  scale_linewidth_manual(values = c(`TRUE` = 1.4, `FALSE` = 0.8), guide = "none") +
  labs(x = "Week of year", y = "Mean Aedes per trap-night", color = "Year") +
  theme_minimal(base_size = 16)

# -----------------------------------------------------------------------------
# 6) Assemble Figure 1 (three-panel composite)
# -----------------------------------------------------------------------------
top_row    <- (p_map  + labs(tag = "A")) | (p_peak + labs(tag = "B"))
Figure1    <- (top_row / (p_ts + labs(tag = "C"))) +
  plot_layout(widths = c(1.1, 2.2), heights = c(1.2, 0.8)) &
  theme(
    plot.tag          = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1)
  )

ggsave(
  file.path(fig_path, "Figure1_operational_context.pdf"),
  plot   = Figure1,
  device = "pdf",
  width  = 8.5, height = 4.6, units = "in",
  dpi    = 300
)

message("✓ 03_eda_figures.R complete — Figure 1 saved.")
