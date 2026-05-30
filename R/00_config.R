# =============================================================================
# 00_config.R
# Global configuration: packages, paths, and shared helper functions
#
# SOURCE THIS SCRIPT FIRST before running any other script.
# All other scripts call source("R/00_config.R") at the top.
# =============================================================================

# -----------------------------------------------------------------------------
# 1) Packages
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(MMWRweek)
library(ggplot2)
library(patchwork)
library(scales)
library(sf)
library(sp)
library(gstat)
library(mgcv)
library(ranger)
library(FNN)
library(RSpectra)
library(spdep)
library(adespatial)
library(gratia)
library(tidytext)
library(fields)
library(here)
library(tibble)
library(beepr)       # optional: plays a sound when long jobs finish

sf::sf_use_s2(FALSE) # suppress s2 geometry warnings for planar operations

# -----------------------------------------------------------------------------
# 2) Paths  (all relative to project root via {here})
# -----------------------------------------------------------------------------
data_raw_path    <- here::here("data", "raw")
data_proc_path   <- here::here("data", "processed")
fig_path         <- here::here("figures")
output_path      <- here::here("output")
boundary_path    <- here::here("inst", "boundaries", "district")
counties_path    <- here::here("inst", "boundaries", "Cali_3310")

# Create output directories if they don't exist yet
dir.create(data_proc_path, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_path,       showWarnings = FALSE, recursive = TRUE)
dir.create(output_path,    showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 3) Shared helper functions
# -----------------------------------------------------------------------------

#' Project lon/lat columns to UTM Zone 11N (EPSG:32611) and add X, Y columns
#'
#' @param df   data frame with longitude and latitude columns
#' @param lon  name of longitude column (default "longitude")
#' @param lat  name of latitude column  (default "latitude")
#' @param epsg target CRS (default 32611 = UTM Zone 11N, meters)
#' @return df with added numeric columns X and Y in meters
to_xy <- function(df, lon = "longitude", lat = "latitude", epsg = 32611) {
  s  <- sf::st_as_sf(df, coords = c(lon, lat), crs = 4326) |>
    sf::st_transform(epsg)
  xy <- sf::st_coordinates(s)
  df$X <- xy[, 1]
  df$Y <- xy[, 2]
  df
}

#' Effort-weighted mean (used for district-week covariate aggregation)
#'
#' @param x numeric vector of covariate values
#' @param w numeric vector of weights (e.g., trap_nights)
#' @return scalar weighted mean, or NA if no finite pairs
w_mean <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  k <- is.finite(x) & is.finite(w)
  if (!any(k)) return(NA_real_)
  sum(x[k] * w[k]) / sum(w[k])
}

#' Mean absolute error
mae <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)

#' Root mean squared error
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

#' R-squared (coefficient of determination)
r2 <- function(y, yhat) {
  ok   <- is.finite(y) & is.finite(yhat)
  y    <- y[ok]; yhat <- yhat[ok]
  ybar <- mean(y)
  1 - sum((y - yhat)^2) / sum((y - ybar)^2)
}

#' Mean bias (predicted - observed)
bias <- function(y, yhat) mean(yhat - y, na.rm = TRUE)

#' Relative MAE (percentage, with small constant c to avoid division by zero)
rel_mae <- function(y, yhat, c = 1) {
  mean(abs(y - yhat) / (y + c) * 100, na.rm = TRUE)
}

#' Relative RMSE (percentage)
rel_rmse <- function(y, yhat, c = 1) {
  sqrt(mean(((yhat - y) / (y + c) * 100)^2, na.rm = TRUE))
}

#' Build kernel PCA spatial features for RFsp (Nyström projection)
#'
#' Computes Gaussian kernel PCA on training coordinates, optionally filters
#' principal components by Moran's I significance, then projects test
#' coordinates via Nyström approximation.
#'
#' @param train       training data frame (must contain x_col, y_col)
#' @param test        test/prediction data frame (must contain x_col, y_col)
#' @param x_col       name of easting column  (default "X", meters)
#' @param y_col       name of northing column (default "Y", meters)
#' @param bandwidth   Gaussian kernel bandwidth in meters (default 500)
#' @param q           number of PCs to compute (default 10)
#' @param k_nn        number of neighbors for Moran graph (default 8)
#' @param moran_filter if TRUE, retain only PCs with significant spatial
#'                    autocorrelation (Moran's I p < 0.05)
#' @return list with elements:
#'   train      (train with PC columns appended),
#'   test       (test  with PC columns appended),
#'   kept_pcs   (character vector of retained PC names),
#'   eigenvalues (numeric vector of eigenvalues),
#'   bandwidth  (bandwidth used)
build_rfsp_features <- function(train, test,
                                x_col = "X", y_col = "Y",
                                bandwidth = 500,
                                q = 10,
                                k_nn = 8,
                                moran_filter = TRUE) {
  stopifnot(all(c(x_col, y_col) %in% names(train)),
            all(c(x_col, y_col) %in% names(test)))

  A <- as.matrix(train[, c(x_col, y_col)])
  B <- as.matrix(test[,  c(x_col, y_col)])

  kfun <- function(U, V, h) exp(-fields::rdist(U, V)^2 / (h^2))

  # Train kernel + double-centering
  K  <- kfun(A, A, bandwidth); K <- (K + t(K)) / 2
  cm <- colMeans(K); rm <- rowMeans(K); m <- mean(K)
  Kc <- K - matrix(cm, nrow(K), ncol(K), byrow = TRUE) -
    matrix(rm, nrow(K), ncol(K), byrow = FALSE) + m
  Kc <- (Kc + t(Kc)) / 2

  # Eigenpairs
  q_eff <- max(1, min(q, nrow(Kc) - 1L))
  eg    <- RSpectra::eigs_sym(Kc, k = q_eff, which = "LA")
  lam   <- pmax(Re(eg$values), 0)
  V_mat <- Re(eg$vectors)
  inv_sqrtLam <- diag(1 / sqrt(pmax(lam, .Machine$double.eps)), q_eff, q_eff)

  # Train scores
  Z_train <- Kc %*% V_mat %*% inv_sqrtLam
  colnames(Z_train) <- paste0("PC", seq_len(q_eff))

  # Nyström projection for test
  K_nm  <- kfun(B, A, bandwidth)
  rm_n  <- rowMeans(K_nm)
  Kc_nm <- K_nm - matrix(cm, nrow(K_nm), ncol(K_nm), byrow = TRUE) -
    matrix(rm_n, nrow(K_nm), ncol(K_nm), byrow = FALSE) + m
  Z_test <- Kc_nm %*% V_mat %*% inv_sqrtLam
  colnames(Z_test) <- paste0("PC", seq_len(q_eff))

  # Optional Moran's I filtering
  keep <- colnames(Z_train)
  if (moran_filter) {
    nn  <- FNN::get.knn(A, k = k_nn)$nn.index
    n   <- nrow(A); W <- matrix(0, n, n)
    for (i in seq_len(n)) W[i, nn[i, ]] <- 1
    rs <- rowSums(W); rs[rs == 0] <- 1
    W  <- W / rs
    lw <- spdep::mat2listw(W, style = "W")
    pv <- sapply(seq_len(q_eff), function(j)
      adespatial::moran.randtest(Z_train[, j], lw, nrepet = 499)$pvalue
    )
    keep <- colnames(Z_train)[pv < 0.05]
    if (!length(keep)) keep <- colnames(Z_train)  # fallback: keep all
  }

  train_sp <- dplyr::bind_cols(train, as.data.frame(Z_train[, keep, drop = FALSE]))
  test_sp  <- dplyr::bind_cols(test,  as.data.frame(Z_test[,  keep, drop = FALSE]))

  list(
    train       = train_sp,
    test        = test_sp,
    kept_pcs    = keep,
    eigenvalues = lam,
    bandwidth   = bandwidth
  )
}

# -----------------------------------------------------------------------------
# 4) Colorblind-safe palette (Okabe-Ito) used across all figures
# -----------------------------------------------------------------------------
cb_pal <- c(
  "Observed"    = "black",
  "Persistence" = "#D55E00",   # orange-red
  "Baseline"    = "#D55E00",
  "GAM"         = "#009E73",   # bluish green
  "RF"          = "#0072B2",   # blue
  "RFsp"        = "#CC79A7"    # reddish purple
)

# -----------------------------------------------------------------------------
# 5) Save session info for reproducibility record
# -----------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           here::here("session_info.txt"))

message("✓ 00_config.R loaded — paths set, helpers defined.")
