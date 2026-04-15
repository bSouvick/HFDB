#!/usr/bin/env Rscript

suppressMessages({
  library(fields)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(doRNG)
})

##
## Command line arguments
##
args <- commandArgs(trailingOnly = TRUE)

mode <- if (length(args) >= 1) args[1] else "all"                   # precompute / subsample / all
betaValue <- if (length(args) >= 2) as.numeric(args[2]) else 7      # block length
cache_dir <- if (length(args) >= 3) args[3] else "cache_variogram_model_ng"

if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cat("Mode      :", mode, "\n")
cat("betaValue :", betaValue, "\n")
cat("cache_dir :", cache_dir, "\n")

start_time <- proc.time()

##
## Utility: cores for HPC
##
get_ncores <- function(max_cores = 6L) {
  slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA))
  if (!is.na(slurm_cores) && slurm_cores >= 1) {
    return(slurm_cores)
  }
  n_cores <- min(max_cores, parallel::detectCores(logical = TRUE))
  if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
  n_cores
}

##
## Global parameters
##
set.seed(2)

nMonteCarlo <- 1000
nBoot       <- 500

Delta  <- 1
xSeq   <- seq(1, 50, by = Delta)
ySeq   <- seq(1, 50, by = Delta)

n1 <- length(xSeq)
n2 <- length(ySeq)
N  <- n1 * n2

grid   <- list(x = xSeq, y = ySeq)
Locs   <- expand.grid(x = xSeq, y = ySeq, KEEP.OUT.ATTRS = FALSE)
ProInd <- as.matrix(Locs)

lag_mat <- rbind(c(1, 0),
                 c(1, 1),
                 c(0, 1))
phi_working <- 2.5

#---smoothing aprameters---
hx <- 0.05
hy <- 0.05
kernel_type <- "gaussian"

## cache file names
lag_tag <- paste(apply(lag_mat, 1, paste, collapse = "_"), collapse = "__")
base_tag <- paste0(
  "variogram_model_ng_n1_", n1,
  "_n2_", n2,
  "_MC_", nMonteCarlo,
  "_B_", nBoot,
  "_phiw_", format(phi_working, trim = TRUE, scientific = FALSE),
  "_lags_", lag_tag
)

full_cache_file <- file.path(cache_dir, paste0("full_precompute_", base_tag, ".rds"))
sub_cache_file  <- function(beta) file.path(cache_dir, paste0("subsample_beta_", beta, "_", base_tag, ".rds"))
res_cache_file  <- function(beta) file.path(cache_dir, paste0("final_results_beta_", beta, "_", base_tag, ".rds"))

##
## Non-Gaussian field
##
rinnov <- function(n) {
  x <- rgamma(n, shape = 0.2, scale = 1)
  (x - 0.2) / sqrt(0.2)
}

A <- matrix(c(
  1.00, 0.30, 0.10,
  0.50, 0.20, 0.10,
  0.25, 0.10, 0.05
), nrow = 3, byrow = TRUE)

A <- A / sqrt(sum(A^2))

simulate_ma22_field <- function(n1, n2, A, innov_fun) {
  eps <- matrix(innov_fun((n1 + 2) * (n2 + 2)), nrow = n1 + 2, ncol = n2 + 2)
  Z <- matrix(0, nrow = n1, ncol = n2)

  for (i in 1:n1) {
    for (j in 1:n2) {
      Z[i, j] <-
        A[1,1] * eps[i + 2, j + 2] +
        A[2,1] * eps[i + 1, j + 2] +
        A[3,1] * eps[i,     j + 2] +
        A[1,2] * eps[i + 2, j + 1] +
        A[2,2] * eps[i + 1, j + 1] +
        A[3,2] * eps[i,     j + 1] +
        A[1,3] * eps[i + 2, j    ] +
        A[2,3] * eps[i + 1, j    ] +
        A[3,3] * eps[i,     j    ]
    }
  }

  Z
}

##
## Working variogram model
##
g0_fun <- function(h, phi = phi_working) {
  hnorm <- sqrt(sum(h^2))
  1 - exp(-hnorm / phi)
}

g0_vec <- apply(lag_mat, 1, g0_fun, phi = phi_working)

cov_ma22_lag <- function(h, A) {
  hx <- h[1]
  hy <- h[2]
  out <- 0

  for (u in 0:2) {
    for (v in 0:2) {
      u2 <- u + hx
      v2 <- v + hy
      if (u2 >= 0 && u2 <= 2 && v2 >= 0 && v2 <= 2) {
        out <- out + A[u + 1, v + 1] * A[u2 + 1, v2 + 1]
      }
    }
  }

  out
}

gamma_true_lag <- function(h, A) {
  1 - cov_ma22_lag(h, A)
}

gamma_true_vec <- apply(lag_mat, 1, gamma_true_lag, A = A)
theta_true <- sum(gamma_true_vec * g0_vec) / sum(g0_vec^2)

##
## Periodogram helpers
##
periodogram_general <- function(Z, locs, omega) {
  Z <- as.numeric(Z)
  locs <- as.matrix(locs)
  omega <- as.numeric(omega)
  n <- nrow(locs)

  phase <- locs %*% omega
  dft <- sum(Z * exp(-1i * phase))
  (2 * pi)^(-2) * (1 / n) * Mod(dft)^2
}

periodogramSub <- function(Z, locs, omega) {
  Z <- as.numeric(Z)
  locs <- as.matrix(locs)
  omega <- as.numeric(omega)
  n <- nrow(locs)

  phase <- locs %*% omega
  dft <- sum(Z * exp(-1i * phase))
  (2 * pi)^(-2) * (1 / n) * Mod(dft)^2
}

##
## Frequency-grid helpers
##
centered_range <- function(n) {
  if (n %% 2 == 0) {
    seq(-n / 2, n / 2 - 1, by = 1)
  } else {
    seq(-(n - 1) / 2, (n - 1) / 2, by = 1)
  }
}

reflect_centered <- function(j, n) {
  if (n %% 2 == 0) {
    ifelse(j == -n / 2, -n / 2, -j)
  } else {
    -j
  }
}

build_freq_objects <- function(n1, n2) {
  j1_range <- centered_range(n1)
  j2_range <- centered_range(n2)

  f1 <- (2 * pi * j1_range) / n1
  f2 <- (2 * pi * j2_range) / n2

  freqGrid_full <- expand.grid(f1 = f1, f2 = f2, KEEP.OUT.ATTRS = FALSE)
  freqGrid_full <- as.matrix(freqGrid_full)

  zero_logical <- (freqGrid_full[, 1] == 0 & freqGrid_full[, 2] == 0)
  freqGrid <- freqGrid_full[!zero_logical, , drop = FALSE]

  Jn <- expand.grid(j1 = j1_range, j2 = j2_range, KEEP.OUT.ATTRS = FALSE)
  Jn$row <- match(Jn$j2, rev(j2_range))
  Jn$col <- match(Jn$j1, j1_range)

  zero_ind <- which(Jn$j1 == 0 & Jn$j2 == 0)

  partner_j1 <- reflect_centered(Jn$j1, n1)
  partner_j2 <- reflect_centered(Jn$j2, n2)

  partner_ind <- match(
    paste(partner_j1, partner_j2),
    paste(Jn$j1, Jn$j2)
  )

  if (anyNA(partner_ind)) stop("Some frequencies do not have symmetric partners.")

  self_ind  <- which(partner_ind == seq_len(nrow(Jn)))
  pos_ind   <- which(seq_len(nrow(Jn)) < partner_ind)
  neg_ind   <- partner_ind[pos_ind]
  resid_ind <- setdiff(self_ind, zero_ind)

  flat_order <- order(Jn$row, Jn$col)
  use_ind <- setdiff(flat_order, zero_ind)

  list(
    f1 = f1,
    f2 = f2,
    freqGrid = freqGrid,
    Jn = Jn,
    zero_ind = zero_ind,
    pos_ind = pos_ind,
    neg_ind = neg_ind,
    resid_ind = resid_ind,
    use_ind = use_ind
  )
}

freq_obj_full <- build_freq_objects(n1, n2)

##
## Kernel smoothing helpers
##
kernel_2d <- function(u1, u2, type = c("epanechnikov", "gaussian")) {
  type <- match.arg(type)

  if (type == "epanechnikov") {
    out <- (9 / 16) * (1 - u1^2) * (1 - u2^2)
    out[abs(u1) > 1 | abs(u2) > 1] <- 0
    out[out < 0] <- 0
    return(out)
  }

  (1 / (2 * pi)) * exp(-0.5 * (u1^2 + u2^2))
}

wrap_diff_scalar <- function(a, b) {
  d <- a - b
  ((d + pi) %% (2 * pi)) - pi
}

build_literal_kernel_weights <- function(f1, f2, hx, hy,
                                         kernel = c("epanechnikov", "gaussian"),
                                         periodic = TRUE) {
  kernel <- match.arg(kernel)

  col_freq <- as.numeric(f1)
  row_freq <- as.numeric(rev(f2))

  nr <- length(row_freq)
  nc <- length(col_freq)

  W <- array(0, dim = c(nr, nc, nr, nc))
  D <- matrix(0, nrow = nr, ncol = nc)

  for (r0 in seq_len(nr)) {
    lambda2 <- row_freq[r0]
    for (c0 in seq_len(nc)) {
      lambda1 <- col_freq[c0]
      den <- 0
      for (r in seq_len(nr)) {
        lamj2 <- row_freq[r]
        d2 <- if (periodic) wrap_diff_scalar(lambda2, lamj2) else (lambda2 - lamj2)
        for (c in seq_len(nc)) {
          lamj1 <- col_freq[c]
          d1 <- if (periodic) wrap_diff_scalar(lambda1, lamj1) else (lambda1 - lamj1)
          w <- kernel_2d(d1 / hx, d2 / hy, type = kernel)
          W[r0, c0, r, c] <- w
          den <- den + w
        }
      }
      D[r0, c0] <- den
    }
  }

  list(W = W, D = D)
}

smooth_periodogram_fast <- function(pMatrix, weight_obj, zero_row = NULL, zero_col = NULL) {
  W <- weight_obj$W
  D <- weight_obj$D

  nr <- dim(W)[1]
  nc <- dim(W)[2]

  if (!all(dim(pMatrix) == c(nr, nc))) {
    stop("pMatrix dimensions do not match weight object.")
  }

  out <- matrix(0, nrow = nr, ncol = nc)

  for (r0 in seq_len(nr)) {
    for (c0 in seq_len(nc)) {
      num <- sum(W[r0, c0, , ] * pMatrix)
      den <- D[r0, c0]
      out[r0, c0] <- if (den > 0) num / den else NA_real_
    }
  }

  if (!is.null(zero_row) && !is.null(zero_col)) {
    out[zero_row, zero_col] <- 0
  }

  out
}

weight_obj <- build_literal_kernel_weights(
  f1 = freq_obj_full$f1,
  f2 = freq_obj_full$f2,
  hx = hx,
  hy = hy,
  kernel = kernel_type,
  periodic = TRUE
)

##
## Estimating-equation objects
##
estimate_gamma_spec <- function(Zvec, locs, freqGrid, lag_mat) {
  pdgramVals <- apply(freqGrid, 1, function(omega) periodogram_general(Zvec, locs, omega))
  psi_mat <- sapply(seq_len(nrow(lag_mat)), function(k) {
    1 - cos(freqGrid %*% lag_mat[k, ])
  })
  gamma_hat <- as.numeric(((2 * pi)^2 / nrow(locs)) * colSums(psi_mat * pdgramVals))
  list(gamma_hat = gamma_hat, pdgramVals = pdgramVals, psi_mat = psi_mat)
}

estimate_theta_from_gamma <- function(gamma_hat_vec, g0_vec) {
  sum(gamma_hat_vec * g0_vec) / sum(g0_vec^2)
}

c_theta <- g0_vec / sum(g0_vec^2)

build_phi_theta <- function(freqGrid, lag_mat, c_theta) {
  psi_mat <- sapply(seq_len(nrow(lag_mat)), function(k) {
    1 - cos(freqGrid %*% lag_mat[k, ])
  })
  as.numeric(psi_mat %*% c_theta)
}

phi_theta_full <- build_phi_theta(freq_obj_full$freqGrid, lag_mat, c_theta)
phi_pair_full  <- phi_theta_full * (phi_theta_full + phi_theta_full)

##
## Full precompute stage
##
precompute_full <- function() {
  cat("Starting full precomputation for variogram model fitting case...\n")

  dataset <- matrix(NA_real_, nrow = N, ncol = nMonteCarlo)
  for (i in seq_len(nMonteCarlo)) {
    dataset[, i] <- c(simulate_ma22_field(n1, n2, A, rinnov))
  }

  n_cores <- get_ncores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  gEstimate <- foreach(
    q = seq_len(nMonteCarlo),
    .combine = c,
    .export = c("dataset", "ProInd", "freq_obj_full", "lag_mat",
                "estimate_gamma_spec", "g0_vec", "estimate_theta_from_gamma",
                "periodogram_general")
  ) %dopar% {
    z <- dataset[, q]
    full_est <- estimate_gamma_spec(z, ProInd, freq_obj_full$freqGrid, lag_mat)
    estimate_theta_from_gamma(full_est$gamma_hat, g0_vec)
  }

  stopCluster(cl)

  mDataVar <- (N / nMonteCarlo) * sum((gEstimate - mean(gEstimate))^2)

  n_cores <- get_ncores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  doRNG::registerDoRNG(123)

  v_hpb    <- matrix(NA_real_, nrow = nBoot, ncol = nMonteCarlo)
  var_vHpb <- numeric(nMonteCarlo)
  boot_ctr <- numeric(nMonteCarlo)

  for (q in seq_len(nMonteCarlo)) {
    z <- dataset[, q]
    full_est <- estimate_gamma_spec(z, ProInd, freq_obj_full$freqGrid, lag_mat)
    pdgramVals_full <- full_est$pdgramVals

    pMatrix <- matrix(0, nrow = length(freq_obj_full$f2), ncol = length(freq_obj_full$f1))
    Jn <- freq_obj_full$Jn
    use_ind <- freq_obj_full$use_ind
    zero_ind <- freq_obj_full$zero_ind

    pMatrix[cbind(Jn$row[use_ind], Jn$col[use_ind])] <- pdgramVals_full
    zr <- Jn$row[zero_ind]
    zc <- Jn$col[zero_ind]
    pMatrix[zr, zc] <- 0

    mat <- smooth_periodogram_fast(
      pMatrix = pMatrix,
      weight_obj = weight_obj,
      zero_row = zr,
      zero_col = zc
    )

    smoothPdgram <- mat[cbind(Jn$row[use_ind], Jn$col[use_ind])]

    v_hpb_boot <- foreach(
      l = seq_len(nBoot),
      .combine = c,
      .export = c("phi_theta_full", "smoothPdgram", "freq_obj_full", "N")
    ) %dorng% {
      expVals <- numeric(nrow(freq_obj_full$Jn))

      epos <- rexp(length(freq_obj_full$pos_ind), rate = 1) - 1
      expVals[freq_obj_full$pos_ind] <- epos
      expVals[freq_obj_full$neg_ind] <- epos

      if (length(freq_obj_full$resid_ind) > 0) {
        expVals[freq_obj_full$resid_ind] <- rexp(length(freq_obj_full$resid_ind), rate = 1) - 1
      }

      expVector <- expVals[freq_obj_full$use_ind]
      vPart <- phi_theta_full * smoothPdgram * expVector
      sqrt(N) * ((2 * pi)^2 / N) * sum(vPart)
    }

    v_hpb[, q] <- v_hpb_boot
    boot_ctr[q] <- mean(v_hpb_boot)
    var_vHpb[q] <- ((2 * pi)^4) * (1 / N) * sum(phi_pair_full * smoothPdgram^2)
  }

  stopCluster(cl)

  saveRDS(
    list(
      dataset = dataset,
      gEstimate = gEstimate,
      mDataVar = mDataVar,
      v_hpb = v_hpb,
      var_vHpb = var_vHpb,
      boot_ctr = boot_ctr,
      theta_true = theta_true,
      gamma_true_vec = gamma_true_vec,
      g0_vec = g0_vec,
      params = list(
        nMonteCarlo = nMonteCarlo,
        nBoot = nBoot,
        n1 = n1, n2 = n2, N = N,
        Delta = Delta,
        lag_mat = lag_mat,
        phi_working = phi_working,
        A = A,
        hx = hx,
        hy = hy,
        kernel_type = kernel_type
      )
    ),
    file = full_cache_file,
    compress = FALSE
  )

  cat("Saved full precompute cache to:\n", full_cache_file, "\n")
}

##
## Subsampling stage for one beta
##
run_subsample_for_beta <- function(betaValue) {
  if (!file.exists(full_cache_file)) {
    stop("Full cache file not found. Run mode='precompute' first.")
  }

  cat("Loading full cache...\n")
  full_obj <- readRDS(full_cache_file)

  dataset   <- full_obj$dataset
  gEstimate <- full_obj$gEstimate
  mDataVar  <- full_obj$mDataVar
  v_hpb     <- full_obj$v_hpb
  var_vHpb  <- full_obj$var_vHpb
  boot_ctr  <- full_obj$boot_ctr
  theta_true <- full_obj$theta_true

  nBlockx <- betaValue
  nBlocky <- betaValue
  b <- nBlockx * nBlocky

  nSubx <- n1 - nBlockx + 1
  nSuby <- n2 - nBlocky + 1
  count <- nSubx * nSuby

  xySubInd <- as.matrix(expand.grid(
    head(xSeq, -(nBlockx - 1)),
    head(ySeq, -(nBlocky - 1)),
    KEEP.OUT.ATTRS = FALSE
  ))

  freq_obj_sub <- build_freq_objects(nBlockx, nBlocky)
  phi_theta_sub <- build_phi_theta(freq_obj_sub$freqGrid, lag_mat, c_theta)
  phi_pair_sub  <- phi_theta_sub * (phi_theta_sub + phi_theta_sub)

  gEst    <- numeric(nMonteCarlo)
  mHatVar <- numeric(nMonteCarlo)
  bias    <- numeric(nMonteCarlo)
  sigma1  <- numeric(nMonteCarlo)

  n_cores <- get_ncores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  for (q in seq_len(nMonteCarlo)) {
    z <- dataset[, q]

    results <- foreach(
      k = seq_len(count),
      .combine = "c",
      .multicombine = TRUE,
      .export = c("nBlockx", "nBlocky", "xySubInd", "ProInd", "z",
                  "freq_obj_sub", "lag_mat", "g0_vec", "estimate_gamma_spec",
                  "estimate_theta_from_gamma", "periodogramSub","periodogram_general")
    ) %dopar% {
      xSub <- xySubInd[k, 1] + 0:(nBlockx - 1)
      ySub <- xySubInd[k, 2] + 0:(nBlocky - 1)

      keep <- ProInd[, 1] %in% xSub & ProInd[, 2] %in% ySub
      locs_sub <- ProInd[keep, , drop = FALSE]
      z_sub <- z[keep]

      sub_est <- estimate_gamma_spec(z_sub, locs_sub, freq_obj_sub$freqGrid, lag_mat)
      theta_sub_k <- estimate_theta_from_gamma(sub_est$gamma_hat, g0_vec)

      list(list(
        theta_sub_k = theta_sub_k,
        pdgramVals = sub_est$pdgramVals
      ))
    }

    theta_sub <- sapply(results, function(x) x$theta_sub_k)
    gEst[q] <- mean(theta_sub)

    Hsub <- sqrt(b) * (theta_sub - gEst[q])
    mHatVar[q] <- ((count - 1) / count) * var(Hsub)

    bias[q] <- sqrt(b) * mean(theta_sub - gEstimate[q])

    allPdgram <- do.call(cbind, lapply(results, function(x) x$pdgramVals))
    pd_var <- apply(allPdgram, 1, var)

    subSigma <- phi_pair_sub * ((count - 1) / count) * pd_var
    sigma1[q] <- ((2 * pi)^4) * (1 / b) * sum(subSigma)
  }

  stopCluster(cl)

  sigma2 <- mHatVar - sigma1
  if (any(sigma2 < -1e-12, na.rm = TRUE)) {
    warning("Some sigma2 entries are negative beyond roundoff; check numerical stability.")
  }
  sigma2[sigma2 < 0] <- 0

  vCorrect <- var_vHpb + sigma2

  T_hpb <- matrix(NA_real_, nrow = nBoot, ncol = nMonteCarlo)
  for (i in seq_len(nMonteCarlo)) {
    if (var_vHpb[i] <= 0) {
      T_hpb[, i] <- v_hpb[, i] + bias[i]
    } else {
      T_hpb[, i] <- (sqrt(vCorrect[i]) * v_hpb[, i] / sqrt(var_vHpb[i])) + bias[i]
    }
  }

  saveRDS(
    list(
      betaValue = betaValue,
      gEst = gEst,
      mHatVar = mHatVar,
      bias = bias,
      sigma1 = sigma1,
      sigma2 = sigma2,
      vCorrect = vCorrect,
      theta_true = theta_true
    ),
    file = sub_cache_file(betaValue),
    compress = FALSE
  )

  z975 <- qnorm(0.975)

  vhpbQuantile <- apply(v_hpb, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
  ThpbQuantile <- apply(T_hpb, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

  CItHPB <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
  CItHPB[,1] <- gEstimate - (ThpbQuantile[2,] / sqrt(N))
  CItHPB[,2] <- gEstimate - (ThpbQuantile[1,] / sqrt(N))

  CIvHPB <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
  CIvHPB[,1] <- gEstimate - (vhpbQuantile[2,] / sqrt(N))
  CIvHPB[,2] <- gEstimate - (vhpbQuantile[1,] / sqrt(N))

  CIt <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
  CIt[,1] <- gEstimate - z975 * sqrt(mHatVar / N)
  CIt[,2] <- gEstimate + z975 * sqrt(mHatVar / N)

  CIv <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
  CIv[,1] <- gEstimate - z975 * sqrt(var_vHpb / N)
  CIv[,2] <- gEstimate + z975 * sqrt(var_vHpb / N)

  conf_tHPB <- mean((theta_true >= CItHPB[,1]) & (theta_true <= CItHPB[,2]))
  conf_vHPB <- mean((theta_true >= CIvHPB[,1]) & (theta_true <= CIvHPB[,2]))
  conf_t    <- mean((theta_true >= CIt[,1]) & (theta_true <= CIt[,2]))
  conf_v    <- mean((theta_true >= CIv[,1]) & (theta_true <= CIv[,2]))

  diag_vec <- c(
    theta_true      = theta_true,
    mean_theta_hat  = mean(gEstimate),
    mc_var_target   = var(sqrt(N) * (gEstimate - theta_true)),
    mean_mHatVar    = mean(mHatVar),
    mean_sigma1     = mean(sigma1),
    mean_sigma2     = mean(sigma2),
    mean_v_boot     = mean(apply(v_hpb, 2, var)),
    mean_v_formula  = mean(var_vHpb),
    mean_v_corr     = mean(vCorrect),
    mean_bias       = mean(bias),
    mean_boot_ctr   = mean(boot_ctr),
    mc_scaled_var   = mDataVar
  )

  out <- list(
    betaValue = betaValue,
    theta_true = theta_true,
    conf_tHPB = conf_tHPB,
    conf_vHPB = conf_vHPB,
    conf_t = conf_t,
    conf_v = conf_v,
    diagnostics = diag_vec,
    sigma1_mean = mean(sigma1),
    sigma2_mean = mean(sigma2),
    bias_mean = mean(bias),
    mHatVar_mean = mean(mHatVar),
    vCorrect_mean = mean(vCorrect)
  )

  saveRDS(out, file = res_cache_file(betaValue), compress = FALSE)

  print(diag_vec)

  result_line <- paste(
    "True param=", theta_true,
    "block_length =", betaValue,
    "| conf for 95% HFDB:", round(conf_tHPB, 3),
    "| conf for 95% FDWB:", round(conf_vHPB, 3),
    "| conf for normal approx.:", round(conf_t, 3),
    "\n"
  )

  cat(result_line)
  cat("Saved subsample cache to:\n", sub_cache_file(betaValue), "\n")
  cat("Saved final result to:\n", res_cache_file(betaValue), "\n")
}

##
## Run selected mode
##
if (mode == "precompute") {
  precompute_full()

} else if (mode == "subsample") {
  run_subsample_for_beta(betaValue)

} else if (mode == "all") {
  if (!file.exists(full_cache_file)) {
    precompute_full()
  } else {
    cat("Full cache already exists. Skipping precompute.\n")
  }
  run_subsample_for_beta(betaValue)

} else {
  stop("Unknown mode. Use one of: precompute, subsample, all")
}

time_taken <- proc.time() - start_time
cat("Elapsed time (sec):", round(time_taken["elapsed"], 3), "\n")
