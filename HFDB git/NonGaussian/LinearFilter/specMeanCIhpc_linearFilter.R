#!/usr/bin/env Rscript

suppressMessages({
  library(fields)
  library(binhf)
  library(cubature)
  library(Bessel)
  library(astsa)
  library(geoR)
  library(ncf)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(matrixStats)
  library(doRNG)
})

##
## Command line arguments
##
args <- commandArgs(trailingOnly = TRUE)

mode <- if (length(args) >= 1) args[1] else "all"               # precompute / subsample / all
betaValue <- if (length(args) >= 2) as.numeric(args[2]) else 6   # block length
cache_dir <- if (length(args) >= 3) args[3] else "cache_spatial_ci_linear_filter"

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

Delta <- 1
lagh  <- c(1, 0)   # change to c(1,1) or c(0,1) as needed

xSeq <- seq(1, 50, by = Delta)
ySeq <- seq(1, 50, by = Delta)

n1 <- length(xSeq)
n2 <- length(ySeq)
N  <- n1 * n2

grid <- list(x = xSeq, y = ySeq)
Locs <- expand.grid(x = xSeq, y = ySeq, KEEP.OUT.ATTRS = FALSE)
ProInd <- as.matrix(Locs)

## cache file names
lag_tag <- paste0("lag_", lagh[1], "_", lagh[2])
base_tag <- paste0("linear_filter_n1_", n1, "_n2_", n2, "_MC_", nMonteCarlo, "_B_", nBoot, "_", lag_tag)
full_cache_file <- file.path(cache_dir, paste0("full_precompute_", base_tag, ".rds"))
sub_cache_file  <- function(beta) file.path(cache_dir, paste0("subsample_beta_", beta, "_", base_tag, ".rds"))
res_cache_file  <- function(beta) file.path(cache_dir, paste0("final_results_beta_", beta, "_", base_tag, ".rds"))

##
## Non-Gaussian linear filter process
##
rinnov <- function(n) {
  x <- rgamma(n, shape = 0.2, scale = 1)
  (x - 0.2) / sqrt(0.2)
}

A <- matrix(c(
  1.0,  0.3,  0.1,
  0.5,  0.2,  0.1,
  0.25, 0.1,  0.05
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

spec_density <- function(omega, A_mat = A) {
  eval_one <- function(w1, w2) {
    Afreq <-
      A_mat[1,1] +
      A_mat[2,1] * exp(-1i * w1) +
      A_mat[3,1] * exp(-2i * w1) +
      A_mat[1,2] * exp(-1i * w2) +
      A_mat[2,2] * exp(-1i * (w1 + w2)) +
      A_mat[3,2] * exp(-1i * (2 * w1 + w2)) +
      A_mat[1,3] * exp(-2i * w2) +
      A_mat[2,3] * exp(-1i * (w1 + 2 * w2)) +
      A_mat[3,3] * exp(-1i * (2 * w1 + 2 * w2))

    (1 / (2 * pi)^2) * Mod(Afreq)^2
  }

  if (is.vector(omega) && length(omega) == 2) {
    return(eval_one(omega[1], omega[2]))
  }

  omega <- as.matrix(omega)
  out <- numeric(nrow(omega))
  for (k in seq_len(nrow(omega))) {
    out[k] <- eval_one(omega[k, 1], omega[k, 2])
  }
  out
}

gamma_true_fun <- function(A, h1, h2) {
  nr <- nrow(A)
  nc <- ncol(A)

  if (h1 >= nr || h2 >= nc) return(0)

  r1 <- 1:(nr - h1)
  r2 <- (1 + h1):nr
  c1 <- 1:(nc - h2)
  c2 <- (1 + h2):nc

  sum(A[r1, c1] * A[r2, c2])
}

gamma_true <- gamma_true_fun(A, lagh[1], lagh[2])

##
## Periodogram functions
##
periodogram <- function(Z, locs, omega) {
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
## Precompute index maps
##
ProcessIndices <- expand.grid(x = xSeq, y = ySeq, KEEP.OUT.ATTRS = FALSE)

idx_x <- ((ProcessIndices$x - xSeq[1]) / Delta) + 1
idx_y <- ((ProcessIndices$y - ySeq[1]) / Delta) + 1

ProcessMap <- data.frame(
  x = ProcessIndices$x,
  y = ProcessIndices$y,
  i = as.integer(idx_x),
  j = as.integer(idx_y),
  stringsAsFactors = FALSE
)

##
## Full-data frequency grid
##
centered_range <- function(n) {
  if (n %% 2 == 0) {
    seq(-n/2, n/2 - 1, by = 1)
  } else {
    seq(-(n - 1)/2, (n - 1)/2, by = 1)
  }
}

reflect_centered <- function(j, n) {
  if (n %% 2 == 0) {
    ifelse(j == -n/2, -n/2, -j)
  } else {
    -j
  }
}

j1_range <- centered_range(n1)
j2_range <- centered_range(n2)

f1 <- (2 * pi * j1_range) / n1
f2 <- (2 * pi * j2_range) / n2

freqGridData_full <- expand.grid(f1 = f1, f2 = f2, KEEP.OUT.ATTRS = FALSE)
freqGridData_full <- as.matrix(freqGridData_full)

zero_freq_logical <- (freqGridData_full[,1] == 0 & freqGridData_full[,2] == 0)
freqGridData <- freqGridData_full[!zero_freq_logical, , drop = FALSE]

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

if (anyNA(partner_ind)) stop("Symmetric partner not found for some frequency.")

self_ind <- which(partner_ind == seq_len(nrow(Jn)))
pos_ind  <- which(seq_len(nrow(Jn)) < partner_ind)
neg_ind  <- partner_ind[pos_ind]
resid_ind <- setdiff(self_ind, zero_ind)

flat_order <- order(Jn$row, Jn$col)
use_ind <- setdiff(flat_order, zero_ind)

if (length(use_ind) != nrow(freqGridData)) {
  stop("Mismatch between use_ind and freqGridData rows.")
}

##
## 2D smoothing helpers
##
kernel_2d <- function(u1, u2, type = c("epanechnikov", "gaussian")) {
  type <- match.arg(type)

  if (type == "epanechnikov") {
    out <- (9/16) * (1 - u1^2) * (1 - u2^2)
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

##
## Smoother setup (depends only on full grid)
##
hx <- 0.05
hy <- 0.05
kernel_type <- "gaussian"

weight_obj <- build_literal_kernel_weights(
  f1 = f1,
  f2 = f2,
  hx = hx,
  hy = hy,
  kernel = kernel_type,
  periodic = TRUE
)

##
## Full precompute stage
##
precompute_full <- function() {
  cat("Starting full precomputation for non-Gaussian linear filter field...\n")

  dataset <- matrix(NA_real_, nrow = N, ncol = nMonteCarlo)
  for (i in seq_len(nMonteCarlo)) {
    dataset[, i] <- c(simulate_ma22_field(n1, n2, A, rinnov))
  }

  ## ---- full-sample estimates ----
  n_cores <- get_ncores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  hLambda_full <- as.vector(freqGridData %*% lagh)
  cos_h_full <- cos(hLambda_full)

  gEstimate <- foreach(
    q = seq_len(nMonteCarlo),
    .combine = c,
    .packages = "fields",
    .export = c("dataset", "grid", "ProcessMap", "freqGridData", "periodogram",
                "lagh", "N", "ProcessIndices")
  ) %dopar% {
    gProcess <- dataset[, q]
    Gprocess <- as.surface(grid, gProcess)
    GprocessData <- as.vector(Gprocess$z[cbind(ProcessMap$i, ProcessMap$j)])

    pdgramVals <- apply(
      freqGridData, 1,
      function(omega) periodogram(GprocessData, as.matrix(ProcessIndices), omega)
    )

    ((2 * pi)^2) * (1 / N) * sum(cos_h_full * pdgramVals)
  }

  mDataVar <- (N / nMonteCarlo) * sum((gEstimate - mean(gEstimate))^2)

  stopCluster(cl)

  ## ---- full-sample bootstrap objects ----
  n_cores <- get_ncores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  doRNG::registerDoRNG(123)

  v_hpb    <- matrix(NA_real_, nrow = nBoot, ncol = nMonteCarlo)
  var_vHpb <- numeric(nMonteCarlo)

  hLambdaData <- lagh[1] * freqGridData[,1] + lagh[2] * freqGridData[,2]
  phiFunc <- cos(hLambdaData)
  phiFunctional <- cos(hLambdaData) * (cos(hLambdaData) + cos(-hLambdaData))

  for (q in seq_len(nMonteCarlo)) {
    gProcess <- as.matrix(dataset[, q])
    Gprocess <- as.surface(grid, gProcess)

    xyDataInd <- cbind(
      ((ProInd[,1] - xSeq[1]) / Delta) + 1,
      ((ProInd[,2] - ySeq[1]) / Delta) + 1
    )
    GprocessData <- Gprocess$z[cbind(xyDataInd[,1], xyDataInd[,2])]

    vals <- apply(
      freqGridData, 1,
      function(omega) periodogram(GprocessData, ProInd, omega)
    )

    pMatrix <- matrix(0, nrow = length(f2), ncol = length(f1))
    pMatrix[cbind(Jn$row[use_ind], Jn$col[use_ind])] <- vals
    pMatrix[Jn$row[zero_ind], Jn$col[zero_ind]] <- 0

    mat <- smooth_periodogram_fast(
      pMatrix = pMatrix,
      weight_obj = weight_obj,
      zero_row = Jn$row[zero_ind],
      zero_col = Jn$col[zero_ind]
    )

    smoothPdgram <- mat[cbind(Jn$row[use_ind], Jn$col[use_ind])]

    v_hpb_boot <- foreach(
      l = seq_len(nBoot),
      .combine = c,
      .export = c("phiFunc", "smoothPdgram", "Jn", "pos_ind", "neg_ind",
                  "resid_ind", "use_ind", "N")
    ) %dorng% {
      expVals <- numeric(nrow(Jn))

      epos <- rexp(length(pos_ind), rate = 1) - 1
      expVals[pos_ind] <- epos
      expVals[neg_ind] <- epos

      if (length(resid_ind) > 0) {
        expVals[resid_ind] <- rexp(length(resid_ind), rate = 1) - 1
      }

      expVector <- expVals[use_ind]
      vPart <- phiFunc * smoothPdgram * expVector
      sqrt(N) * ((2 * pi)^2 / N) * sum(vPart)
    }

    v_hpb[, q] <- v_hpb_boot
    var_vHpb[q] <- ((2 * pi)^4) * (1 / N) * sum(phiFunctional * smoothPdgram^2)
  }

  stopCluster(cl)

  saveRDS(
    list(
      dataset = dataset,
      gEstimate = gEstimate,
      mDataVar = mDataVar,
      v_hpb = v_hpb,
      var_vHpb = var_vHpb,
      gamma_true = gamma_true,
      params = list(
        nMonteCarlo = nMonteCarlo,
        nBoot = nBoot,
        n1 = n1, n2 = n2, N = N,
        Delta = Delta,
        lagh = lagh,
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

  dataset    <- full_obj$dataset
  gEstimate  <- full_obj$gEstimate
  v_hpb      <- full_obj$v_hpb
  var_vHpb   <- full_obj$var_vHpb
  gamma_true <- full_obj$gamma_true

  ## ---- block setup depending on beta ----
  nBlockx <- betaValue
  nBlocky <- betaValue
  b <- betaValue^2

  nSubx <- n1 - nBlockx + 1
  nSuby <- n2 - nBlocky + 1
  count <- nSubx * nSuby

  if (nSubx < 1 || nSuby < 1) {
    stop("Block size too large for the grid.")
  }

  j1_sub <- seq(floor(-(nBlockx - 1)/2), (nBlockx - floor(nBlockx/2)), by = 1)
  j2_sub <- seq(floor(-(nBlocky - 1)/2), (nBlocky - floor(nBlocky/2)), by = 1)

  j1_sub <- (2 * pi * j1_sub) / nBlockx
  j2_sub <- (2 * pi * j2_sub) / nBlocky

  freqGrid_sub_full <- expand.grid(j1 = j1_sub, j2 = j2_sub, KEEP.OUT.ATTRS = FALSE)
  freqGrid_sub_full <- as.matrix(freqGrid_sub_full)

  freqGrid_sub <- freqGrid_sub_full[
    !(freqGrid_sub_full[,1] == 0 & freqGrid_sub_full[,2] == 0),
    , drop = FALSE
  ]

  xySubInd <- as.matrix(expand.grid(
    head(xSeq, -(nBlockx - 1)),
    head(ySeq, -(nBlocky - 1)),
    KEEP.OUT.ATTRS = FALSE
  ))

  hLambda <- freqGrid_sub %*% lagh
  cos_h <- cos(hLambda)
  cos_pair <- cos_h * (cos_h + cos(-hLambda))

  ## ---- outputs ----
  gEst    <- numeric(nMonteCarlo)
  mHatVar <- numeric(nMonteCarlo)
  bias    <- numeric(nMonteCarlo)
  sigma1  <- numeric(nMonteCarlo)

  ## ---- parallel backend ----
  n_cores <- get_ncores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  for (q in seq_len(nMonteCarlo)) {
    gProcess <- as.matrix(dataset[, q])
    Gprocess <- as.surface(grid, gProcess)

    results <- foreach(
      k = seq_len(count),
      .packages = "fields",
      .export = c("freqGrid_sub", "Gprocess", "ProcessIndices", "Delta",
                  "nBlockx", "nBlocky", "xySubInd", "b", "cos_h",
                  "periodogramSub", "cos_pair", "xSeq", "ySeq")
    ) %dopar% {
      xSub <- xySubInd[k,1] + 0:(nBlockx - 1)
      ySub <- xySubInd[k,2] + 0:(nBlocky - 1)

      keep <- ProcessIndices[,1] %in% xSub & ProcessIndices[,2] %in% ySub
      xySubProcess <- ProcessIndices[keep, , drop = FALSE]

      xInd <- ((xySubProcess[,1] - xSeq[1]) / Delta) + 1
      yInd <- ((xySubProcess[,2] - ySeq[1]) / Delta) + 1
      GprocessSub <- Gprocess$z[cbind(xInd, yInd)]

      pdgramVals <- apply(
        freqGrid_sub, 1,
        function(omega) periodogramSub(GprocessSub, xySubProcess, omega)
      )

      mHatSub_k <- ((2 * pi)^2) * (1 / b) * sum(cos_h * pdgramVals)

      list(mHatSub_k = mHatSub_k, pdgramVals = pdgramVals)
    }

    mHatSub <- sapply(results, function(x) x$mHatSub_k)
    gEst[q] <- mean(mHatSub)

    Hsub <- sqrt(b) * (mHatSub - gEst[q])
    mHatVar[q] <- (count - 1) / count * var(Hsub)

    biasSub <- mHatSub - gEstimate[q]
    bias[q] <- sqrt(b) * mean(biasSub)

    allPdgram <- do.call(cbind, lapply(results, function(x) x$pdgramVals))
    pd_var <- apply(allPdgram, 1, var)

    subSigma <- cos_pair * ((count - 1) / count) * pd_var
    sigma1[q] <- ((2 * pi)^4) * (1 / b) * sum(subSigma)
  }

  stopCluster(cl)

  sigma2 <- mHatVar - sigma1
  if (any(sigma2 < -1e-12, na.rm = TRUE)) {
    warning("Some sigma2 values are negative beyond roundoff.")
  }
  sigma2[sigma2 < 0] <- 0

  saveRDS(
    list(
      betaValue = betaValue,
      b = b,
      count = count,
      gEst = gEst,
      mHatVar = mHatVar,
      bias = bias,
      sigma1 = sigma1,
      sigma2 = sigma2
    ),
    file = sub_cache_file(betaValue),
    compress = FALSE
  )

  ## ---- corrected CI calculation ----
  varCorrection <- sigma2
  vCorrect <- var_vHpb + varCorrection

  T_hpb <- matrix(NA_real_, nrow = nBoot, ncol = nMonteCarlo)
  for (i in seq_len(nMonteCarlo)) {
    T_hpb[, i] <- (sqrt(vCorrect[i]) * v_hpb[, i] / sqrt(var_vHpb[i])) + bias[i]
  }

  vhpbQuantile <- apply(v_hpb, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
  ThpbQuantile <- apply(T_hpb, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

  CItHPB <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
  CItHPB[,1] <- gEstimate - (ThpbQuantile[2, ] / sqrt(N))
  CItHPB[,2] <- gEstimate - (ThpbQuantile[1, ] / sqrt(N))

  CIvHPB <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
  CIvHPB[,1] <- gEstimate - (vhpbQuantile[2, ] / sqrt(N))
  CIvHPB[,2] <- gEstimate - (vhpbQuantile[1, ] / sqrt(N))

  z975 <- qnorm(0.975)

  CIt <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
  CIt[,1] <- gEstimate - z975 * sqrt(mHatVar / N)
  CIt[,2] <- gEstimate + z975 * sqrt(mHatVar / N)

  CIv <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
  CIv[,1] <- gEstimate - z975 * sqrt(var_vHpb / N)
  CIv[,2] <- gEstimate + z975 * sqrt(var_vHpb / N)

  conf_tHPB <- mean((gamma_true >= CItHPB[,1]) & (gamma_true <= CItHPB[,2]))
  conf_vHPB <- mean((gamma_true >= CIvHPB[,1]) & (gamma_true <= CIvHPB[,2]))
  conf_t    <- mean((gamma_true >= CIt[,1])    & (gamma_true <= CIt[,2]))
  conf_v    <- mean((gamma_true >= CIv[,1])    & (gamma_true <= CIv[,2]))

  out <- list(
    betaValue = betaValue,
    gamma_true = gamma_true,
    conf_tHPB = conf_tHPB,
    conf_vHPB = conf_vHPB,
    conf_t = conf_t,
    conf_v = conf_v,
    sigma1_mean = mean(sigma1),
    sigma2_mean = mean(sigma2),
    bias_mean = mean(bias),
    mHatVar_mean = mean(mHatVar),
    vCorrect_mean = mean(vCorrect)
  )

  saveRDS(out, file = res_cache_file(betaValue), compress = FALSE)

  result_line <- paste(
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
