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

mode <- if (length(args) >= 1) args[1] else "all"  # precompute / subsample / all
cache_dir <- if (length(args) >= 2) args[2] else "cache_spatial_ci_linear_filter_bw"

if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cat("Mode      :", mode, "\n")
cat("cache_dir :", cache_dir, "\n")

start_time <- proc.time()

## 
## Utility: cores for HPC
## 
get_ncores <- function(max_cores = 64L) {
  slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA))
  if (!is.na(slurm_cores) && slurm_cores >= 1) {
    return(slurm_cores)
  }
  n_cores <- min(max_cores, parallel::detectCores(logical = TRUE))
  if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
  n_cores
}

make_hpc_cluster <- function(n_cores = get_ncores()) {
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  cl
}

## 
## Global parameters
##
set.seed(2)

nMonteCarlo <- 1000
nBoot       <- 500

Delta <- 1
lagh  <- c(1, 0)

n1 <- 50
n2 <- 50
N  <- n1 * n2

block_side <- 10
nBlockx <- block_side
nBlocky <- block_side
b <- nBlockx * nBlocky

bandwidth_grid <- c(0.05, 0.10, 0.15, 0.20, 0.25)
kernel_type <- "gaussian"

xSeq <- seq(1, n1, by = Delta)
ySeq <- seq(1, n2, by = Delta)

grid <- list(x = xSeq, y = ySeq)
Locs <- expand.grid(x = xSeq, y = ySeq, KEEP.OUT.ATTRS = FALSE)
ProInd <- as.matrix(Locs)
ProcessIndices <- ProInd

## 
## Cache file names
##
bw_tag <- function(h) {
  paste0("bw_", gsub("\\.", "p", formatC(h, format = "f", digits = 2)))
}

lag_tag <- paste0("lag_", lagh[1], "_", lagh[2])
base_tag <- paste0(
  "linear_filter_n1_", n1,
  "_n2_", n2,
  "_MC_", nMonteCarlo,
  "_B_", nBoot,
  "_block_", block_side,
  "_", lag_tag
)

common_cache_file <- file.path(cache_dir, paste0("common_precompute_", base_tag, ".rds"))
bootstrap_cache_file <- function(h) file.path(cache_dir, paste0("bootstrap_", bw_tag(h), "_", base_tag, ".rds"))
subsample_cache_file <- file.path(cache_dir, paste0("subsample_fixedblock_", base_tag, ".rds"))
final_result_file <- function(h) file.path(cache_dir, paste0("final_result_", bw_tag(h), "_", base_tag, ".rds"))
summary_result_file <- file.path(cache_dir, paste0("summary_bandwidth_scan_", base_tag, ".rds"))

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
periodogram_general <- function(Z, locs, omega) {
  Z <- as.numeric(Z)
  locs <- as.matrix(locs)
  omega <- as.numeric(omega)
  n <- nrow(locs)

  phase <- locs %*% omega
  dft <- sum(Z * exp(-1i * phase))
  (2 * pi)^(-2) * (1 / n) * Mod(dft)^2
}

##
## Index maps
##
idx_x <- ((ProcessIndices[, 1] - xSeq[1]) / Delta) + 1
idx_y <- ((ProcessIndices[, 2] - ySeq[1]) / Delta) + 1

ProcessMap <- data.frame(
  x = ProcessIndices[, 1],
  y = ProcessIndices[, 2],
  i = as.integer(idx_x),
  j = as.integer(idx_y),
  stringsAsFactors = FALSE
)

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

build_full_freq_objects <- function(n1, n2) {
  j1_range <- centered_range(n1)
  j2_range <- centered_range(n2)

  f1 <- (2 * pi * j1_range) / n1
  f2 <- (2 * pi * j2_range) / n2

  freqGrid_full <- as.matrix(expand.grid(f1 = f1, f2 = f2, KEEP.OUT.ATTRS = FALSE))
  zero_freq_logical <- (freqGrid_full[, 1] == 0 & freqGrid_full[, 2] == 0)
  freqGrid <- freqGrid_full[!zero_freq_logical, , drop = FALSE]

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

  list(
    j1_range = j1_range,
    j2_range = j2_range,
    f1 = f1,
    f2 = f2,
    freqGridData_full = freqGrid_full,
    freqGridData = freqGrid,
    Jn = Jn,
    zero_ind = zero_ind,
    pos_ind = pos_ind,
    neg_ind = neg_ind,
    resid_ind = resid_ind,
    use_ind = use_ind
  )
}

full_freq_obj <- build_full_freq_objects(n1, n2)
freqGridData <- full_freq_obj$freqGridData
Jn_full <- full_freq_obj$Jn
use_ind_full <- full_freq_obj$use_ind
zero_ind_full <- full_freq_obj$zero_ind

if (length(use_ind_full) != nrow(freqGridData)) {
  stop("Mismatch between full-grid indexing and nonzero frequencies.")
}

build_sub_freq_objects <- function(nBlockx, nBlocky) {
  j1_sub <- seq(floor(-(nBlockx - 1) / 2), (nBlockx - floor(nBlockx / 2)), by = 1)
  j2_sub <- seq(floor(-(nBlocky - 1) / 2), (nBlocky - floor(nBlocky / 2)), by = 1)

  f1_sub <- (2 * pi * j1_sub) / nBlockx
  f2_sub <- (2 * pi * j2_sub) / nBlocky

  freqGrid_sub_full <- as.matrix(expand.grid(j1 = f1_sub, j2 = f2_sub, KEEP.OUT.ATTRS = FALSE))
  freqGrid_sub <- freqGrid_sub_full[
    !(freqGrid_sub_full[, 1] == 0 & freqGrid_sub_full[, 2] == 0),
    , drop = FALSE
  ]

  list(
    freqGrid_sub = freqGrid_sub,
    f1_sub = f1_sub,
    f2_sub = f2_sub
  )
}

sub_freq_obj <- build_sub_freq_objects(nBlockx, nBlocky)
freqGrid_sub <- sub_freq_obj$freqGrid_sub

xySubInd <- as.matrix(expand.grid(
  head(xSeq, -(nBlockx - 1)),
  head(ySeq, -(nBlocky - 1)),
  KEEP.OUT.ATTRS = FALSE
))

nSubx <- n1 - nBlockx + 1
nSuby <- n2 - nBlocky + 1
count <- nSubx * nSuby

if (count != nrow(xySubInd)) {
  stop("Subsampling index count mismatch.")
}

##
## 2D smoothing helpers
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

vector_to_periodogram_matrix <- function(vals, Jn, use_ind, zero_ind) {
  pMatrix <- matrix(0, nrow = max(Jn$row), ncol = max(Jn$col))
  pMatrix[cbind(Jn$row[use_ind], Jn$col[use_ind])] <- vals
  pMatrix[Jn$row[zero_ind], Jn$col[zero_ind]] <- 0
  pMatrix
}

##
## Common precompute stage: data + full periodograms + gEstimate
##
precompute_common <- function() {
  if (file.exists(common_cache_file)) {
    cat("Common cache already exists. Using existing file:\n", common_cache_file, "\n")
    return(invisible(NULL))
  }

  cat("Starting common precomputation for 50x50 non-Gaussian linear filter field...\n")

  dataset <- matrix(NA_real_, nrow = N, ncol = nMonteCarlo)
  for (i in seq_len(nMonteCarlo)) {
    dataset[, i] <- c(simulate_ma22_field(n1, n2, A, rinnov))
  }

  hLambda_full <- as.vector(freqGridData %*% lagh)
  cos_h_full <- cos(hLambda_full)

  n_cores <- get_ncores()
  cl <- make_hpc_cluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)

  per_mc <- foreach(
    q = seq_len(nMonteCarlo),
    .packages = "fields",
    .export = c("dataset", "grid", "ProcessMap", "freqGridData", "periodogram_general",
                "lagh", "N", "ProcessIndices")
  ) %dopar% {
    gProcess <- dataset[, q]
    Gprocess <- as.surface(grid, gProcess)
    GprocessData <- as.vector(Gprocess$z[cbind(ProcessMap$i, ProcessMap$j)])

    pdgramVals <- apply(
      freqGridData, 1,
      function(omega) periodogram_general(GprocessData, ProcessIndices, omega)
    )

    gEstimate_q <- ((2 * pi)^2) * (1 / N) * sum(cos_h_full * pdgramVals)

    list(
      pdgramVals = pdgramVals,
      gEstimate = gEstimate_q
    )
  }

  pdgram_full <- do.call(cbind, lapply(per_mc, `[[`, "pdgramVals"))
  gEstimate <- vapply(per_mc, `[[`, numeric(1), "gEstimate")
  mDataVar <- (N / nMonteCarlo) * sum((gEstimate - mean(gEstimate))^2)

  saveRDS(
    list(
      dataset = dataset,
      pdgram_full = pdgram_full,
      gEstimate = gEstimate,
      mDataVar = mDataVar,
      gamma_true = gamma_true,
      params = list(
        nMonteCarlo = nMonteCarlo,
        nBoot = nBoot,
        n1 = n1,
        n2 = n2,
        N = N,
        Delta = Delta,
        lagh = lagh,
        block_side = block_side,
        b = b,
        bandwidth_grid = bandwidth_grid,
        kernel_type = kernel_type,
        A = A
      )
    ),
    file = common_cache_file,
    compress = FALSE
  )

  cat("Saved common cache to:\n", common_cache_file, "\n")
}

##
## Bootstrap precompute for all bandwidths
##
precompute_bootstrap_for_all_bandwidths <- function() {
  if (!file.exists(common_cache_file)) {
    stop("Common cache file not found. Run precompute_common() first.")
  }

  common_obj <- readRDS(common_cache_file)
  pdgram_full <- common_obj$pdgram_full

  hLambdaData <- lagh[1] * freqGridData[, 1] + lagh[2] * freqGridData[, 2]
  phiFunc <- cos(hLambdaData)
  phiFunctional <- cos(hLambdaData) * (cos(hLambdaData) + cos(-hLambdaData))

  for (h in bandwidth_grid) {
    out_file <- bootstrap_cache_file(h)

    if (file.exists(out_file)) {
      cat("Bootstrap cache already exists for bandwidth", h, ":\n", out_file, "\n")
      next
    }

    cat("Starting bootstrap precompute for bandwidth =", h, "\n")
    weight_obj <- build_literal_kernel_weights(
      f1 = full_freq_obj$f1,
      f2 = full_freq_obj$f2,
      hx = h,
      hy = h,
      kernel = kernel_type,
      periodic = TRUE
    )

    pMatrix_list <- vector("list", nMonteCarlo)
    for (q in seq_len(nMonteCarlo)) {
      pMatrix_list[[q]] <- vector_to_periodogram_matrix(
        vals = pdgram_full[, q],
        Jn = Jn_full,
        use_ind = use_ind_full,
        zero_ind = zero_ind_full
      )
    }

    n_cores <- get_ncores()
    cl <- make_hpc_cluster(n_cores)
    doRNG::registerDoRNG(123)

    per_mc_bw <- foreach(
      q = seq_len(nMonteCarlo),
      .packages = character(),
      .export = c("pMatrix_list", "weight_obj", "smooth_periodogram_fast", "Jn_full",
                  "use_ind_full", "zero_ind_full", "phiFunc", "phiFunctional", "nBoot",
                  "N", "full_freq_obj")
    ) %dorng% {
      pMatrix <- pMatrix_list[[q]]

      smoothed_matrix <- smooth_periodogram_fast(
        pMatrix = pMatrix,
        weight_obj = weight_obj,
        zero_row = Jn_full$row[zero_ind_full],
        zero_col = Jn_full$col[zero_ind_full]
      )

      smoothPdgram <- smoothed_matrix[cbind(Jn_full$row[use_ind_full], Jn_full$col[use_ind_full])]

      v_hpb_boot <- numeric(nBoot)
      for (l in seq_len(nBoot)) {
        expVals <- numeric(nrow(Jn_full))

        epos <- rexp(length(full_freq_obj$pos_ind), rate = 1) - 1
        expVals[full_freq_obj$pos_ind] <- epos
        expVals[full_freq_obj$neg_ind] <- epos

        if (length(full_freq_obj$resid_ind) > 0) {
          expVals[full_freq_obj$resid_ind] <- rexp(length(full_freq_obj$resid_ind), rate = 1) - 1
        }

        expVector <- expVals[use_ind_full]
        vPart <- phiFunc * smoothPdgram * expVector
        v_hpb_boot[l] <- sqrt(N) * ((2 * pi)^2 / N) * sum(vPart)
      }

      var_vHpb_q <- ((2 * pi)^4) * (1 / N) * sum(phiFunctional * smoothPdgram^2)

      list(
        v_hpb = v_hpb_boot,
        var_vHpb = var_vHpb_q
      )
    }

    stopCluster(cl)

    v_hpb <- do.call(cbind, lapply(per_mc_bw, `[[`, "v_hpb"))
    var_vHpb <- vapply(per_mc_bw, `[[`, numeric(1), "var_vHpb")

    saveRDS(
      list(
        bandwidth = h,
        hx = h,
        hy = h,
        v_hpb = v_hpb,
        var_vHpb = var_vHpb
      ),
      file = out_file,
      compress = FALSE
    )

    cat("Saved bootstrap cache for bandwidth", h, "to:\n", out_file, "\n")
  }
}

##
## Fixed-block subsampling stage (independent of bandwidth)
##
run_fixed_block_subsampling <- function() {
  if (file.exists(subsample_cache_file)) {
    cat("Fixed-block subsample cache already exists. Using existing file:\n", subsample_cache_file, "\n")
    return(invisible(NULL))
  }

  if (!file.exists(common_cache_file)) {
    stop("Common cache file not found. Run precompute stage first.")
  }

  common_obj <- readRDS(common_cache_file)
  dataset <- common_obj$dataset
  gEstimate <- common_obj$gEstimate

  hLambda_sub <- freqGrid_sub %*% lagh
  cos_h_sub <- cos(hLambda_sub)
  cos_pair_sub <- cos_h_sub * (cos_h_sub + cos(-hLambda_sub))

  gEst <- numeric(nMonteCarlo)
  mHatVar <- numeric(nMonteCarlo)
  bias <- numeric(nMonteCarlo)
  sigma1 <- numeric(nMonteCarlo)

  n_cores <- get_ncores()
  cl <- make_hpc_cluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)

  for (q in seq_len(nMonteCarlo)) {
    gProcess <- as.matrix(dataset[, q])
    Gprocess <- as.surface(grid, gProcess)

    results <- foreach(
      k = seq_len(count),
      .packages = "fields",
      .export = c("freqGrid_sub", "Gprocess", "ProcessIndices", "Delta",
                  "nBlockx", "nBlocky", "xySubInd", "b", "cos_h_sub",
                  "periodogram_general", "cos_pair_sub", "xSeq", "ySeq")
    ) %dopar% {
      xSub <- xySubInd[k, 1] + 0:(nBlockx - 1)
      ySub <- xySubInd[k, 2] + 0:(nBlocky - 1)

      keep <- ProcessIndices[, 1] %in% xSub & ProcessIndices[, 2] %in% ySub
      xySubProcess <- ProcessIndices[keep, , drop = FALSE]

      xInd <- ((xySubProcess[, 1] - xSeq[1]) / Delta) + 1
      yInd <- ((xySubProcess[, 2] - ySeq[1]) / Delta) + 1
      GprocessSub <- Gprocess$z[cbind(xInd, yInd)]

      pdgramVals <- apply(
        freqGrid_sub, 1,
        function(omega) periodogram_general(GprocessSub, xySubProcess, omega)
      )

      mHatSub_k <- ((2 * pi)^2) * (1 / b) * sum(cos_h_sub * pdgramVals)

      list(
        mHatSub_k = mHatSub_k,
        pdgramVals = pdgramVals
      )
    }

    mHatSub <- vapply(results, `[[`, numeric(1), "mHatSub_k")
    gEst[q] <- mean(mHatSub)

    Hsub <- sqrt(b) * (mHatSub - gEst[q])
    mHatVar[q] <- (count - 1) / count * var(Hsub)

    biasSub <- mHatSub - gEstimate[q]
    bias[q] <- sqrt(b) * mean(biasSub)

    allPdgram <- do.call(cbind, lapply(results, `[[`, "pdgramVals"))
    pd_var <- apply(allPdgram, 1, var)

    subSigma <- cos_pair_sub * ((count - 1) / count) * pd_var
    sigma1[q] <- ((2 * pi)^4) * (1 / b) * sum(subSigma)
  }

  sigma2 <- mHatVar - sigma1
  if (any(sigma2 < -1e-12, na.rm = TRUE)) {
    warning("Some sigma2 values are negative beyond roundoff.")
  }
  sigma2[sigma2 < 0] <- 0

  saveRDS(
    list(
      block_side = block_side,
      b = b,
      count = count,
      gEst = gEst,
      mHatVar = mHatVar,
      bias = bias,
      sigma1 = sigma1,
      sigma2 = sigma2
    ),
    file = subsample_cache_file,
    compress = FALSE
  )

  cat("Saved fixed-block subsample cache to:\n", subsample_cache_file, "\n")
}

##
## Finalize results for each bandwidth
##
finalize_results_all_bandwidths <- function() {
  if (!file.exists(common_cache_file)) {
    stop("Common cache file not found.")
  }
  if (!file.exists(subsample_cache_file)) {
    stop("Subsample cache file not found.")
  }

  common_obj <- readRDS(common_cache_file)
  sub_obj <- readRDS(subsample_cache_file)

  gEstimate <- common_obj$gEstimate
  gamma_true <- common_obj$gamma_true
  bias <- sub_obj$bias
  sigma1 <- sub_obj$sigma1
  sigma2 <- sub_obj$sigma2
  mHatVar <- sub_obj$mHatVar

  summary_list <- vector("list", length(bandwidth_grid))

  for (i in seq_along(bandwidth_grid)) {
    h <- bandwidth_grid[i]
    boot_file <- bootstrap_cache_file(h)

    if (!file.exists(boot_file)) {
      stop("Missing bootstrap cache for bandwidth ", h, ".")
    }

    boot_obj <- readRDS(boot_file)
    v_hpb <- boot_obj$v_hpb
    var_vHpb <- boot_obj$var_vHpb

    varCorrection <- sigma2
    vCorrect <- var_vHpb + varCorrection

    T_hpb <- matrix(NA_real_, nrow = nBoot, ncol = nMonteCarlo)
    for (j in seq_len(nMonteCarlo)) {
      T_hpb[, j] <- (sqrt(vCorrect[j]) * v_hpb[, j] / sqrt(var_vHpb[j])) + bias[j]
    }

    vhpbQuantile <- apply(v_hpb, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
    ThpbQuantile <- apply(T_hpb, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

    CItHPB <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
    CItHPB[, 1] <- gEstimate - (ThpbQuantile[2, ] / sqrt(N))
    CItHPB[, 2] <- gEstimate - (ThpbQuantile[1, ] / sqrt(N))

    CIvHPB <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
    CIvHPB[, 1] <- gEstimate - (vhpbQuantile[2, ] / sqrt(N))
    CIvHPB[, 2] <- gEstimate - (vhpbQuantile[1, ] / sqrt(N))

    z975 <- qnorm(0.975)

    CIt <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
    CIt[, 1] <- gEstimate - z975 * sqrt(mHatVar / N)
    CIt[, 2] <- gEstimate + z975 * sqrt(mHatVar / N)

    CIv <- matrix(NA_real_, nrow = nMonteCarlo, ncol = 2)
    CIv[, 1] <- gEstimate - z975 * sqrt(var_vHpb / N)
    CIv[, 2] <- gEstimate + z975 * sqrt(var_vHpb / N)

    conf_tHPB <- mean((gamma_true >= CItHPB[, 1]) & (gamma_true <= CItHPB[, 2]))
    conf_vHPB <- mean((gamma_true >= CIvHPB[, 1]) & (gamma_true <= CIvHPB[, 2]))
    conf_t    <- mean((gamma_true >= CIt[, 1])    & (gamma_true <= CIt[, 2]))
    conf_v    <- mean((gamma_true >= CIv[, 1])    & (gamma_true <= CIv[, 2]))

    out <- list(
      bandwidth = h,
      hx = h,
      hy = h,
      block_side = block_side,
      lagh = lagh,
      gamma_true = gamma_true,
      conf_tHPB = conf_tHPB,
      conf_vHPB = conf_vHPB,
      conf_t = conf_t,
      conf_v = conf_v,
      sigma1_mean = mean(sigma1),
      sigma2_mean = mean(sigma2),
      bias_mean = mean(bias),
      mHatVar_mean = mean(mHatVar),
      var_vHpb_mean = mean(var_vHpb),
      vCorrect_mean = mean(vCorrect)
    )

    saveRDS(out, file = final_result_file(h), compress = FALSE)
    summary_list[[i]] <- out

    result_line <- paste(
      "bandwidth =", sprintf("%.2f", h),
      "| block =", block_side,
      "| lag = (", lagh[1], ",", lagh[2], ")",
      "| conf for 95% HFDB:", round(conf_tHPB, 3),
      "| conf for 95% FDWB:", round(conf_vHPB, 3),
      "| conf for normal approx.:", round(conf_t, 3)
    )
    cat(result_line, "\n")
  }

  summary_df <- do.call(rbind, lapply(summary_list, function(x) {
    data.frame(
      bandwidth = x$bandwidth,
      conf_tHPB = x$conf_tHPB,
      conf_vHPB = x$conf_vHPB,
      conf_t = x$conf_t,
      conf_v = x$conf_v,
      sigma1_mean = x$sigma1_mean,
      sigma2_mean = x$sigma2_mean,
      bias_mean = x$bias_mean,
      mHatVar_mean = x$mHatVar_mean,
      var_vHpb_mean = x$var_vHpb_mean,
      vCorrect_mean = x$vCorrect_mean
    )
  }))

  saveRDS(
    list(
      summary = summary_df,
      details = summary_list
    ),
    file = summary_result_file,
    compress = FALSE
  )

  cat("Saved bandwidth summary to:\n", summary_result_file, "\n")
  print(summary_df)
}

##
## Run selected mode
##
if (mode == "precompute") {
  precompute_common()
  precompute_bootstrap_for_all_bandwidths()

} else if (mode == "subsample") {
  run_fixed_block_subsampling()
  finalize_results_all_bandwidths()

} else if (mode == "all") {
  precompute_common()
  precompute_bootstrap_for_all_bandwidths()
  run_fixed_block_subsampling()
  finalize_results_all_bandwidths()

} else {
  stop("Unknown mode. Use one of: precompute, subsample, all")
}

time_taken <- proc.time() - start_time
cat("Elapsed time (sec):", round(time_taken["elapsed"], 3), "\n")

