#!/usr/bin/env Rscript
#Spherical covariance + quadratic non-Gaussian construction

suppressMessages({
  library(fields)
  library(parallel)
  library(matrixStats)
  library(doParallel)
  library(foreach)
  library(doRNG)
})

args <- commandArgs(trailingOnly = TRUE)
psiR <- if (length(args) >= 1) as.numeric(args[1]) else 1
out_dir <- if (length(args) >= 2) args[2] else file.path(getwd(), paste0("psiR_", format(psiR, trim = TRUE, scientific = FALSE)))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

get_ncores <- function(max_cores = 6L) {
  slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA))
  if (!is.na(slurm_cores) && slurm_cores >= 1L) return(slurm_cores)
  nc <- min(max_cores, parallel::detectCores(logical = TRUE))
  if (is.na(nc) || nc < 1L) nc <- 1L
  nc
}

start_time <- proc.time()
set.seed(2)

nMonteCarlo <- 1000
nBoot       <- 500
Delta       <- 1
n1          <- 50
n2          <- 50
xSeq        <- seq_len(n1)
ySeq        <- seq_len(n2)
N           <- n1 * n2
betaValue   <- 11
nBlockx     <- betaValue
nBlocky     <- betaValue
b           <- nBlockx * nBlocky
hx          <- 0.15
hy          <- 0.15
kernel_type <- "gaussian"
sigma2_par  <- 3
phi_par     <- 6
eta_par     <- 2
psi_A       <- 0
alpha_levels <- c(0.10, 0.05)
alpha_ng    <- as.numeric(Sys.getenv("ALPHA_NG", unset = "0.5"))
if (is.na(alpha_ng)) alpha_ng <- 0.5
h_i <- c(1, 0)
h_j <- c(0, 1)

make_aniso_B <- function(psi_A, psi_R) {
  R <- matrix(c(cos(psi_A), sin(psi_A),
                -sin(psi_A), cos(psi_A)), nrow = 2, byrow = TRUE)
  Tmat <- diag(c(1, psi_R))
  t(R) %*% t(Tmat) %*% Tmat %*% R
}

ani_dist <- function(h, B) {
  sqrt(drop(t(h) %*% B %*% h))
}

# Spherical covariance model.
cov_spherical <- function(r, sigma2 = 3, phi = 6, eta = 2) {
  if (r == 0) return(sigma2 + eta)
  if (r <= phi) return(sigma2 * (1 - 1.5 * (r / phi) + 0.5 * (r / phi)^3))
  0
}

build_cov_matrix <- function(coords, B, sigma2, phi, eta) {
  n <- nrow(coords)
  out <- matrix(0, n, n)
  for (i in seq_len(n)) {
    out[i, i] <- sigma2 + eta
    if (i < n) {
      for (j in (i + 1):n) {
        h <- as.numeric(coords[i, ] - coords[j, ])
        val <- cov_spherical(ani_dist(h, B), sigma2 = sigma2, phi = phi, eta = eta)
        out[i, j] <- val
        out[j, i] <- val
      }
    }
  }
  out
}


simulate_nongaussian_field <- function(SC, n1, n2, alpha = alpha_ng) {
  # Step 1: Gaussian spherical field
  eps <- rnorm(n1 * n2)
  g <- as.numeric(t(SC) %*% eps)
  Gmat <- matrix(g, nrow = n1, ncol = n2)

  # Step 2: quadratic non-Gaussian transform
  G2_centered <- Gmat^2 - mean(Gmat^2)
  Xmat <- Gmat + alpha * G2_centered

  # Step 3: center and standardize for comparability
  Xmat <- Xmat - mean(Xmat)
  Xmat <- Xmat / sd(c(Xmat))

  Xmat
}

compute_variogram <- function(Zmat, lagx, lagy) {
  nr <- nrow(Zmat)
  nc <- ncol(Zmat)
  x_to <- nr - lagx
  y_to <- nc - lagy
  dif <- Zmat[1:x_to, 1:y_to, drop = FALSE] -
    Zmat[(1 + lagx):nr, (1 + lagy):nc, drop = FALSE]
  mean(dif^2)
}

compute_vario_diff <- function(Zmat) {
  gamma_hi <- compute_variogram(Zmat, h_i[1], h_i[2])
  gamma_hj <- compute_variogram(Zmat, h_j[1], h_j[2])
  gamma_hj - gamma_hi
}

periodogram_generic <- function(Z, locs, omega) {
  Z <- as.numeric(Z)
  locs <- as.matrix(locs)
  omega <- as.numeric(omega)
  phase <- locs %*% omega
  dft <- sum(Z * exp(-1i * phase))
  (2 * pi)^(-2) * (1 / nrow(locs)) * Mod(dft)^2
}

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
  if (!all(dim(pMatrix) == c(nr, nc))) stop("pMatrix dimensions do not match weight object.")
  out <- matrix(0, nrow = nr, ncol = nc)
  for (r0 in seq_len(nr)) {
    for (c0 in seq_len(nc)) {
      num <- sum(W[r0, c0, , ] * pMatrix)
      den <- D[r0, c0]
      out[r0, c0] <- if (den > 0) num / den else 0
    }
  }
  if (!is.null(zero_row) && !is.null(zero_col)) out[zero_row, zero_col] <- 0
  out
}

centered_range <- function(n) {
  if (n %% 2 == 0) seq(-n / 2, n / 2 - 1) else seq(-(n - 1) / 2, (n - 1) / 2)
}

coords <- as.matrix(expand.grid(x = xSeq, y = ySeq, KEEP.OUT.ATTRS = FALSE))
ProInd <- coords
T_count_full  <- (n1 - 1) * n2
T_count_block <- (nBlockx - 1) * nBlocky

j1_full <- centered_range(n1)
j2_full <- centered_range(n2)
f1 <- (2 * pi * j1_full) / n1
f2 <- (2 * pi * j2_full) / n2
freqGridData_full <- as.matrix(expand.grid(f1 = f1, f2 = f2, KEEP.OUT.ATTRS = FALSE))
zero_full <- which(freqGridData_full[, 1] == 0 & freqGridData_full[, 2] == 0)
use_full  <- setdiff(seq_len(nrow(freqGridData_full)), zero_full)
freqGridData <- freqGridData_full[use_full, , drop = FALSE]

psi_full <- 2 * cos(as.vector(freqGridData %*% h_j)) - 2 * cos(as.vector(freqGridData %*% h_i))
cos_pair_full <- psi_full * (psi_full + psi_full)

Jn <- expand.grid(j1 = j1_full, j2 = j2_full)
Jn$row <- max(j2_full) - Jn$j2 + 1
Jn$col <- Jn$j1 - min(j1_full) + 1
zero_ind <- which(Jn$j1 == 0 & Jn$j2 == 0)
use_ind <- setdiff(seq_len(nrow(Jn)), zero_ind)

self_sym <- which((Jn$j1 == 0 & Jn$j2 == 0) |
                    ((n1 %% 2 == 0) & abs(Jn$j1) == n1 / 2) |
                    ((n2 %% 2 == 0) & abs(Jn$j2) == n2 / 2))
pairable_ind <- setdiff(use_ind, self_sym)
pos_ind <- pairable_ind[Jn$j1[pairable_ind] > 0 |
                          (Jn$j1[pairable_ind] == 0 & Jn$j2[pairable_ind] > 0)]
neg_ind <- match(paste(-Jn$j1[pos_ind], -Jn$j2[pos_ind]), paste(Jn$j1, Jn$j2))
if (anyNA(neg_ind)) stop("Failed to match some negative-frequency partners in full grid.")
resid_ind <- setdiff(use_ind, c(pos_ind, neg_ind))

weight_obj <- build_literal_kernel_weights(
  f1 = f1,
  f2 = f2,
  hx = hx,
  hy = hy,
  kernel = kernel_type,
  periodic = TRUE
)

j1_sub <- centered_range(nBlockx)
j2_sub <- centered_range(nBlocky)
f1_sub <- (2 * pi * j1_sub) / nBlockx
f2_sub <- (2 * pi * j2_sub) / nBlocky
freqGrid_sub_full <- as.matrix(expand.grid(f1 = f1_sub, f2 = f2_sub, KEEP.OUT.ATTRS = FALSE))
zero_sub <- which(freqGrid_sub_full[, 1] == 0 & freqGrid_sub_full[, 2] == 0)
use_sub  <- setdiff(seq_len(nrow(freqGrid_sub_full)), zero_sub)
freqGrid_sub <- freqGrid_sub_full[use_sub, , drop = FALSE]

psi_sub <- 2 * cos(as.vector(freqGrid_sub %*% h_j)) - 2 * cos(as.vector(freqGrid_sub %*% h_i))
cos_pair_sub <- psi_sub * (psi_sub + psi_sub)

xySubInd <- as.matrix(expand.grid(
  1:(n1 - nBlockx + 1),
  1:(n2 - nBlocky + 1),
  KEEP.OUT.ATTRS = FALSE
))
count <- nrow(xySubInd)

grid <- list(x = xSeq, y = ySeq)
ProcessIndices <- coords

B <- make_aniso_B(psi_A = psi_A, psi_R = psiR)
cov_matrix <- build_cov_matrix(coords, B, sigma2 = sigma2_par,
                               phi = phi_par, eta = eta_par)
cov_matrix <- cov_matrix + diag(1e-10, nrow(cov_matrix))
SC <- chol(cov_matrix)

n_cores <- get_ncores()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
doRNG::registerDoRNG(123)

mc_results <- foreach(
  q = seq_len(nMonteCarlo),
  .packages = c("fields", "matrixStats"),
  .export = c(
    "SC", "n1", "n2", "N", "nBoot", "grid", "ProcessIndices", "ProInd",
    "Delta", "xSeq", "ySeq", "T_count_full", "T_count_block",
    "freqGridData", "freqGrid_sub", "psi_full", "psi_sub",
    "cos_pair_full", "cos_pair_sub", "Jn", "use_ind", "zero_ind",
    "pos_ind", "neg_ind", "resid_ind", "weight_obj", "xySubInd", "count",
    "nBlockx", "nBlocky", "b", "compute_vario_diff", "compute_variogram","periodogram_generic",
    "smooth_periodogram_fast", "simulate_nongaussian_field", "alpha_ng"
  )
) %dorng% {
  Zmat <- simulate_nongaussian_field(SC, n1 = n1, n2 = n2)
  Zvec <- as.vector(Zmat)

  d_full <- compute_vario_diff(Zmat)
  T_obs_spat <- T_count_full * d_full^2

  raw_pd <- apply(freqGridData, 1, function(omega) periodogram_generic(Zvec, ProInd, omega))
  mHatData <- ((2 * pi)^2) * (1 / N) * sum(psi_full * raw_pd)
  T_obs_freq <- abs(sqrt(N) * mHatData)

  Gprocess <- as.surface(grid, Zvec)

  sub_results <- vector("list", count)
  for (k in seq_len(count)) {
    xSub <- xySubInd[k, 1] + 0:(nBlockx - 1)
    ySub <- xySubInd[k, 2] + 0:(nBlocky - 1)

    keep <- ProcessIndices[, 1] %in% xSub & ProcessIndices[, 2] %in% ySub
    xySubProcess <- ProcessIndices[keep, , drop = FALSE]

    xInd <- ((xySubProcess[, 1] - xSeq[1]) / Delta) + 1
    yInd <- ((xySubProcess[, 2] - ySeq[1]) / Delta) + 1
    GprocessSub <- Gprocess$z[cbind(xInd, yInd)]
    Zblock <- matrix(GprocessSub, nrow = nBlockx, ncol = nBlocky)

    d_block <- compute_vario_diff(Zblock)
    T_sub_block <- T_count_block * d_block^2

    pdgramVals_sub <- apply(
      freqGrid_sub, 1,
      function(omega) periodogram_generic(GprocessSub, xySubProcess, omega)
    )

    mHatSub_k <- ((2 * pi)^2) * (1 / b) * sum(psi_sub * pdgramVals_sub)

    sub_results[[k]] <- list(
      mHatSub_k = mHatSub_k,
      pdgramVals = pdgramVals_sub,
      T_sub_block = T_sub_block
    )
  }

  mHatSub <- sapply(sub_results, function(x) x$mHatSub_k)
  gEst_q <- mean(mHatSub)
  Hsub <- sqrt(b) * (mHatSub - gEst_q)
  mHatVar_q <- (count - 1) / count * var(Hsub)
  bias_q <- sqrt(b) * mean(mHatSub - mHatData)

  allPdgram <- do.call(cbind, lapply(sub_results, function(x) x$pdgramVals))
  pd_var <- apply(allPdgram, 1, var)
  subSigma <- cos_pair_sub * ((count - 1) / count) * pd_var
  sigma1_q <- ((2 * pi)^4) * (1 / b) * sum(subSigma)

  T_sub_vals <- sapply(sub_results, function(x) x$T_sub_block)
  p_sub_q <- mean(T_sub_vals >= T_obs_spat)

  sigma2_q <- mHatVar_q - sigma1_q
  if (!is.finite(sigma2_q) || sigma2_q < 0) sigma2_q <- max(0, sigma2_q, na.rm = TRUE)

  pMatrix <- matrix(0, nrow = n2, ncol = n1)
  pMatrix[cbind(Jn$row[use_ind], Jn$col[use_ind])] <- raw_pd
  zr <- Jn$row[zero_ind]
  zc <- Jn$col[zero_ind]
  pMatrix[zr, zc] <- 0

  smMatrix <- smooth_periodogram_fast(
    pMatrix = pMatrix,
    weight_obj = weight_obj,
    zero_row = zr,
    zero_col = zc
  )
  smoothPdgram <- smMatrix[cbind(Jn$row[use_ind], Jn$col[use_ind])]

  v_hpb <- numeric(nBoot)
  for (l in seq_len(nBoot)) {
    expVals <- numeric(nrow(Jn))
    if (length(pos_ind) > 0) {
      epos <- rexp(length(pos_ind), rate = 1) - 1
      expVals[pos_ind] <- epos
      expVals[neg_ind] <- epos
    }
    if (length(resid_ind) > 0) {
      expVals[resid_ind] <- rexp(length(resid_ind), rate = 1) - 1
    }
    expVector <- expVals[use_ind]
    v_hpb[l] <- sqrt(N) * ((2 * pi)^2) * (1 / N) * sum(psi_full * smoothPdgram * expVector)
  }

  var_vHpb_q <- ((2 * pi)^4) * (1 / N) * sum(cos_pair_full * smoothPdgram^2)

  v_hpb_centered <- v_hpb - mean(v_hpb)
  p_fdb_q <- mean(abs(v_hpb_centered) >= T_obs_freq)

  vCorrect_q <- var_vHpb_q + sigma2_q
  if (!is.finite(var_vHpb_q) || var_vHpb_q <= 0 || !is.finite(vCorrect_q) || vCorrect_q <= 0) {
    p_hfdb_q <- NA_real_
  } else {
    T_hpb <- (sqrt(vCorrect_q) * v_hpb / sqrt(var_vHpb_q)) + bias_q
    T_hpb_centered <- T_hpb - mean(T_hpb)
    p_hfdb_q <- mean(abs(T_hpb_centered) >= T_obs_freq)
  }

  c(
    p_sub = p_sub_q,
    p_fdb = p_fdb_q,
    p_hfdb = p_hfdb_q,
    sigma1 = sigma1_q,
    sigma2 = sigma2_q,
    bias = bias_q,
    gEstimate = mHatData,
    gEst = gEst_q,
    mHatVar = mHatVar_q,
    var_vHpb = var_vHpb_q,
    vCorrect = vCorrect_q,
    T_obs_freq = T_obs_freq
  )
}

stopCluster(cl)

res_mat <- do.call(rbind, mc_results)

pValue_sub  <- as.numeric(res_mat[, "p_sub"])
pValue_fdb  <- as.numeric(res_mat[, "p_fdb"])
pValue_hfdb <- as.numeric(res_mat[, "p_hfdb"])
sigma1      <- as.numeric(res_mat[, "sigma1"])
sigma2_vec  <- as.numeric(res_mat[, "sigma2"])
bias        <- as.numeric(res_mat[, "bias"])
gEstimate   <- as.numeric(res_mat[, "gEstimate"])
gEst        <- as.numeric(res_mat[, "gEst"])
mHatVar     <- as.numeric(res_mat[, "mHatVar"])
var_vHpb    <- as.numeric(res_mat[, "var_vHpb"])
vCorrect    <- as.numeric(res_mat[, "vCorrect"])
T_obs_freq  <- as.numeric(res_mat[, "T_obs_freq"])

summary_out <- data.frame(
  psi_R = psiR,
  method = c("Spatial Subsampling", "FDWB", "HFDB"),
  reject_10 = c(
    mean(pValue_sub  < alpha_levels[1], na.rm = TRUE),
    mean(pValue_fdb  < alpha_levels[1], na.rm = TRUE),
    mean(pValue_hfdb < alpha_levels[1], na.rm = TRUE)
  ),
  reject_05 = c(
    mean(pValue_sub  < alpha_levels[2], na.rm = TRUE),
    mean(pValue_fdb  < alpha_levels[2], na.rm = TRUE),
    mean(pValue_hfdb < alpha_levels[2], na.rm = TRUE)
  ),
  mean_pvalue = c(
    mean(pValue_sub, na.rm = TRUE),
    mean(pValue_fdb, na.rm = TRUE),
    mean(pValue_hfdb, na.rm = TRUE)
  ),
  mean_sigma1 = c(mean(sigma1, na.rm = TRUE), mean(sigma1, na.rm = TRUE), mean(sigma1, na.rm = TRUE)),
  mean_sigma2 = c(mean(sigma2_vec, na.rm = TRUE), mean(sigma2_vec, na.rm = TRUE), mean(sigma2_vec, na.rm = TRUE)),
  mean_bias = c(mean(bias, na.rm = TRUE), mean(bias, na.rm = TRUE), mean(bias, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

result_obj <- list(
  summary = summary_out,
  pValue_sub = pValue_sub,
  pValue_fdb = pValue_fdb,
  pValue_hfdb = pValue_hfdb,
  sigma1 = sigma1,
  sigma2 = sigma2_vec,
  bias = bias,
  gEstimate = gEstimate,
  gEst = gEst,
  mHatVar = mHatVar,
  var_vHpb = var_vHpb,
  vCorrect = vCorrect,
  T_obs_freq = T_obs_freq,
  params = list(
    psiR = psiR,
    nMonteCarlo = nMonteCarlo,
    nBoot = nBoot,
    Delta = Delta,
    n1 = n1,
    n2 = n2,
    betaValue = betaValue,
    hx = hx,
    hy = hy,
    kernel_type = kernel_type,
    sigma2_par = sigma2_par,
    phi_par = phi_par,
    eta_par = eta_par,
    psi_A = psi_A,
    alpha_ng = alpha_ng,
    h_i = h_i,
    h_j = h_j,
    seed = 2
  )
)

write.csv(summary_out, file.path(out_dir, "summary.csv"), row.names = FALSE)
saveRDS(result_obj, file.path(out_dir, "results.rds"), compress = FALSE)

elapsed <- proc.time() - start_time
run_info <- c(
  paste0("psiR=", psiR),
  paste0("nMonteCarlo=", nMonteCarlo),
  paste0("nBoot=", nBoot),
  paste0("seed=2"),
  paste0("ncores=", n_cores),
  paste0("elapsed_seconds=", round(unname(elapsed["elapsed"]), 3))
)
writeLines(run_info, file.path(out_dir, "run_info.txt"))

print(summary_out)
cat("\nElapsed time:\n")
print(elapsed)

