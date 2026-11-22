# =============================================================================
# 04_sim_generator.R
# -----------------------------------------------------------------------------
# Simulated multi-block data with PLANTED block-pair structure.
#
# Model: each block X_b = sum_f  f * w_{f,b}^T  + E,   E ~ N(0, 1) per variable.
#   - A factor f is either CLASS-related (f = rho*y + sqrt(1-rho^2)*e, y = +-1)
#     or a NUISANCE (f ~ N(0,1), unrelated to the class).
#   - w_{f,b} is sparse: k_support variables, entries +-s_{f,b}, so s is the
#     per-variable signal SD on the support (s^2 = per-variable SNR there).
#   - Loadings are drawn ONCE per dataset and shared by its train and test sets,
#     so test BER is a true out-of-sample error (no CV needed for evaluation).
#
# Planted design ("oracle"): a block pair is coupled iff both blocks load on a
# shared CLASS factor. Pairs linked only through a nuisance factor should be
# decoupled: coupling them pulls components toward class-irrelevant variance.
# =============================================================================

## ---- Scenario specifications -------------------------------------------------

# 3 blocks (d = 3). Class signal strong in B, weak in A; strong nuisance in A, C.
# Planted optimum: A-B high, A-C low, B-C low  -> OFF the 1-D diagonal.
spec_3block <- function(sA = 0.3, sB = 0.6, sU = 1.0, p = 200, k = 20, rho = 0.8)
  list(name = sprintf("3blk_sA%.2f_sB%.2f_sU%.2f", sA, sB, sU),
       blocks = c(A = p, B = p, C = p), k_support = k, rho = rho,
       factors = list(
         class = list(class = TRUE,  load = c(A = sA, B = sB, C = 0)),
         nuis  = list(class = FALSE, load = c(A = sU, B = 0,  C = sU))))

# Control: every block carries the same class factor, no nuisance. All pairs are
# equivalent, so the planted design is the full design: ON the diagonal. Here
# d-dimensional search must NOT beat the sweep; if it does, that is selection bias.
spec_3block_null <- function(s = 0.3, p = 200, k = 20, rho = 0.8)
  list(name = sprintf("3blk_control_s%.2f", s),
       blocks = c(A = p, B = p, C = p), k_support = k, rho = rho,
       factors = list(
         class = list(class = TRUE, load = c(A = s, B = s, C = s))))

# 5 blocks (d = 10). Class factor in A (weak), B (strong), D (weak);
# nuisance 1 in A, C; nuisance 2 in D, E.
# Planted optimum: A-B, A-D, B-D high; all other 7 pairs low.
spec_5block <- function(sW = 0.3, sB = 0.6, sU = 1.0, p = 150, k = 20, rho = 0.8)
  list(name = sprintf("5blk_sW%.2f_sB%.2f_sU%.2f", sW, sB, sU),
       blocks = c(A = p, B = p, C = p, D = p, E = p), k_support = k, rho = rho,
       factors = list(
         class = list(class = TRUE,  load = c(A = sW, B = sB, C = 0,  D = sW, E = 0)),
         nuis1 = list(class = FALSE, load = c(A = sU, B = 0,  C = sU, D = 0,  E = 0)),
         nuis2 = list(class = FALSE, load = c(A = 0,  B = 0,  C = 0,  D = sU, E = sU))))

## ---- Loadings (drawn once per dataset) ---------------------------------------
draw_loadings <- function(spec) {
  W <- lapply(spec$factors, function(f) {
    out <- lapply(names(spec$blocks), function(b) {
      p <- spec$blocks[[b]]; s <- f$load[[b]]; w <- numeric(p)
      if (s > 0) {
        idx <- sample(p, spec$k_support)
        w[idx] <- s * sample(c(-1, 1), spec$k_support, replace = TRUE)
      }
      w
    })
    names(out) <- names(spec$blocks); out
  })
  names(W) <- names(spec$factors); W
}

## ---- One sample set given loadings -------------------------------------------
gen_set <- function(spec, W, n) {
  y  <- rep(c(-1, 1), length.out = n)[sample(n)]          # balanced 2 classes
  Fm <- sapply(names(spec$factors), function(fn) {
    if (spec$factors[[fn]]$class) spec$rho * y + sqrt(1 - spec$rho^2) * rnorm(n)
    else rnorm(n)
  })
  Fm <- matrix(Fm, nrow = n, dimnames = list(NULL, names(spec$factors)))
  X <- lapply(names(spec$blocks), function(b) {
    p <- spec$blocks[[b]]
    M <- matrix(rnorm(n * p), n, p)
    for (fn in names(spec$factors)) M <- M + outer(Fm[, fn], W[[fn]][[b]])
    colnames(M) <- paste0(b, "_", seq_len(p))
    M
  })
  names(X) <- names(spec$blocks)
  list(X = X, Y = factor(ifelse(y > 0, "pos", "neg")))
}

## ---- A full dataset: train + large independent test --------------------------
simulate_dataset <- function(spec, n_train = 100, n_test = 2000) {
  W <- draw_loadings(spec)
  list(train = gen_set(spec, W, n_train), test = gen_set(spec, W, n_test),
       spec = spec)
}

## ---- Planted design ----------------------------------------------------------
# Logical vector over block pairs, in the SAME order as M[upper.tri(M)].
pair_index <- function(bn) {
  idx <- which(upper.tri(diag(length(bn))), arr.ind = TRUE)
  idx <- idx[order(idx[, 2], idx[, 1]), , drop = FALSE]    # column-major, as upper.tri
  list(idx = idx, names = paste(bn[idx[, 1]], bn[idx[, 2]], sep = "-"))
}

oracle_pairs <- function(spec) {
  bn <- names(spec$blocks); pi <- pair_index(bn)
  hi <- apply(pi$idx, 1, function(ij) any(vapply(spec$factors, function(f)
    f$class && f$load[[bn[ij[1]]]] > 0 && f$load[[bn[ij[2]]]] > 0, logical(1))))
  setNames(hi, pi$names)
}

# Coordinate vector (logit space) for a logical pair pattern.
pattern_coord <- function(hi, L = LOGIT_RANGE) ifelse(hi, L, -L)
