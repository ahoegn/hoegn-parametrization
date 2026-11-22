# =============================================================================
# 01_parametrization_core.R
# -----------------------------------------------------------------------------
# Extracted from 1_public_hoegn_parametrization_v5.R. Contains ONLY the parts
# that belong to the parametrization paper. Everything criterion-related
# (coranking, q_local/q_global, envelope/maxmin, legacy metrics, Procrustes,
# k_split) has been REMOVED (separate project).
#
# The objective is now PLUGGABLE: any function  objective(design, Xd, Yd) -> number
# (lower = better by convention here, e.g. nested-CV BER from 02_ber_objective.R).
#
# THE METHOD
#   A block.splsda design matrix is symmetric, hollow, bounded in [0,1]. Its only
#   free parameters are the d = nb*(nb-1)/2 off-diagonal entries. The Hoegn
#   Parametrization maps those onto unconstrained R^d via elementwise logit.
#   Scope: X-X sub-block only. mixOmics adds a Y row/column at weight 1.
# =============================================================================

suppressPackageStartupMessages({ library(mixOmics); library(caret) })

## ---- Single-threaded BLAS ----------------------------------------------------
# The parallel scripts run many small fits at once. A multithreaded BLAS
# (OpenBLAS/MKL) inside every worker oversubscribes the CPU and slows everything
# down. Environment variables set here come too late for a BLAS that R loaded at
# startup, so the thread count is set at runtime via RhpcBLASctl when available.
# Forked workers inherit the setting. (Launching with OPENBLAS_NUM_THREADS=1, or
# putting that line in ~/.Renviron, works too.)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
} else {
  message("Note: BLAS threads not limited. Run install.packages('RhpcBLASctl'), ",
          "or launch with OPENBLAS_NUM_THREADS=1.")
}

## ---- Fixed model settings (never tuned) -------------------------------------
NCOMP       <- 2
KEEPX_FRAC  <- 0.5
LOGIT_RANGE <- 6     # box half-width L: entries in ~[0.0025, 0.9975]. A parameter
                     # of the procedure; report it. Uniform sampling on [-6, 6]
                     # puts only ~37% of draws in [0.1, 0.9].

make_keepX <- function(X, ncomp = NCOMP, frac = KEEPX_FRAC) {
  k <- lapply(X, function(m) rep(max(1L, round(ncol(m) * frac)), ncomp))
  names(k) <- names(X); k
}

## ---- Coordinate <-> design matrix -------------------------------------------
lg <- function(x) { x <- pmin(pmax(x, 1e-6), 1 - 1e-6); log(x / (1 - x)) }
il <- function(x) 1 / (1 + exp(-x))

design_to_coord <- function(M) lg(M[upper.tri(M)])   # diagonal is NOT a parameter

coord_to_design <- function(v, bn) {
  n <- length(bn)
  M <- matrix(0, n, n, dimnames = list(bn, bn))
  M[upper.tri(M)] <- il(v)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  diag(M) <- 0
  M
}

## ---- Fitting with SGCCA convergence counter ---------------------------------
# block.splsda warns rather than errors on non-convergence. Counted and reported.
.CONV <- new.env(parent = emptyenv()); .CONV$nonconv <- 0L; .CONV$total <- 0L

fit_block <- function(Xd, Yd, des, keepX = make_keepX(Xd)) {
  .CONV$total <- .CONV$total + 1L
  withCallingHandlers(
    tryCatch(suppressMessages(
      block.splsda(Xd, Yd, ncomp = NCOMP, keepX = keepX,
                   design = des, scale = TRUE)),
      error = function(e) NULL),
    warning = function(w) {
      if (grepl("converg", conditionMessage(w), ignore.case = TRUE)) {
        .CONV$nonconv <- .CONV$nonconv + 1L
        invokeRestart("muffleWarning")
      }
    })
}

report_convergence <- function() {
  cat(sprintf("SGCCA non-converged fits: %d / %d (%.1f%%)\n",
              .CONV$nonconv, .CONV$total,
              100 * .CONV$nonconv / max(1L, .CONV$total)))
}

## ---- Candidate generation ----------------------------------------------------
# Space-filling (maximum-dissimilarity) sampling in the logit box, plus named
# reference points. NOTE: this is maxDissim, not LHS; earlier code called it LHS.
# space = "logit" samples uniformly in [-L, L]^d; space = "unit" samples
# uniformly in [0,1]^d and maps through logit.
generate_candidates <- function(d, n_cand = 150, coord_emp = NULL,
                                L = LOGIT_RANGE, space = c("logit", "unit")) {
  space <- match.arg(space)
  draw <- function(m) {
    if (space == "logit") matrix(runif(m * d, -L, L), ncol = d)
    else matrix(pmin(pmax(lg(runif(m * d)), -L), L), ncol = d)
  }
  pool  <- draw(2000)
  idx   <- caret::maxDissim(draw(200), pool, n = n_cand)
  cand  <- pool[idx, , drop = FALSE]
  tag   <- rep("search", n_cand)
  refs  <- list(null = rep(-L, d), full = rep(L, d))
  if (!is.null(coord_emp)) refs <- c(list(empirical = coord_emp), refs)
  cand  <- rbind(cand, do.call(rbind, refs))
  tag   <- c(tag, names(refs))
  list(coord = cand, tag = tag)
}

## ---- 1-D baseline: the uniform sweep c * J ------------------------------------
# This IS the restriction the paper argues against: all off-diagonals equal.
uniform_sweep <- function(d, n = 50, L = LOGIT_RANGE)
  matrix(rep(seq(-L, L, length.out = n), each = d), ncol = d, byrow = TRUE)

## ---- Diagnostic slice null -> empirical -> full (figure only) ----------------
# Routed through the empirical point: each cell traces 0 -> M_emp[i,j] -> 1.
slice_path <- function(coord_emp, n = 120, L = LOGIT_RANGE, a_ref = 0.5) {
  d <- length(coord_emp)
  alpha <- seq(0, 1, length.out = n)
  sc <- vapply(alpha, function(a) {
    if (a <= a_ref) (1 - a / a_ref) * rep(-L, d) + (a / a_ref) * coord_emp
    else { b <- (a - a_ref) / (1 - a_ref); (1 - b) * coord_emp + b * rep(L, d) }
  }, numeric(d))
  list(alpha = alpha, coord = if (d == 1) matrix(sc, ncol = 1) else t(sc))
}

## ---- Evaluate a set of coordinates under any objective -----------------------
# objective: function(design_matrix, Xd, Yd) -> numeric(1)
score_coords <- function(coords, objective, Xd, Yd, verbose = TRUE) {
  bn <- names(Xd)
  vapply(seq_len(nrow(coords)), function(i) {
    if (verbose && i %% 10 == 0) cat(sprintf("    candidate %d / %d\n", i, nrow(coords)))
    objective(coord_to_design(coords[i, ], bn), Xd, Yd)
  }, numeric(1))
}

## ---- Nested selection (objective-agnostic) ------------------------------------
# Selection happens inside each outer training fold; the outer test fold is
# then used ONLY to evaluate the selected design. Returns per-fold selected
# coordinates (for stability) and outer-fold performance.
# outer_eval: function(design, Xtr, Ytr, Xte, Yte) -> numeric(1)
nested_select <- function(X, Y, objective, outer_eval, d, coord_emp = NULL,
                          outer_folds = 5, n_cand = 150, seed = 42) {
  set.seed(seed)
  folds <- caret::createFolds(Y, k = outer_folds, list = TRUE, returnTrain = FALSE)
  bn <- names(X)
  res <- lapply(seq_along(folds), function(f) {
    te <- folds[[f]]; tr <- setdiff(seq_along(Y), te)
    Xtr <- lapply(X, function(b) b[tr, , drop = FALSE])
    Xte <- lapply(X, function(b) b[te, , drop = FALSE])
    Ytr <- droplevels(Y[tr]); Yte <- Y[te]
    if (nlevels(Ytr) < nlevels(Y)) return(NULL)
    cand <- generate_candidates(d, n_cand, coord_emp)
    sc   <- score_coords(cand$coord, objective, Xtr, Ytr, verbose = FALSE)
    # Ties go to the FIRST candidate in order (search points, then empirical,
    # null, full). On flat objectives (e.g. BER at its floor) this ordering
    # decides the selection; report it.
    j    <- which.min(replace(sc, is.na(sc), Inf))
    des  <- coord_to_design(cand$coord[j, ], bn)
    list(coord = cand$coord[j, ], inner = sc[j],
         outer = outer_eval(des, Xtr, Ytr, Xte, Yte))
  })
  Filter(Negate(is.null), res)
}

# Per-fold coordinate SD in logit space: low SD = stable selected design.
selection_stability <- function(res) {
  Cm <- do.call(rbind, lapply(res, `[[`, "coord"))
  list(median = apply(Cm, 2, median), sd = apply(Cm, 2, sd), coords = Cm)
}
# CAUTION: the coordinatewise fold median can land OFF a ridge of near-optimal
# designs. Evaluate the median itself; do not assume it inherits the folds' scores.
