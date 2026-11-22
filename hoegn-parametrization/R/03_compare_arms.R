# =============================================================================
# 03_compare_arms.R
# -----------------------------------------------------------------------------
# Does d-dimensional search beat the 1-D sweep? Nested-CV comparison of search
# arms under BER (the demonstration objective).
#
# Outer CV: each arm SELECTS a design on the outer training fold (by inner-CV
# BER), then the selected design is scored ONCE on the outer test fold.
# All arms share the same outer folds (paired) and, within a training fold,
# the same inner folds (ber_objective reseeds identically for every candidate).
#
# Arms
#   sweep        1-D uniform sweep c*(J-I), N_SWEEP points      (the restriction)
#   search       maxDissim in logit space, N_SEARCH + 3 refs    (the method)
#   search_unit  same, sampled uniformly in [0,1]^d             (density check)
#   grid         exhaustive grid over GRID_LEVELS^d             (ground truth, d<=3)
#   empirical    fixed mean-|cor| design, no selection
#   null / full  fixed box-edge designs (theta = -L / +L)
#
# Reads : public_datasets.rds          (00_data_loader.R)
# Writes: arms_<dataset>.csv           one row per repeat x fold x arm
#         arms_<dataset>_landscape.rds inner-BER of every grid point per fold
#
# Run from the folder that contains R/ and public_datasets.rds.
# Set QUICK <- TRUE first: a smoke test that finishes in a few minutes.
# =============================================================================

source("R/02_ber_objective.R")
suppressPackageStartupMessages(library(parallel))

## ---- CONFIG -----------------------------------------------------------------
DATASET       <- "nutrimouse"    # "tcga" | "nutrimouse" | "mogonet_brca"
NCOMP_RUN     <- 2                 # registered value is 2; 4 = sensitivity analysis
QUICK         <- FALSE             # TRUE = tiny settings for a smoke test
OUTER_FOLDS   <- 5
OUTER_REPEATS <- 5                 # repeats 1-3 use the SAME outer folds as the ncomp = 2 run
INNER_FOLDS   <- 5
N_SWEEP       <- 15
N_SEARCH      <- 60
GRID_LEVELS   <- c(0.0025, 0.1, 0.25, 0.5, 0.75, 0.9975)   # 6^3 = 216
INCLUDE_UNIT  <- TRUE             # density check already done in the ncomp = 2 run
FLAT_TOL      <- 0.02              # BER tolerance for the flatness diagnostic
NCORES        <- max(1L, detectCores() - 1L)
SEED          <- 42

if (QUICK) {
  OUTER_FOLDS <- 3; INNER_FOLDS <- 3; N_SWEEP <- 4; N_SEARCH <- 6
  GRID_LEVELS <- c(0.1, 0.9)
}

# Number of components: overrides the global NCOMP from 01_parametrization_core.R.
# fit_block(), make_keepX() and get_class() all read NCOMP at call time.
NCOMP <- NCOMP_RUN
TAG   <- paste0(DATASET, if (NCOMP != 2) paste0("_ncomp", NCOMP) else "",
                if (QUICK) "_quick" else "")
PARTIAL_CSV <- sprintf("arms_%s_partial.csv", TAG)
PARTIAL_RDS <- sprintf("arms_%s_partial_landscape.rds", TAG)

## ---- Data -------------------------------------------------------------------
D  <- readRDS("public_datasets.rds")[[DATASET]]
X  <- D$X; Y <- D$Y; bn <- names(X)
nb <- length(X); d <- nb * (nb - 1) / 2
L  <- LOGIT_RANGE
coord_emp <- design_to_coord(D$design)

idx <- which(upper.tri(diag(nb)), arr.ind = TRUE)          # same order as M[upper.tri(M)]
PAIR_NAMES <- paste(bn[idx[, 1]], bn[idx[, 2]], sep = "-")

grid_coords <- as.matrix(expand.grid(rep(list(lg(GRID_LEVELS)), d)))
dimnames(grid_coords) <- NULL

## ---- Convergence tally (fork-safe) ------------------------------------------
# fit_block() counts into .CONV, but increments inside forked workers are lost.
# Every call below returns its own delta, and we sum deltas here instead.
TALLY <- c(nonconv = 0, total = 0)
with_delta <- function(expr_fun) {
  n0 <- .CONV$nonconv; t0 <- .CONV$total
  v  <- expr_fun()
  c(value = v, nonconv = .CONV$nonconv - n0, total = .CONV$total - t0)
}

## ---- Parallel scoring of a candidate set ------------------------------------
score_par <- function(coords, Xd, Yd) {
  out <- mclapply(seq_len(nrow(coords)), function(i)
    tryCatch(with_delta(function()
      ber_objective(coord_to_design(coords[i, ], bn), Xd, Yd,
                    inner_folds = INNER_FOLDS)),
      error = function(e) c(value = NA, nonconv = 0, total = 0)),
    mc.cores = NCORES)
  M <- do.call(rbind, lapply(out, function(r)
    if (inherits(r, "try-error") || length(r) != 3) c(NA, 0, 0) else r))
  TALLY <<- TALLY + c(sum(M[, 2]), sum(M[, 3]))
  list(val = M[, 1], fits = sum(M[, 3]))
}

# Ties go to the FIRST candidate in order. On a flat landscape this decides.
select_from <- function(coords, Xtr, Ytr) {
  s <- score_par(coords, Xtr, Ytr)
  j <- which.min(replace(s$val, is.na(s$val), Inf))
  list(coord = coords[j, ], inner = s$val[j], fits = s$fits, all = s$val)
}

outer_score <- function(coord, Xtr, Ytr, Xte, Yte) {
  r <- with_delta(function()
    ber_outer_eval(coord_to_design(coord, bn), Xtr, Ytr, Xte, Yte))
  TALLY <<- TALLY + r[c("nonconv", "total")]
  unname(r["value"])
}

## ---- Cost estimate ----------------------------------------------------------
n_sel <- N_SWEEP + (N_SEARCH + 3) * (1 + INCLUDE_UNIT) + nrow(grid_coords)
est   <- OUTER_REPEATS * OUTER_FOLDS * (n_sel * INNER_FOLDS + 7)
cat(sprintf("%s: %d blocks, d = %d | ncomp %d | grid %d points | cores %d\n",
            DATASET, nb, d, NCOMP, nrow(grid_coords), NCORES))
cat(sprintf("Estimated model fits: %d\n\n", est))

## ---- Main loop --------------------------------------------------------------
# Resume: folds already saved in the partial file are skipped. Outer folds and
# candidate draws are seeded per (repeat, fold), so a resumed run is identical.
rows <- list(); landscape <- list(); done <- character(0)
if (file.exists(PARTIAL_CSV)) {
  prev <- read.csv(PARTIAL_CSV, check.names = FALSE)
  rows <- list(prev); done <- unique(paste(prev$rep, prev$fold))
  if (file.exists(PARTIAL_RDS)) landscape <- readRDS(PARTIAL_RDS)
  cat(sprintf("Resuming: %d fold(s) already done, skipping them.\n\n", length(done)))
}

for (r in seq_len(OUTER_REPEATS)) {
  set.seed(SEED + r)
  folds <- caret::createFolds(Y, k = OUTER_FOLDS, list = TRUE, returnTrain = FALSE)

  for (f in seq_along(folds)) {
    if (paste(r, f) %in% done) next
    t_start <- Sys.time()
    te  <- folds[[f]]; tr <- setdiff(seq_along(Y), te)
    Xtr <- lapply(X, function(b) b[tr, , drop = FALSE])
    Xte <- lapply(X, function(b) b[te, , drop = FALSE])
    Ytr <- droplevels(Y[tr]); Yte <- Y[te]
    if (nlevels(Ytr) < nlevels(Y)) { cat("  skipping fold (missing class)\n"); next }

    set.seed(SEED + 1000 * r + f)          # candidate draws reproducible per fold
    arms <- list(
      sweep  = uniform_sweep(d, N_SWEEP, L),
      search = generate_candidates(d, N_SEARCH, coord_emp, space = "logit")$coord,
      grid   = grid_coords)
    if (INCLUDE_UNIT)
      arms$search_unit <- generate_candidates(d, N_SEARCH, coord_emp, space = "unit")$coord

    sel <- lapply(arms, select_from, Xtr = Xtr, Ytr = Ytr)
    landscape[[sprintf("r%d_f%d", r, f)]] <- sel$grid$all

    fixed <- list(empirical = coord_emp, null = rep(-L, d), full = rep(L, d))
    for (nm in names(fixed))
      sel[[nm]] <- list(coord = fixed[[nm]], inner = NA_real_, fits = 0L)

    for (nm in names(sel)) {
      ob  <- outer_score(sel[[nm]]$coord, Xtr, Ytr, Xte, Yte)
      ent <- il(sel[[nm]]$coord)
      rows[[length(rows) + 1]] <- data.frame(
        rep = r, fold = f, arm = nm,
        t(setNames(round(ent, 4), PAIR_NAMES)),
        inner_ber = sel[[nm]]$inner, outer_ber = ob, n_fits = sel[[nm]]$fits,
        check.names = FALSE)
    }

    g <- sel$grid$all
    cat(sprintf("repeat %d fold %d/%d  (%.1f min) | grid inner BER range %.3f-%.3f, %.0f%% within %.2f of min\n",
                r, f, OUTER_FOLDS, as.numeric(difftime(Sys.time(), t_start, units = "mins")),
                min(g, na.rm = TRUE), max(g, na.rm = TRUE),
                100 * mean(g <= min(g, na.rm = TRUE) + FLAT_TOL, na.rm = TRUE), FLAT_TOL))

    # Save after every fold, so an interruption loses at most one fold.
    write.csv(do.call(rbind, rows), PARTIAL_CSV, row.names = FALSE)
    saveRDS(landscape, PARTIAL_RDS)
  }
}

res <- do.call(rbind, rows)
write.csv(res, sprintf("arms_%s.csv", TAG), row.names = FALSE)
saveRDS(list(landscape = landscape, grid_coords = grid_coords,
             grid_levels = GRID_LEVELS, pair_names = PAIR_NAMES),
        sprintf("arms_%s_landscape.rds", TAG))

## ---- Summary ----------------------------------------------------------------
ORDER <- intersect(c("sweep", "search", "search_unit", "grid",
                     "empirical", "null", "full"), unique(res$arm))

cat("\n== OUTER BER by arm (lower is better) ==\n")
tb <- do.call(rbind, lapply(ORDER, function(a) {
  x <- res$outer_ber[res$arm == a]
  data.frame(arm = a, mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
             fits_per_fold = mean(res$n_fits[res$arm == a]))
}))
print(format(tb, digits = 3), row.names = FALSE)

res$key <- paste(res$rep, res$fold)
w <- reshape(res[, c("key", "arm", "outer_ber")], idvar = "key",
             timevar = "arm", direction = "wide")
for (ref in c("sweep", "empirical")) {
  cat(sprintf("\n== Paired vs %s (negative = better than %s) ==\n", toupper(ref), ref))
  for (a in setdiff(ORDER, ref)) {
    dif <- w[[paste0("outer_ber.", a)]] - w[[paste0("outer_ber.", ref)]]
    cat(sprintf("  %-12s %+.4f   better in %d / %d folds, tied in %d\n",
                a, mean(dif, na.rm = TRUE), sum(dif < 0, na.rm = TRUE),
                sum(!is.na(dif)), sum(dif == 0, na.rm = TRUE)))
  }
}

cat("\n== Selection stability (SD of selected entries across folds, [0,1] scale) ==\n")
for (a in c("sweep", "search", "search_unit", "grid")) {
  if (!a %in% res$arm) next
  s <- apply(res[res$arm == a, PAIR_NAMES, drop = FALSE], 2, sd)
  cat(sprintf("  %-12s %s\n", a, paste(sprintf("%s %.3f", PAIR_NAMES, s), collapse = " | ")))
}

cat("\n== Grid landscape flatness (inner BER, per fold) ==\n")
for (k in names(landscape)) {
  g <- landscape[[k]]
  cat(sprintf("  %s  min %.3f  sd %.3f  share within %.2f of min: %.0f%%\n",
              k, min(g, na.rm = TRUE), sd(g, na.rm = TRUE), FLAT_TOL,
              100 * mean(g <= min(g, na.rm = TRUE) + FLAT_TOL, na.rm = TRUE)))
}
cat("-> a large share near the minimum means the objective barely distinguishes designs.\n")

cat(sprintf("\nSGCCA non-converged fits: %d / %d (%.2f%%)\n",
            TALLY["nonconv"], TALLY["total"],
            100 * TALLY["nonconv"] / max(1, TALLY["total"])))
cat(sprintf("Written: arms_%s.csv, arms_%s_landscape.rds\n", TAG, TAG))

## ---- Compare with the registered ncomp = 2 run on the SAME outer folds ------
ref_file <- sprintf("arms_%s.csv", DATASET)
if (NCOMP != 2 && file.exists(ref_file)) {
  ref <- read.csv(ref_file, check.names = FALSE)
  ref <- ref[paste(ref$rep, ref$fold) %in% paste(res$rep, res$fold), ]
  cat(sprintf("\n== ncomp %d vs ncomp 2 (same %d outer folds): mean outer BER ==\n",
              NCOMP, length(unique(paste(res$rep, res$fold)))))
  for (a in intersect(ORDER, unique(ref$arm))) {
    cat(sprintf("  %-12s ncomp2 %.4f  ->  ncomp%d %.4f\n", a,
                mean(ref$outer_ber[ref$arm == a], na.rm = TRUE), NCOMP,
                mean(res$outer_ber[res$arm == a], na.rm = TRUE)))
  }
}
