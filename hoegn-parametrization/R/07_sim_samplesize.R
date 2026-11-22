# =============================================================================
# 07_sim_samplesize.R
# -----------------------------------------------------------------------------
# Same nested simulation as 06, across TRAINING SAMPLE SIZES, with REPEATED inner
# CV. Question: does the gain of d-dimensional search over the 1-D sweep grow,
# and does the selection optimism (winner's curse) shrink, as n increases?
#
# Differences from 06:
#   - each scenario is run at n_train = 100, 200, 400
#   - inner CV is repeated INNER_REPEATS times (averaging lowers the noise in
#     each candidate's score, which directly reduces the winner's curse)
#   - each (scenario, n) cell is saved to its own CSV as soon as it finishes;
#     rerunning the script SKIPS cells already on disk, so it can be stopped and
#     resumed. Delete sim_n/ to start over.
#
# Original 06 description:
#
# Per simulated dataset (train n = 100, independent test n = 1000):
#   1. Inner 5-fold CV on the TRAINING set scores every candidate design (BER).
#   2. Each arm picks its design from its own candidates:
#        sweep       best of N_SWEEP uniform designs c*(J-I)       (1-D)
#        search      best of maxDissim candidates in logit space   (d-D)
#                    (+ empirical, null, full as reference candidates)
#        search_1se  among search candidates within 1 SE of the best,
#                    the one with the LEAST total coupling (guards against
#                    the winner's curse seen on TCGA)
#   3. The chosen design is refit on the full training set and scored ONCE on
#      the independent test set -> true out-of-sample BER.
#   Fixed references (no selection): empirical, null, full, oracle (planted).
#
# All arms share the same inner folds within a dataset (paired comparison).
# Parallel over datasets.
#
# Writes: sim_n/<scenario>_n<n_train>.csv, then sim_samplesize.csv (all cells)
# Set QUICK <- TRUE first for a smoke test (~3-5 min).
# Full run: several hours. Run from a terminal:
#   Rscript R/07_sim_samplesize.R 2>&1 | tee sim_samplesize_log.txt
# =============================================================================

source("R/02_ber_objective.R")
source("R/04_sim_generator.R")
suppressPackageStartupMessages(library(parallel))

## ---- CONFIG -----------------------------------------------------------------
QUICK         <- FALSE
N_TRAIN_GRID  <- c(100, 200, 400)
N_TEST        <- 1000
INNER_FOLDS   <- 5
INNER_REPEATS <- 3
N_SWEEP     <- 15
NCORES      <- max(1L, detectCores() - 1L)
# RNGkind("L'Ecuyer-CMRG"); set.seed(7)

# Ordered cheapest first, so the most informative results arrive early.
SCENARIOS <- list(
  list(spec = spec_3block(sA = 0.30, sB = 0.60, sU = 1.0), n_data = 50, n_search = 60),
  list(spec = spec_3block_null(s = 0.30),                  n_data = 30, n_search = 60),
  list(spec = spec_5block(sW = 0.30, sB = 0.30, sU = 1.0), n_data = 30, n_search = 100)
)

if (QUICK) {
  SCENARIOS <- lapply(SCENARIOS[c(1, 3)], function(s) { s$n_data <- 2; s$n_search <- 6; s })
  N_TRAIN_GRID <- c(100, 200); N_SWEEP <- 4; INNER_FOLDS <- 3; INNER_REPEATS <- 2
}
OUT_DIR <- if (QUICK) "sim_n_quick" else "sim_n"
dir.create(OUT_DIR, showWarnings = FALSE)

L <- LOGIT_RANGE

## ---- Helpers -----------------------------------------------------------------
# Empirical design from the TRAINING data: mean |cor| between blocks.
emp_coord <- function(X) {
  nb <- length(X); M <- matrix(0, nb, nb)
  for (i in 1:(nb - 1)) for (j in (i + 1):nb) M[i, j] <- mean(abs(cor(X[[i]], X[[j]])))
  lg(M[upper.tri(M)])
}

# Inner-CV BER of every candidate: rows = candidates, cols = (repeat, fold).
# FID is an n x INNER_REPEATS matrix of fold labels, drawn once per dataset and
# shared by all candidates/arms. No set.seed() here: reseeding inside a worker
# would repeat datasets.
inner_scores <- function(coords, X, Y, bn, FID) {
  K <- max(FID); R <- ncol(FID)
  t(vapply(seq_len(nrow(coords)), function(i) {
    des <- coord_to_design(coords[i, ], bn)
    unlist(lapply(seq_len(R), function(r) vapply(seq_len(K), function(k) {
      te <- which(FID[, r] == k); tr <- setdiff(seq_along(Y), te)
      Ytr <- droplevels(Y[tr])
      if (nlevels(Ytr) < nlevels(Y)) return(NA_real_)
      ber_outer_eval(des, lapply(X, function(b) b[tr, , drop = FALSE]), Ytr,
                     lapply(X, function(b) b[te, , drop = FALSE]), Y[te])
    }, numeric(1))))
  }, numeric(K * R)))
}

# Ties go to the first candidate in order.
pick_min <- function(S) {
  m <- rowMeans(S, na.rm = TRUE)
  which.min(replace(m, is.na(m), Inf))
}

# One-standard-error rule: among candidates within 1 SE of the best mean,
# take the one with the smallest total coupling (sum of design entries).
pick_1se <- function(S, coords) {
  m  <- rowMeans(S, na.rm = TRUE); m[is.na(m)] <- Inf
  j  <- which.min(m)
  # SE from the K folds of ONE repeat (repeats are correlated, so pooling all
  # K*R scores would understate it).
  se <- sd(S[j, seq_len(INNER_FOLDS)], na.rm = TRUE) / sqrt(INNER_FOLDS)
  ok <- which(m <= m[j] + se)
  ok[which.min(rowSums(il(coords[ok, , drop = FALSE])))]
}

## ---- One dataset ---------------------------------------------------------------
run_dataset <- function(sc, id, n_train) {
  spec <- sc$spec; bn <- names(spec$blocks); d <- length(bn) * (length(bn) - 1) / 2
  pn   <- pair_index(bn)$names
  ds   <- simulate_dataset(spec, n_train, N_TEST)
  Xtr  <- ds$train$X; Ytr <- ds$train$Y
  FID  <- sapply(seq_len(INNER_REPEATS), function(r)
    caret::createFolds(Ytr, k = INNER_FOLDS, list = FALSE))

  ce <- emp_coord(Xtr)
  sw <- uniform_sweep(d, N_SWEEP, L)
  se <- generate_candidates(d, sc$n_search, ce, space = "logit")$coord

  S_sw <- inner_scores(sw, Xtr, Ytr, bn, FID)
  S_se <- inner_scores(se, Xtr, Ytr, bn, FID)

  j_sw <- pick_min(S_sw); j_se <- pick_min(S_se); j_1s <- pick_1se(S_se, se)
  sel <- list(
    sweep      = list(c = sw[j_sw, ], inner = mean(S_sw[j_sw, ], na.rm = TRUE), fits = length(S_sw)),
    search     = list(c = se[j_se, ], inner = mean(S_se[j_se, ], na.rm = TRUE), fits = length(S_se)),
    search_1se = list(c = se[j_1s, ], inner = mean(S_se[j_1s, ], na.rm = TRUE), fits = length(S_se)),
    empirical  = list(c = ce,                          inner = NA, fits = 0),
    null       = list(c = rep(-L, d),                  inner = NA, fits = 0),
    full       = list(c = rep(L, d),                   inner = NA, fits = 0),
    oracle     = list(c = pattern_coord(oracle_pairs(spec)), inner = NA, fits = 0))

  do.call(rbind, lapply(names(sel), function(a) {
    s <- sel[[a]]
    tb <- ber_outer_eval(coord_to_design(s$c, bn), Xtr, Ytr, ds$test$X, ds$test$Y)
    data.frame(scenario = spec$name, n_train = n_train, dataset = id, arm = a,
               test_ber = tb, inner_ber = s$inner, n_fits = s$fits,
               t(setNames(round(il(s$c), 4), pn)), check.names = FALSE)
  }))
}

## ---- Run -----------------------------------------------------------------------
ARMS <- c("sweep", "search", "search_1se", "empirical", "null", "full", "oracle")

summarise_cell <- function(R, hi) {
  w <- reshape(R[, c("dataset", "arm", "test_ber")], idvar = "dataset",
               timevar = "arm", direction = "wide")
  cat("  test BER (mean) | optimism = test - inner:\n")
  for (a in ARMS) {
    x <- R$test_ber[R$arm == a]; inn <- R$inner_ber[R$arm == a]
    cat(sprintf("    %-11s %.4f%s\n", a, mean(x, na.rm = TRUE),
                if (all(is.na(inn))) "" else
                  sprintf("   optimism %+.4f", mean(x - inn, na.rm = TRUE))))
  }
  cat("  paired vs SWEEP (negative = better than sweep):\n")
  for (a in setdiff(ARMS, "sweep")) {
    dif <- w[[paste0("test_ber.", a)]] - w$test_ber.sweep
    n   <- sum(!is.na(dif)); se <- sd(dif, na.rm = TRUE) / sqrt(n)
    cat(sprintf("    %-11s %+.4f  (SE %.4f, z = %5.1f)  better %2d / %d, tied %d\n",
                a, mean(dif, na.rm = TRUE), se, mean(dif, na.rm = TRUE) / se,
                sum(dif < 0, na.rm = TRUE), n, sum(dif == 0, na.rm = TRUE)))
  }
  if (any(hi) && any(!hi)) {
    cat("  structure recovery (planted pairs coupled more than the rest):\n")
    for (a in c("search", "search_1se")) {
      E <- as.matrix(R[R$arm == a, names(hi), drop = FALSE])
      rec <- rowMeans(E[, hi, drop = FALSE]) > rowMeans(E[, !hi, drop = FALSE])
      cat(sprintf("    %-11s %d / %d datasets\n", a, sum(rec), length(rec)))
    }
  }
}

for (sc in SCENARIOS) for (n_train in N_TRAIN_GRID) {
  hi   <- oracle_pairs(sc$spec)
  file <- file.path(OUT_DIR, sprintf("%s_n%d.csv", sc$spec$name, n_train))
  cat(sprintf("\n=== %s | n_train %d | %d datasets | planted coupled: %s ===\n",
              sc$spec$name, n_train, sc$n_data, paste(names(which(hi)), collapse = ", ")))
  if (file.exists(file)) {
    cat("  (already on disk, skipping -- delete the file to rerun)\n")
    summarise_cell(read.csv(file, check.names = FALSE), hi)
    next
  }
  t0 <- Sys.time()
  # Cell-specific seed: reproducible per cell, independent of which cells ran before.
  set.seed(match(sc$spec$name, vapply(SCENARIOS, function(z) z$spec$name, "")) * 10000 + n_train)
  out <- mclapply(seq_len(sc$n_data), function(i)
    tryCatch(run_dataset(sc, i, n_train), error = function(e) conditionMessage(e)),
    mc.cores = NCORES)
  bad <- !vapply(out, is.data.frame, logical(1))
  if (any(bad)) cat(sprintf("  !! %d dataset(s) failed: %s\n", sum(bad), unlist(out[bad])[1]))
  R <- do.call(rbind, out[!bad])
  write.csv(R, file, row.names = FALSE)
  cat(sprintf("  (%.1f min) -> %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), file))
  summarise_cell(R, hi)
}

## ---- Combine all cells + the headline table ------------------------------------
files <- list.files(OUT_DIR, pattern = "\\.csv$", full.names = TRUE)
cells <- lapply(files, read.csv, check.names = FALSE)
cols  <- unique(unlist(lapply(cells, names)))
res   <- do.call(rbind, lapply(cells, function(R) {
  for (m in setdiff(cols, names(R))) R[[m]] <- NA
  R[, cols]
}))
write.csv(res, sprintf("sim_samplesize%s.csv", if (QUICK) "_quick" else ""), row.names = FALSE)

cat("\n== HEADLINE: search - sweep (test BER) and optimism, by scenario and n ==\n")
for (scn in unique(res$scenario)) for (n in sort(unique(res$n_train[res$scenario == scn]))) {
  R <- res[res$scenario == scn & res$n_train == n, ]
  w <- reshape(R[, c("dataset", "arm", "test_ber")], idvar = "dataset",
               timevar = "arm", direction = "wide")
  dif <- w$test_ber.search - w$test_ber.sweep
  se  <- sd(dif, na.rm = TRUE) / sqrt(sum(!is.na(dif)))
  orc <- mean(w$test_ber.oracle - w$test_ber.sweep, na.rm = TRUE)
  opt <- function(a) with(R[R$arm == a, ], mean(test_ber - inner_ber, na.rm = TRUE))
  cat(sprintf("  %-28s n=%3d  search-sweep %+.4f (z %5.1f) | oracle-sweep %+.4f | optimism sweep %+.3f search %+.3f\n",
              scn, n, mean(dif, na.rm = TRUE), mean(dif, na.rm = TRUE) / se, orc,
              opt("sweep"), opt("search")))
}
cat(sprintf("\nWritten: sim_samplesize%s.csv\n", if (QUICK) "_quick" else ""))
