# =============================================================================
# 06_sim_nested.R
# -----------------------------------------------------------------------------
# Full simulation: does d-dimensional SEARCH beat the 1-D SWEEP when each
# method has to choose its own design from the training data?
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
# Writes: sim_nested.csv  (one row per scenario x dataset x arm)
# Set QUICK <- TRUE first for a smoke test (~2-3 min).
# =============================================================================

source("R/02_ber_objective.R")
source("R/04_sim_generator.R")
suppressPackageStartupMessages(library(parallel))

## ---- CONFIG -----------------------------------------------------------------
QUICK       <- FALSE
N_TRAIN     <- 100
N_TEST      <- 1000
INNER_FOLDS <- 5
N_SWEEP     <- 15
NCORES      <- max(1L, detectCores() - 1L)
RNGkind("L'Ecuyer-CMRG"); set.seed(7)

SCENARIOS <- list(
  list(spec = spec_3block(sA = 0.30, sB = 0.60, sU = 1.0), n_data = 30, n_search = 60),
  list(spec = spec_3block(sA = 0.30, sB = 0.30, sU = 1.0), n_data = 30, n_search = 60),
  list(spec = spec_3block(sA = 0.30, sB = 0.30, sU = 1.5), n_data = 30, n_search = 60),
  list(spec = spec_3block_null(s = 0.30),                  n_data = 30, n_search = 60),
  list(spec = spec_5block(sW = 0.30, sB = 0.30, sU = 1.0), n_data = 20, n_search = 150)
)

if (QUICK) {
  SCENARIOS <- lapply(SCENARIOS[c(2, 4)], function(s) { s$n_data <- 2; s$n_search <- 6; s })
  N_SWEEP <- 4; INNER_FOLDS <- 3
}

L <- LOGIT_RANGE

## ---- Helpers -----------------------------------------------------------------
# Empirical design from the TRAINING data: mean |cor| between blocks.
emp_coord <- function(X) {
  nb <- length(X); M <- matrix(0, nb, nb)
  for (i in 1:(nb - 1)) for (j in (i + 1):nb) M[i, j] <- mean(abs(cor(X[[i]], X[[j]])))
  lg(M[upper.tri(M)])
}

# Inner-CV BER of every candidate: rows = candidates, cols = inner folds.
# Folds (fid) are drawn once per dataset and shared by all candidates/arms.
# No set.seed() here: reseeding inside a worker would repeat datasets.
inner_scores <- function(coords, X, Y, bn, fid) {
  t(vapply(seq_len(nrow(coords)), function(i) {
    des <- coord_to_design(coords[i, ], bn)
    vapply(seq_len(max(fid)), function(k) {
      te <- which(fid == k); tr <- setdiff(seq_along(Y), te)
      Ytr <- droplevels(Y[tr])
      if (nlevels(Ytr) < nlevels(Y)) return(NA_real_)
      ber_outer_eval(des, lapply(X, function(b) b[tr, , drop = FALSE]), Ytr,
                     lapply(X, function(b) b[te, , drop = FALSE]), Y[te])
    }, numeric(1))
  }, numeric(max(fid))))
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
  se <- sd(S[j, ], na.rm = TRUE) / sqrt(sum(!is.na(S[j, ])))
  ok <- which(m <= m[j] + se)
  ok[which.min(rowSums(il(coords[ok, , drop = FALSE])))]
}

## ---- One dataset ---------------------------------------------------------------
run_dataset <- function(sc, id) {
  spec <- sc$spec; bn <- names(spec$blocks); d <- length(bn) * (length(bn) - 1) / 2
  pn   <- pair_index(bn)$names
  ds   <- simulate_dataset(spec, N_TRAIN, N_TEST)
  Xtr  <- ds$train$X; Ytr <- ds$train$Y
  fid  <- caret::createFolds(Ytr, k = INNER_FOLDS, list = FALSE)

  ce <- emp_coord(Xtr)
  sw <- uniform_sweep(d, N_SWEEP, L)
  se <- generate_candidates(d, sc$n_search, ce, space = "logit")$coord

  S_sw <- inner_scores(sw, Xtr, Ytr, bn, fid)
  S_se <- inner_scores(se, Xtr, Ytr, bn, fid)

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
    data.frame(scenario = spec$name, dataset = id, arm = a,
               test_ber = tb, inner_ber = s$inner, n_fits = s$fits,
               t(setNames(round(il(s$c), 4), pn)), check.names = FALSE)
  }))
}

## ---- Run -----------------------------------------------------------------------
all_rows <- list()
for (sc in SCENARIOS) {
  t0 <- Sys.time()
  hi <- oracle_pairs(sc$spec)
  cat(sprintf("\n=== %s | %d datasets | planted coupled: %s ===\n", sc$spec$name,
              sc$n_data, paste(names(which(hi)), collapse = ", ")))

  out <- mclapply(seq_len(sc$n_data), function(i)
    tryCatch(run_dataset(sc, i), error = function(e) conditionMessage(e)),
    mc.cores = NCORES)
  bad <- !vapply(out, is.data.frame, logical(1))
  if (any(bad)) cat(sprintf("  !! %d dataset(s) failed: %s\n", sum(bad), unlist(out[bad])[1]))
  R <- do.call(rbind, out[!bad]); all_rows[[sc$spec$name]] <- R
  cat(sprintf("  (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  ARMS <- c("sweep", "search", "search_1se", "empirical", "null", "full", "oracle")
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
    cat("  structure recovery (selected design couples planted pairs more than the rest):\n")
    for (a in c("search", "search_1se")) {
      E <- as.matrix(R[R$arm == a, names(hi), drop = FALSE])
      rec <- rowMeans(E[, hi, drop = FALSE]) > rowMeans(E[, !hi, drop = FALSE])
      cat(sprintf("    %-11s %d / %d datasets\n", a, sum(rec), length(rec)))
    }
  }
}

res <- do.call(rbind, lapply(all_rows, function(R) {
  missing <- setdiff(unique(unlist(lapply(all_rows, names))), names(R))
  for (m in missing) R[[m]] <- NA
  R[, unique(unlist(lapply(all_rows, names)))]
}))
write.csv(res, "sim_nested.csv", row.names = FALSE)
cat("\nWritten: sim_nested.csv\n")
