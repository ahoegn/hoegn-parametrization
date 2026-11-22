# =============================================================================
# 05_sim_pilot.R
# -----------------------------------------------------------------------------
# PILOT: is the planted scenario informative at all?
#
# No search here. For each simulated dataset we fit a set of FIXED designs on
# the training set and score them on a large independent test set:
#   uniform_c   the 1-D sweep, at 9 values of c (null ... full)
#   oracle      the planted design (coupled pairs high, the rest low)
#   anti        the complement of the planted design
#
# The scenario is informative iff the oracle beats the BEST uniform design.
# "Best uniform" is chosen on the TEST set, i.e. the sweep is given an unfair
# advantage. If the oracle still wins, there is real structure off the diagonal
# that no sweep can reach -- and the full nested simulation is worth running.
# On the control scenario the oracle IS a uniform design, so the gap must be ~0.
#
# Run from the project folder (the one containing R/). Takes a few minutes.
# =============================================================================

source("R/02_ber_objective.R")
source("R/04_sim_generator.R")
suppressPackageStartupMessages(library(parallel))

## ---- CONFIG -----------------------------------------------------------------
N_DATASETS <- 50
N_TRAIN    <- 100
N_TEST     <- 1000
NCORES     <- max(1L, detectCores() - 1L)
RNGkind("L'Ecuyer-CMRG"); set.seed(42)      # reproducible across forked workers

# Pilot v2. v1 showed: coupling HELPS when partners share class signal (control:
# 0.228 -> 0.203 from null to full) and HURTS when a partner carries nuisance.
# The off-diagonal case pays off most when both effects are strong, so B is no
# longer much stronger than A (a strong B predicts well alone and dilutes the gain).
SCENARIOS <- list(
  spec_3block(sA = 0.30, sB = 0.60, sU = 1.0),   # v1 scenario 2, rerun with n = 50
  spec_3block(sA = 0.30, sB = 0.30, sU = 1.0),   # balanced class signal, strong nuisance
  spec_3block(sA = 0.30, sB = 0.30, sU = 1.5),   # balanced, stronger nuisance
  spec_3block_null(s = 0.30),                    # control: planted design ON the diagonal
  spec_5block(sW = 0.30, sB = 0.30, sU = 1.0)    # d = 10
)

C_SWEEP <- c(0.0025, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.9975)

## ---- One dataset -> test BER of every fixed design ---------------------------
run_one <- function(spec) {
  ds <- simulate_dataset(spec, N_TRAIN, N_TEST)
  bn <- names(spec$blocks); d <- length(bn) * (length(bn) - 1) / 2
  hi <- oracle_pairs(spec)

  coords <- c(
    setNames(lapply(C_SWEEP, function(c) rep(lg(c), d)), sprintf("uniform_%g", C_SWEEP)),
    list(oracle = pattern_coord(hi), anti = pattern_coord(!hi)))

  vapply(coords, function(v)
    ber_outer_eval(coord_to_design(v, bn), ds$train$X, ds$train$Y,
                   ds$test$X, ds$test$Y),
    numeric(1))
}

## ---- Run ---------------------------------------------------------------------
results <- list()
for (spec in SCENARIOS) {
  t0 <- Sys.time()
  cat(sprintf("\n%s  | planted coupled pairs: %s\n", spec$name,
              paste(names(which(oracle_pairs(spec))), collapse = ", ")))
  out <- mclapply(seq_len(N_DATASETS), function(i) run_one(spec), mc.cores = NCORES)
  ok  <- vapply(out, is.numeric, logical(1))
  if (any(!ok)) cat(sprintf("  !! %d dataset(s) failed\n", sum(!ok)))
  B <- do.call(rbind, out[ok])
  results[[spec$name]] <- B

  uni   <- grep("^uniform_", colnames(B), value = TRUE)
  m     <- colMeans(B, na.rm = TRUE)
  # Fair baseline: the single sweep level that is best ON AVERAGE across datasets.
  # (Taking the per-dataset minimum over 9 noisy test BERs is biased low.)
  fixed <- uni[which.min(m[uni])]
  gap   <- B[, "oracle"] - B[, fixed]                    # negative = oracle wins
  se    <- sd(gap, na.rm = TRUE) / sqrt(sum(!is.na(gap)))
  best_uni <- apply(B[, uni, drop = FALSE], 1, min, na.rm = TRUE)

  cat(sprintf("  (%.1f min, %d datasets)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(B)))
  cat("  mean test BER per design:\n")
  print(round(m, 3))
  cat(sprintf("  best fixed sweep level                 : %s (%.3f)\n", fixed, m[fixed]))
  cat(sprintf("  oracle - best fixed sweep level        : %+.4f  (SE %.4f, z = %.1f)\n",
              mean(gap, na.rm = TRUE), se, mean(gap, na.rm = TRUE) / se))
  cat(sprintf("  oracle better in                       : %d / %d datasets (tied %d)\n",
              sum(gap < 0, na.rm = TRUE), sum(!is.na(gap)), sum(gap == 0, na.rm = TRUE)))
  cat(sprintf("  oracle - anti                          : %+.4f\n",
              mean(B[, "oracle"] - B[, "anti"], na.rm = TRUE)))
  cat(sprintf("  [optimistic sweep bound, test-tuned per dataset: %.3f]\n", mean(best_uni)))
}

saveRDS(results, "sim_pilot_v2_results.rds")
cat("\nSaved: sim_pilot_v2_results.rds\n")
cat("Informative: 'oracle - best fixed sweep level' clearly negative (z < -2) on the\n")
cat("off-diagonal scenarios, and ~0 on the control (where the oracle IS a sweep level).\n")
