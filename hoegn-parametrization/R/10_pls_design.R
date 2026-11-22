# =============================================================================
# 10_pls_design.R
# -----------------------------------------------------------------------------
# Adds one fixed reference design to the real-data comparison (03): the
# PLS-based design recommended by mixOmics. For each block pair, a one-component
# PLS is fitted between the two blocks, and the design entry is the correlation
# between the two component scores. Computed on the OUTER TRAINING fold only.
#
# A fixed design needs no selection, so 03 does not have to be rerun: this script
# rebuilds the SAME outer folds (same seeds as 03), fits the PLS design, and
# scores it on each outer test fold.
#
# Check: it also re-scores the null design. Its BER must match the null rows in
# the existing arms_<tag>.csv; if it does, the folds are identical and the new
# arm is paired with all others.
#
# Reads : public_datasets.rds, arms_<tag>.csv (for the check and the summary)
# Writes: arms_<tag>_pls.csv  (rows: arm = "pls"; entries, outer_ber)
# Run from the project root:  Rscript R/10_pls_design.R
# =============================================================================

source("R/02_ber_objective.R")

RUNS <- list(                       # tag, dataset, ncomp, outer repeats (as in 03)
  list(tag = "nutrimouse",          ds = "nutrimouse",   ncomp = 2, reps = 5),
  list(tag = "tcga",                ds = "tcga",         ncomp = 2, reps = 5),
  list(tag = "mogonet_brca",        ds = "mogonet_brca", ncomp = 2, reps = 5),
  list(tag = "mogonet_brca_ncomp4", ds = "mogonet_brca", ncomp = 4, reps = 3))
OUTER_FOLDS <- 5
SEED        <- 42                   # must equal SEED in 03
L           <- LOGIT_RANGE

ALL <- readRDS("public_datasets.rds")

# PLS design: entry (i, j) = cor of first PLS components of blocks i and j.
pls_design <- function(Xd) {
  bn <- names(Xd); nb <- length(Xd)
  M <- matrix(0, nb, nb, dimnames = list(bn, bn))
  for (i in 1:(nb - 1)) for (j in (i + 1):nb) {
    fit <- suppressMessages(suppressWarnings(
      mixOmics::pls(Xd[[i]], Xd[[j]], ncomp = 1, scale = TRUE)))
    M[i, j] <- M[j, i] <- abs(cor(fit$variates$X[, 1], fit$variates$Y[, 1]))
  }
  M
}

for (run in RUNS) {
  NCOMP <- run$ncomp                # fit_block() and get_class() read NCOMP
  D  <- ALL[[run$ds]]; X <- D$X; Y <- D$Y; bn <- names(X); nb <- length(X)
  d  <- nb * (nb - 1) / 2
  idx <- which(upper.tri(diag(nb)), arr.ind = TRUE)
  PAIR_NAMES <- paste(bn[idx[, 1]], bn[idx[, 2]], sep = "-")
  cat(sprintf("\n== %s (ncomp %d) ==\n", run$tag, NCOMP))

  rows <- list(); check <- list()
  for (r in seq_len(run$reps)) {
    set.seed(SEED + r)              # identical to 03
    folds <- caret::createFolds(Y, k = OUTER_FOLDS, list = TRUE, returnTrain = FALSE)
    for (f in seq_along(folds)) {
      te  <- folds[[f]]; tr <- setdiff(seq_along(Y), te)
      Xtr <- lapply(X, function(b) b[tr, , drop = FALSE])
      Xte <- lapply(X, function(b) b[te, , drop = FALSE])
      Ytr <- droplevels(Y[tr]); Yte <- Y[te]
      if (nlevels(Ytr) < nlevels(Y)) next

      M  <- pls_design(Xtr)
      ob <- ber_outer_eval(M, Xtr, Ytr, Xte, Yte)
      rows[[length(rows) + 1]] <- data.frame(
        rep = r, fold = f, arm = "pls",
        t(setNames(round(M[upper.tri(M)], 4), PAIR_NAMES)),
        inner_ber = NA, outer_ber = ob, n_fits = 0, check.names = FALSE)

      nul <- ber_outer_eval(coord_to_design(rep(-L, d), bn), Xtr, Ytr, Xte, Yte)
      check[[length(check) + 1]] <- data.frame(rep = r, fold = f, null_new = nul)
      cat(sprintf("  rep %d fold %d  pls BER %.4f | entries %s\n", r, f, ob,
                  paste(sprintf("%.3f", M[upper.tri(M)]), collapse = " ")))
    }
  }
  res <- do.call(rbind, rows)
  write.csv(res, sprintf("arms_%s_pls.csv", run$tag), row.names = FALSE)

  ## ---- Fold check and summary against the existing 03 results ---------------
  old_file <- sprintf("arms_%s.csv", run$tag)
  if (!file.exists(old_file)) { cat("  (no", old_file, "- skipping check)\n"); next }
  old <- read.csv(old_file, check.names = FALSE, na.strings = "NA")
  chk <- merge(do.call(rbind, check),
               old[old$arm == "null", c("rep", "fold", "outer_ber")], by = c("rep", "fold"))
  ok  <- isTRUE(all.equal(chk$null_new, chk$outer_ber, tolerance = 1e-8))
  cat(sprintf("  FOLD CHECK (null design re-scored): %s (max diff %.2e)\n",
              if (ok) "OK, folds identical" else "MISMATCH, do not merge",
              max(abs(chk$null_new - chk$outer_ber))))

  sw <- old[old$arm == "sweep", c("rep", "fold", "outer_ber")]
  em <- old[old$arm == "empirical", c("rep", "fold", "outer_ber")]
  m  <- merge(merge(res[, c("rep", "fold", "outer_ber")], sw, by = c("rep", "fold"),
                    suffixes = c("", ".sweep")),
              em, by = c("rep", "fold"), suffixes = c("", ".emp"))
  for (ref in c("sweep", "emp")) {
    dif <- m$outer_ber - m[[paste0("outer_ber.", ref)]]
    cat(sprintf("  pls vs %-5s: mean %+.4f | better/tied/worse %d/%d/%d\n", ref,
                mean(dif), sum(dif < -1e-9), sum(abs(dif) <= 1e-9), sum(dif > 1e-9)))
  }
  cat(sprintf("  pls mean outer BER %.4f (SD %.4f); median entries %s\n",
              mean(res$outer_ber), sd(res$outer_ber),
              paste(sprintf("%s %.3f", PAIR_NAMES,
                            apply(res[, PAIR_NAMES, drop = FALSE], 2, median)),
                    collapse = " | ")))
}
