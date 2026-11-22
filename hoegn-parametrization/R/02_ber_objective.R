# =============================================================================
# 02_ber_objective.R
# -----------------------------------------------------------------------------
# BER machinery extracted from 5_public_ber_selection_rules.R, rewritten as a
# pluggable objective for 01_parametrization_core.R. This is the PRIMARY
# demonstration objective (plan step 2), pending Aaron's confirmation.
#
# Two functions matter:
#   ber_objective(design, Xd, Yd)              inner-CV BER, used for SEARCH
#   ber_outer_eval(design, Xtr, Ytr, Xte, Yte)  held-out BER, used for EVALUATION
# =============================================================================

source("R/01_parametrization_core.R")

# Primary distance is "centroids.dist" (the default below). Report all three
# (max.dist, centroids.dist, mahalanobis.dist) in a robustness table.

ber <- function(truth, pred) {
  truth <- factor(truth)
  pred  <- factor(as.character(pred), levels = levels(truth))
  mean(tapply(seq_along(truth), truth,
              function(i) mean(pred[i] != truth[i])), na.rm = TRUE)
}

get_class <- function(pr, dist, comp) {
  for (p in list(c("WeightedVote", dist), c("MajorityVote", dist))) {
    x <- pr[[p[1]]]; if (is.null(x)) next
    x <- x[[p[2]]];  if (is.null(x)) next
    if (is.matrix(x)) return(as.character(x[, min(comp, ncol(x))]))
  }
  NULL
}

# Held-out BER of one design: fit on train, predict test.
ber_outer_eval <- function(design, Xtr, Ytr, Xte, Yte, dist = "centroids.dist") {
  fit <- fit_block(Xtr, Ytr, design)
  if (is.null(fit)) return(NA_real_)
  pr <- tryCatch(predict(fit, newdata = Xte, dist = dist), error = function(e) NULL)
  if (is.null(pr)) return(NA_real_)
  cl <- get_class(pr, dist, NCOMP)
  if (is.null(cl)) NA_real_ else ber(Yte, cl)
}

# Inner-CV BER: the search objective. Uses only the data it is given, so inside
# nested_select() it never sees the outer test fold.
ber_objective <- function(design, Xd, Yd, inner_folds = 5, n_repeat = 1,
                          seed = 1, dist = "centroids.dist") {
  vals <- c()
  for (r in seq_len(n_repeat)) {
    set.seed(seed + r)
    fid <- caret::createFolds(Yd, k = inner_folds, list = FALSE)
    for (k in seq_len(inner_folds)) {
      te <- which(fid == k); tr <- setdiff(seq_along(Yd), te)
      Xtr <- lapply(Xd, function(b) b[tr, , drop = FALSE])
      Xte <- lapply(Xd, function(b) b[te, , drop = FALSE])
      Ytr <- droplevels(Yd[tr])
      if (nlevels(Ytr) < nlevels(Yd)) next
      vals <- c(vals, ber_outer_eval(design, Xtr, Ytr, Xte, Yd[te], dist))
    }
  }
  mean(vals, na.rm = TRUE)
}

# Paired comparison helper: shared folds across designs so differences are paired.
# designs: named list of design matrices. Returns long data.frame of BER.
paired_ber <- function(designs, X, Y, folds = 10, n_repeat = 10, seed = 42,
                       distances = c("max.dist", "centroids.dist", "mahalanobis.dist")) {
  rows <- list()
  for (r in seq_len(n_repeat)) {
    set.seed(seed + r)
    fid <- caret::createFolds(Y, k = folds, list = FALSE)
    for (k in seq_len(folds)) {
      te <- which(fid == k); tr <- setdiff(seq_along(Y), te)
      Xtr <- lapply(X, function(b) b[tr, , drop = FALSE])
      Xte <- lapply(X, function(b) b[te, , drop = FALSE])
      Ytr <- droplevels(Y[tr])
      if (nlevels(Ytr) < nlevels(Y)) next
      for (dn in names(designs)) for (dist in distances) {
        rows[[length(rows) + 1]] <- data.frame(
          rep = r, fold = k, design = dn, distance = dist,
          BER = ber_outer_eval(designs[[dn]], Xtr, Ytr, Xte, Y[te], dist))
      }
    }
  }
  do.call(rbind, rows)
}
