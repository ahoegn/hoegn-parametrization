# =============================================================================
# 00_data_loader.R   (KEPT AS-IS from 0_public_data_loader.R)
# -----------------------------------------------------------------------------
# Loads the public multiblock datasets and reports:
#   - block dimensions and row alignment
#   - class balance and the maximum viable fold count
#   - the empirical design matrix (mean |inter-block correlation|)
#
# Writes: public_datasets.rds  (X / Y / design / folds per dataset)
#
# TODO (plan step 5): add one 4-5 block dataset (d = 6-10) here.
# =============================================================================
setwd("/run/media/aaron/Verbatim/Projects/hoegn_parametrization_kept")

library(mixOmics)
cat("mixOmics version:", as.character(packageVersion("mixOmics")), "\n\n")

## ---- Shared helpers ---------------------------------------------------------
# NOTE: this is the "empirical" design used throughout: mean absolute
# inter-block correlation. It is NOT the PLS-component correlation design from
# the mixOmics literature. Describe it accurately in Methods.
compute_design_matrix <- function(blocks) {
  nb <- length(blocks)
  M  <- matrix(NA_real_, nb, nb, dimnames = list(names(blocks), names(blocks)))
  for (i in seq_len(nb)) for (j in i:nb) {
    m <- mean(abs(cor(blocks[[i]], blocks[[j]]))); M[i, j] <- m; M[j, i] <- m
  }
  round(M, 3)
}

# Largest k such that every class still has >= min_per_class members per fold
max_viable_folds <- function(Y, min_per_class = 2) {
  smallest <- min(table(Y))
  max(2L, floor(smallest / min_per_class))
}

describe <- function(label, X, Y) {
  cat("\n===============================================================\n")
  cat(label, "\n")
  cat("===============================================================\n")
  cat("n subjects :", nrow(X[[1]]), "\n")
  cat("blocks     :", paste(sprintf("%s(%d)", names(X), sapply(X, ncol)), collapse = ", "), "\n")
  stopifnot(all(sapply(X, nrow) == length(Y)))
  cat("class table:\n"); print(table(Y))
  kmax <- max_viable_folds(Y)
  cat(sprintf("max viable folds (>=2 per class): %d\n", kmax))
  if (kmax < 5) cat("  !! too few per class for stable k-fold - reconsider this target\n")
  M <- compute_design_matrix(X)
  cat("empirical design matrix:\n"); print(M)
  invisible(list(X = X, Y = Y, design = M, folds = min(kmax, 10L)))
}

## ---- breast.TCGA ------------------------------------------------------------
# Training arm only; the shipped train/test split is DIABLO's, not ours.
data(breast.TCGA)
tcga_X <- list(
  mRNA    = as.matrix(breast.TCGA$data.train$mrna),
  miRNA   = as.matrix(breast.TCGA$data.train$mirna),
  Protein = as.matrix(breast.TCGA$data.train$protein)
)
tcga_Y <- factor(breast.TCGA$data.train$subtype)
tcga <- describe("breast.TCGA  (3 blocks, target = PAM50 subtype)", tcga_X, tcga_Y)

## ---- nutrimouse -------------------------------------------------------------
# TARGET = genotype (2 levels, 20/20). Diet (5 levels, 8/class) breaks k-fold.
data(nutrimouse)
nutri_X <- list(
  Gene  = as.matrix(nutrimouse$gene),
  Lipid = as.matrix(nutrimouse$lipid)
)
nutri_Y <- factor(nutrimouse$genotype)
nutri <- describe("nutrimouse  (2 blocks, target = genotype)", nutri_X, nutri_Y)

## ---- Free design-matrix parameters ------------------------------------------
# nutrimouse : 2 blocks -> d = 1 (d-D search reduces to the 1-D sweep: sanity check)
# breast.TCGA: 3 blocks -> d = 3 (exhaustive grid feasible: ground truth for optimisers)
cat("\nFree off-diagonal entries -- nutrimouse: 1 | breast.TCGA: 3\n")

saveRDS(list(tcga = tcga, nutrimouse = nutri), "public_datasets.rds")
cat("\nSaved: public_datasets.rds\n")
