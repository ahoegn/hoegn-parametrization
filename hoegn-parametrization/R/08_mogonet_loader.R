# =============================================================================
# 08_mogonet_loader.R
# -----------------------------------------------------------------------------
# Loads the MOGONET BRCA benchmark (Wang et al. 2021, Nat Commun;
# https://github.com/txWang/MOGONET) and adds it to public_datasets.rds as
# "mogonet_brca", in the same format as tcga / nutrimouse, so that
# 03_compare_arms.R runs on it unchanged (DATASET <- "mogonet_brca").
#
# Expected layout (git clone https://github.com/txWang/MOGONET.git data/MOGONET):
#   data/MOGONET/BRCA/{1,2,3}_{tr,te}.csv   numeric, no header, comma-separated
#   data/MOGONET/BRCA/{1,2,3}_featname.csv  feature names
#   data/MOGONET/BRCA/labels_{tr,te}.csv    class codes 0..4 (PAM50 subtypes)
# Views (MOGONET BRCA): 1 = mRNA expression, 2 = DNA methylation, 3 = miRNA.
#
# MOGONET's own train/test split is NOT used: train and test are merged and our
# nested CV generates its own splits, exactly as for the other datasets.
#
# Also prints the time of ONE model fit on the full data, to size the run.
# =============================================================================

suppressPackageStartupMessages(library(mixOmics))

DIR <- "data/MOGONET/BRCA"
stopifnot(dir.exists(DIR))

read_num <- function(f) as.matrix(read.csv(file.path(DIR, f), header = FALSE))
read_names <- function(f) {
  x <- read.csv(file.path(DIR, f), header = FALSE, stringsAsFactors = FALSE)[[1]]
  make.unique(as.character(x))
}

VIEWS <- c(`1` = "mRNA", `2` = "Methylation", `3` = "miRNA")

X <- lapply(names(VIEWS), function(v) {
  M  <- rbind(read_num(sprintf("%s_tr.csv", v)), read_num(sprintf("%s_te.csv", v)))
  fn <- read_names(sprintf("%s_featname.csv", v))
  if (length(fn) == ncol(M)) colnames(M) <- fn
  else { warning(sprintf("view %s: %d names for %d columns", v, length(fn), ncol(M)))
         colnames(M) <- paste0(VIEWS[[v]], "_", seq_len(ncol(M))) }
  M
})
names(X) <- VIEWS

lab <- c(read_num("labels_tr.csv")[, 1], read_num("labels_te.csv")[, 1])
Y   <- factor(paste0("subtype_", as.integer(round(lab))))
n   <- length(Y)
rn  <- sprintf("s%03d", seq_len(n))
X   <- lapply(X, function(m) { rownames(m) <- rn; m })

## ---- Checks ----------------------------------------------------------------
stopifnot(all(vapply(X, nrow, 1L) == n))
cat(sprintf("MOGONET BRCA: n = %d | %s\n", n,
            paste(sprintf("%s %d", names(X), vapply(X, ncol, 1L)), collapse = " | ")))
cat("classes:\n"); print(table(Y))
cat(sprintf("missing values: %d | zero-variance features: %d\n",
            sum(vapply(X, function(m) sum(is.na(m)), 1L)),
            sum(vapply(X, function(m) sum(apply(m, 2, sd) == 0), 1L))))

# Drop zero-variance features (block.splsda scales columns).
X <- lapply(X, function(m) m[, apply(m, 2, sd) > 0, drop = FALSE])

## ---- Empirical design (same definition as 00_data_loader.R) ----------------
nb <- length(X)
M  <- matrix(0, nb, nb, dimnames = list(names(X), names(X)))
for (i in 1:(nb - 1)) for (j in (i + 1):nb) {
  M[i, j] <- M[j, i] <- mean(abs(cor(X[[i]], X[[j]])))
}
diag(M) <- 1
M <- round(M, 3)
cat("empirical design (mean |cor|):\n"); print(M)

max_folds <- max(2L, floor(min(table(Y)) / 2))
cat(sprintf("smallest class: %d -> max viable folds (>= 2 per class): %d\n",
            min(table(Y)), max_folds))

## ---- Save alongside tcga / nutrimouse ---------------------------------------
P <- if (file.exists("public_datasets.rds")) readRDS("public_datasets.rds") else list()
P$mogonet_brca <- list(X = X, Y = Y, design = M, folds = min(max_folds, 10L),
                       source = "https://github.com/txWang/MOGONET (BRCA), train+test merged")
saveRDS(P, "public_datasets.rds")
cat("Saved: public_datasets.rds  (element 'mogonet_brca')\n")

## ---- Timing: one fit on the full data ---------------------------------------
keepX <- lapply(X, function(m) rep(max(1L, round(ncol(m) * 0.5)), 2))
des   <- matrix(0.1, nb, nb, dimnames = list(names(X), names(X))); diag(des) <- 0
t1 <- system.time(suppressWarnings(suppressMessages(
  fit <- block.splsda(X, Y, ncomp = 2, keepX = keepX, design = des, scale = TRUE))))
t2 <- system.time(suppressWarnings(
  pr <- predict(fit, newdata = lapply(X, function(m) m[1:175, ]), dist = "centroids.dist")))
cat(sprintf("\nONE fit (n = %d): %.1f s | predict 175 samples: %.1f s\n",
            n, t1["elapsed"], t2["elapsed"]))
cat("-> send this timing; it decides the settings for the 03 run.\n")
