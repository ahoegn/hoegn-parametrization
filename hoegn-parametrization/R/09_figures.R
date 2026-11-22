# =============================================================================
# 09_figures.R
# -----------------------------------------------------------------------------
# Paper figures. Reads only the result files written by 03 and 07; fits nothing.
#
#   Fig1  The parametrization: design matrix, logit map, and the design cube
#         (sweep = the diagonal, search = points filling the cube). No data.
#   Fig2  Simulation: search, 1-SE search and planted design minus the sweep,
#         by scenario and training size (mean, 95% CI over datasets).
#   Fig3  Selection optimism: (A) simulations by procedure and training size,
#         (B) real data against the number of candidates.
#   Fig4  Real data: (A) per-fold outer BER minus the sweep,
#         (B) the design entries each procedure selected, per fold.
#
# Inputs (project root): sim_samplesize.csv, arms_tcga.csv,
#   arms_mogonet_brca.csv, arms_mogonet_brca_ncomp4.csv
# Writes: figures/Fig1..Fig4 as .tiff (PLOS: LZW, 600 dpi, 7.5 in wide) and .pdf
# Needs: ggplot2, patchwork, caret.   install.packages(c("ggplot2", "patchwork"))
# Run from the project root:  Rscript R/09_figures.R
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

RES <- "."              # folder holding the result CSVs
OUT <- "figures"
DPI <- 600              # PLOS ONE accepts 300-600; 600 also meets Elsevier's 500
dir.create(OUT, showWarnings = FALSE)

rd <- function(f) read.csv(file.path(RES, f), na.strings = "NA",
                           stringsAsFactors = FALSE, check.names = FALSE)

## ---- Shared style --------------------------------------------------------------
# One colour per procedure, used in every figure. Fixed designs are grey.
COL <- c(sweep = "#6b6b6b", search = "#2a78d6", grid = "#eb6834",
         search_unit = "#1baf7a", search_1se = "#eda100", oracle = "#e87ba4",
         empirical = "#9a9a9a", null = "#9a9a9a", full = "#9a9a9a")
LAB <- c(sweep = "Sweep", search = "Search", grid = "Grid",
         search_unit = "Search on [0,1]", search_1se = "Search (1-SE)",
         oracle = "Planted", empirical = "Empirical", null = "Null", full = "Full")
SCEN <- c("3blk_control_s0.30"        = "Control",
          "3blk_sA0.30_sB0.60_sU1.00" = "Mixed 1",
          "5blk_sW0.30_sB0.30_sU1.00" = "Five-block")
REAL <- c("arms_tcga.csv"                = "breast.TCGA",
          "arms_mogonet_brca.csv"        = "MOGONET, ncomp 2",
          "arms_mogonet_brca_ncomp4.csv" = "MOGONET, ncomp 4")

theme_paper <- theme_classic(base_size = 8) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(size = 8, face = "bold"),
        legend.title     = element_blank(),
        legend.key.size  = unit(3, "mm"),
        plot.tag         = element_text(face = "bold"))
theme_set(theme_paper)

save_fig <- function(p, name, h) {
  ggsave(file.path(OUT, paste0(name, ".tiff")), p, width = 7.5, height = h,
         units = "in", dpi = DPI, compression = "lzw")
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = 7.5, height = h,
         units = "in", device = cairo_pdf)   # cairo: handles the en dash
  cat("wrote", name, "\n")
}

## ---- Helpers -------------------------------------------------------------------
# Paired differences vs the sweep within one dataset/fold set.
# d: long data with columns id, arm, <y>.
diff_vs_sweep <- function(d, arms, y) {
  sw <- d[d$arm == "sweep", c("id", y)]
  do.call(rbind, lapply(intersect(arms, unique(d$arm)), function(a) {
    x <- merge(d[d$arm == a, c("id", y)], sw, by = "id", suffixes = c("", ".sw"))
    data.frame(arm = a, id = x$id, diff = x[[y]] - x[[paste0(y, ".sw")]])
  }))
}
mean_se <- function(dd) do.call(rbind, lapply(split(dd$diff, dd$arm), function(v)
  data.frame(m = mean(v), se = sd(v) / sqrt(length(v)), n = length(v))))
arm_factor <- function(a, levels) factor(a, levels = levels, labels = LAB[levels])

sim <- rd("sim_samplesize.csv")
sim$scen <- factor(SCEN[sim$scenario], levels = c("Mixed 1", "Five-block", "Control"))
sim$id   <- sim$dataset

real_list <- lapply(names(REAL), function(f) {
  d <- rd(f)
  d$id <- paste(d$rep, d$fold)
  d$dataset <- REAL[[f]]
  d
})
# Datasets have different block pairs: pad missing pair columns with NA.
all_cols <- unique(unlist(lapply(real_list, names)))
real <- do.call(rbind, lapply(real_list, function(d) {
  d[setdiff(all_cols, names(d))] <- NA
  d[all_cols]
}))
real$dataset <- factor(real$dataset, levels = unname(REAL))
pair_cols <- setdiff(names(real), c("rep", "fold", "arm", "inner_ber", "outer_ber",
                                    "n_fits", "id", "dataset"))

## ---- Fig 1: the parametrization ---------------------------------------------------
L <- 6
# (A) the design matrix: hollow, symmetric, three free entries
mat <- expand.grid(col = 1:3, row = 1:3)
k <- c("1,2" = 1, "1,3" = 2, "2,3" = 3)
mat$lab <- mapply(function(r, c) {
  if (r == c) return("0")
  sprintf("sigma(theta[%d])", k[[paste(min(r, c), max(r, c), sep = ",")]])
}, mat$row, mat$col)
mat$free <- mat$row != mat$col
blocks <- c("X[1]", "X[2]", "X[3]")
f1a <- ggplot(mat, aes(col, -row)) +
  geom_tile(aes(fill = free), colour = "white", linewidth = 1) +
  geom_text(aes(label = lab), parse = TRUE, size = 3) +
  scale_fill_manual(values = c(`FALSE` = "#eeeeee", `TRUE` = "#cde2fb"), guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = parse(text = blocks), position = "top",
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = -(1:3), labels = parse(text = blocks), expand = c(0, 0)) +
  coord_equal() + labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

# (B) the logit map, with the 15 sweep levels (equally spaced in theta)
curve <- data.frame(theta = seq(-L, L, length.out = 300))
curve$C <- plogis(curve$theta)
lev <- data.frame(theta = seq(-L, L, length.out = 15)); lev$C <- plogis(lev$theta)
f1b <- ggplot(curve, aes(theta, C)) +
  annotate("rect", xmin = -L, xmax = L, ymin = 0.1, ymax = 0.9,
           fill = "#f2f2f2") +
  geom_line(linewidth = 0.6) +
  geom_point(data = lev, size = 1.3, colour = COL[["sweep"]]) +
  scale_x_continuous(breaks = c(-6, -3, 0, 3, 6)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = expression(theta), y = expression(C[ij] == sigma(theta)))

# (C) the design cube for three blocks, oblique projection
proj <- function(x, y, z) data.frame(px = x + 0.45 * y, py = z + 0.3 * y)
v <- expand.grid(x = 0:1, y = 0:1, z = 0:1)
edges <- do.call(rbind, lapply(seq_len(nrow(v)), function(i) {
  do.call(rbind, lapply(seq_len(nrow(v)), function(j) {
    if (j <= i || sum(abs(unlist(v[i, ]) - unlist(v[j, ]))) != 1) return(NULL)
    a <- proj(v$x[i], v$y[i], v$z[i]); b <- proj(v$x[j], v$y[j], v$z[j])
    data.frame(x = a$px, y = a$py, xend = b$px, yend = b$py)
  }))
}))
set.seed(1)
pool <- matrix(runif(2000 * 3, -L, L), ncol = 3)
cand <- plogis(pool[caret::maxDissim(matrix(runif(200 * 3, -L, L), ncol = 3),
                                     pool, n = 60), ])
cand_p  <- proj(cand[, 1], cand[, 2], cand[, 3])
sweep_p <- proj(lev$C, lev$C, lev$C)
ends <- rbind(cbind(proj(plogis(-L), plogis(-L), plogis(-L)), lab = "null"),
              cbind(proj(plogis(L), plogis(L), plogis(L)), lab = "full"))
axl <- rbind(cbind(proj(0.5, 0, -0.08), lab = "C[12]"),
             cbind(proj(1.1, 0.5, -0.06), lab = "C[13]"),   # right depth edge
             cbind(proj(-0.07, 0, 0.5), lab = "C[23]"))
f1c <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
               colour = "#bbbbbb", linewidth = 0.3) +
  geom_point(data = cand_p, aes(px, py), colour = COL[["search"]], size = 0.9,
             alpha = 0.8) +
  geom_line(data = sweep_p, aes(px, py), colour = COL[["sweep"]], linewidth = 0.5) +
  geom_point(data = sweep_p, aes(px, py), colour = COL[["sweep"]], size = 1.1) +
  geom_text(data = ends, aes(px, py, label = lab), size = 2.6,
            nudge_x = c(-0.07, 0.07), nudge_y = c(-0.05, 0.05)) +
  geom_text(data = axl, aes(px, py, label = lab), parse = TRUE, size = 2.6) +
  annotate("text", x = 1.45, y = -0.17, hjust = 1, size = 2.4,
           colour = COL[["search"]], label = "search candidates") +
  annotate("text", x = 1.45, y = -0.27, hjust = 1, size = 2.4,
           colour = COL[["sweep"]], label = "sweep (diagonal)") +
  coord_equal() + theme_void()

fig1 <- (f1a | f1b | f1c) + plot_layout(widths = c(0.8, 1, 1.2)) +
  plot_annotation(tag_levels = "A")
save_fig(fig1, "Fig1", 2.4)

## ---- Fig 2: simulation, procedures minus the sweep --------------------------------
f2 <- do.call(rbind, lapply(split(sim, list(sim$scen, sim$n_train), drop = TRUE),
                            function(g) {
                              s <- mean_se(diff_vs_sweep(g, c("search", "search_1se", "oracle"), "test_ber"))
                              cbind(scen = g$scen[1], n_train = g$n_train[1], arm = rownames(s), s)
                            }))
f2$arm <- arm_factor(f2$arm, c("search", "search_1se", "oracle"))
fig2 <- ggplot(f2, aes(factor(n_train), m, colour = arm)) +
  geom_hline(yintercept = 0, colour = COL[["sweep"]], linewidth = 0.4) +
  geom_pointrange(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se),
                  position = position_dodge(width = 0.6), size = 0.25,
                  linewidth = 0.5) +
  facet_wrap(~scen) +
  scale_colour_manual(values = unname(COL[c("search", "search_1se", "oracle")])) +
  labs(x = "Training samples", y = "Test BER minus sweep") +
  theme(legend.position = "top")
save_fig(fig2, "Fig2", 2.8)

## ---- Fig 3: selection optimism ----------------------------------------------------
sim$opt <- sim$test_ber - sim$inner_ber
o_sim <- aggregate(opt ~ scen + n_train + arm,
                   data = sim[sim$arm %in% c("sweep", "search", "search_1se"), ],
                   FUN = mean)
o_sim$arm <- arm_factor(o_sim$arm, c("sweep", "search", "search_1se"))
f3a <- ggplot(o_sim, aes(factor(n_train), opt, colour = arm)) +
  geom_hline(yintercept = 0, colour = "#bbbbbb", linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.6), size = 1.4) +
  facet_wrap(~scen) +
  scale_colour_manual(values = unname(COL[c("sweep", "search", "search_1se")])) +
  labs(x = "Training samples", y = "Optimism (test minus inner BER)")

r3 <- real[real$dataset != REAL[["arms_mogonet_brca_ncomp4.csv"]] &
             real$arm %in% c("sweep", "search", "grid"), ]
r3$opt <- r3$outer_ber - r3$inner_ber
o_real <- aggregate(cbind(opt, n_fits) ~ dataset + arm, data = r3, FUN = mean)
o_real$cand <- o_real$n_fits / 5          # fits per outer fold / inner folds
o_real$dataset <- sub(", ncomp 2", "", as.character(o_real$dataset))
o_real <- o_real[order(o_real$dataset, o_real$cand), ]
lab3 <- o_real[o_real$arm == "grid", ]
f3b <- ggplot(o_real, aes(cand, opt)) +
  geom_hline(yintercept = 0, colour = "#bbbbbb", linewidth = 0.3) +
  geom_line(aes(group = dataset), colour = "#cccccc", linewidth = 0.4) +
  geom_point(aes(colour = arm, shape = dataset), size = 1.8) +
  geom_text(data = lab3, aes(label = dataset), hjust = -0.15, size = 2.4) +
  scale_x_log10(limits = c(10, 2000), breaks = c(15, 63, 216, 343)) +
  scale_colour_manual(values = COL[c("sweep", "search", "grid")],
                      labels = LAB[c("sweep", "search", "grid")]) +
  scale_shape_manual(values = c(16, 15), guide = "none") +
  labs(x = "Candidates per outer fold (log scale)", y = NULL)
fig3 <- (f3a | f3b) + plot_layout(widths = c(1.8, 1)) +
  plot_annotation(tag_levels = "A") & theme(legend.position = "top")
save_fig(fig3, "Fig3", 2.8)

## ---- Fig 4: real data -----------------------------------------------------------
# (A) paired per-fold differences. The full design is left out (breast.TCGA:
# +0.106, it would flatten the panel); its values are in Table W.
arms4 <- c("search", "search_unit", "grid", "empirical", "null")
d4 <- do.call(rbind, lapply(split(real, real$dataset), function(g)
  cbind(dataset = g$dataset[1], diff_vs_sweep(g, arms4, "outer_ber"))))
d4$arm <- arm_factor(d4$arm, arms4)
m4 <- aggregate(diff ~ dataset + arm, d4, mean)
f4a <- ggplot(d4, aes(arm, diff)) +
  geom_hline(yintercept = 0, colour = COL[["sweep"]], linewidth = 0.4) +
  geom_point(aes(colour = arm), size = 0.9, alpha = 0.7,
             position = position_jitter(width = 0.18, height = 0, seed = 1)) +
  geom_errorbar(data = m4, aes(ymin = diff, ymax = diff), width = 0.5,
                linewidth = 0.5) +
  facet_wrap(~dataset, scales = "free") +
  scale_colour_manual(values = unname(COL[arms4]), guide = "none") +
  labs(x = NULL, y = "Outer BER minus sweep") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

# (B) selected design entries, on the logit scale the search works in
sel <- real[real$arm %in% c("sweep", "search", "grid"), ]
long <- do.call(rbind, lapply(pair_cols, function(p) {
  x <- sel[!is.na(sel[[p]]), c("dataset", "arm")]
  if (!nrow(x)) return(NULL)
  cbind(x, pair = sub("-", "\u2013\n", p), C = sel[[p]][!is.na(sel[[p]])])
}))
long$arm <- arm_factor(long$arm, c("sweep", "search", "grid"))
brk <- c(0.0025, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.9975)
f4b <- ggplot(long, aes(pair, qlogis(C), colour = arm)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0,
                                             dodge.width = 0.75, seed = 1),
             size = 0.9, alpha = 0.7) +
  facet_wrap(~dataset, scales = "free_x") +
  scale_y_continuous(breaks = qlogis(brk), labels = brk) +
  scale_colour_manual(values = unname(COL[c("sweep", "search", "grid")])) +
  labs(x = NULL, y = "Selected entry (logit scale)") +
  theme(legend.position = "top")
fig4 <- (f4a / f4b) + plot_annotation(tag_levels = "A")
save_fig(fig4, "Fig4", 6)