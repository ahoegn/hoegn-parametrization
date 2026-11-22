[README.md](https://github.com/user-attachments/files/35238525/README.md)
<div align="center">

# Hoegn Parametrization

### Searching the DIABLO design matrix, and knowing when it pays off

[![DOI](https://zenodo.org/badge/DOI/ZENODO_DOI.svg)](https://doi.org/ZENODO_DOI)
[![Preregistration](https://img.shields.io/badge/OSF-preregistered-blue)](https://doi.org/10.17605/OSF.IO/EVTHG)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![R](https://img.shields.io/badge/R-4.6.1-276DC3?logo=r)
![mixOmics](https://img.shields.io/badge/mixOmics-6.36.0-orange)

</div>

---

Code and results for:

> **Hoegn A.** *Tuning the design matrix of DIABLO: a reparametrization for searching block-pair weights, and when it pays off.* (2026, in preparation.)

## The idea in one picture

DIABLO (multi-block sPLS-DA, [mixOmics](https://mixomics.org)) weights every pair of data blocks through a design matrix **C**. In practice, all pairs get the same value, so tuning C means moving along a single line:

```
          full (1,1,1)                         full (1,1,1)
             ●                                    ● · · ·
            /                                   · · · · ·
           /     uniform sweep:               · · · · · ·    logit parametrization:
          /      one value for all pairs      · · · · · ·    every pair has its own weight
         /       (1 of d dimensions)          · · · · · ·    (the whole d-dimensional space)
        ●                                     ● · · ·
   null (0,0,0)                          null (0,0,0)
```

This repository gives each block pair its own coordinate through an elementwise logit,

$$C_{ij} = \sigma(\theta_k), \qquad \sigma(x) = \frac{1}{1+e^{-x}}, \qquad \theta \in \mathbb{R}^{d},\; d = \tfrac{K(K-1)}{2},$$

which turns the design matrix into an ordinary, unconstrained hyperparameter that any optimiser can search, for any objective.

## What we found

| | Result |
|---|---|
| ✅ **Mixed block-pair structure** | The search beat the uniform sweep on held-out BER, by up to 0.029 at n = 400 (better in 45 of 50 simulated datasets). |
| ➖ **Uniform structure** | Search and sweep were equivalent: no false gains. |
| ⚠️ **The cost** | More candidates mean more selection optimism. It shrinks with sample size and with repeated inner cross-validation. |
| 📉 **Real data** | No gain over the sweep on nutrimouse, breast.TCGA and MOGONET BRCA. The PLS-based design recommended by mixOmics was clearly worse than a tuned uniform value (BER 0.150 vs 0.029 on breast.TCGA). |

**Rule of thumb.** Two blocks: sweep the single entry. Three or more blocks: sweep one shared value by default; search the full space when you expect some pairs to need different coupling and have a few hundred samples. Always report the nested, held-out error.

## Repository structure

```
.
├── R/                  analysis scripts, run in numerical order
├── results/
│   ├── csv/            per-fold and per-dataset results (all tables and figures are built from these)
│   └── RDS/            pilot results and grid landscapes
├── supplementary/      Supplementary File S1 (deviations) and Table S2
├── README.md
├── LICENSE             MIT
└── session_info.txt    R and package versions used for all results
```

## Scripts

| Script | What it does | Paper |
|---|---|---|
| `00_data_loader.R` | Loads nutrimouse and breast.TCGA, computes the empirical design, writes `public_datasets.rds` | Methods: Real data |
| `01_parametrization_core.R` | Logit map θ ↔ C, candidate generation (maxDissim), uniform sweep, model fitting with non-convergence counter | Methods: parametrization, search procedures |
| `02_ber_objective.R` | Balanced error rate as a pluggable objective; inner-CV scoring and held-out evaluation | Methods: objective and nested evaluation |
| `03_compare_arms.R` | Nested comparison of sweep, search, grid and fixed designs on a real dataset | Table 4, Figs 3B and 4 |
| `04_sim_generator.R` | Simulated blocks with planted block-pair structure | Methods: Simulation |
| `05_sim_pilot.R` | Fixed designs on simulated data | Table 1 |
| `06_sim_nested.R` | First nested comparison, single inner CV | Supplementary Table S2 |
| `07_sim_samplesize.R` | Main simulation: n = 100, 200, 400, repeated inner CV | Tables 2 and 3, Figs 2 and 3A |
| `08_mogonet_loader.R` | Loads MOGONET BRCA (train and test merged) into `public_datasets.rds` | Methods: Real data |
| `09_figures.R` | Builds Figs 1–4 (TIFF and PDF) from the result files | Figs 1–4 |
| `10_pls_design.R` | PLS-based reference design, scored on the same outer folds | Table 4 (PLS rows) |

## Reproducing the results

**Requirements.** R ≥ 4.6, plus:

```r
install.packages(c("BiocManager", "caret", "patchwork"))
BiocManager::install("mixOmics")
```

**MOGONET BRCA data.** Clone the original repository into `data/`:

```bash
git clone https://github.com/txWang/MOGONET.git data/MOGONET
```

**Run.** From the repository root:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript R/00_data_loader.R
Rscript R/08_mogonet_loader.R
Rscript R/03_compare_arms.R        # set DATASET (and NCOMP_RUN) at the top of the script
Rscript R/05_sim_pilot.R
Rscript R/06_sim_nested.R
Rscript R/07_sim_samplesize.R      # longest run: several hours
Rscript R/10_pls_design.R
Rscript R/09_figures.R
```

The scripts read and write their files in the working directory. The files in `results/` are the outputs used in the paper; to rebuild the figures from them without rerunning the analyses, copy them to the repository root before running `09_figures.R`. `public_datasets.rds` is not included: `00` and `08` regenerate it from the original sources.

The single-threaded BLAS settings matter: the scripts run many small fits in parallel, and a multithreaded BLAS inside each worker slows everything down. All random draws are seeded. Long runs save each finished cell or fold and skip them on restart.

## Preregistration and deviations

The parametrization and evaluation plan were preregistered on OSF on 10 October 2026: [doi.org/10.17605/OSF.IO/EVTHG](https://doi.org/10.17605/OSF.IO/EVTHG). Every addition, change and omission relative to that plan is listed, with reasons, in `supplementary/S1_File.odt`.

## Citation

If you use this code, please cite the paper and this archive:

```bibtex
@software{hoegn_parametrization_2026,
  author  = {Hoegn, Aaron},
  title   = {Hoegn Parametrization: searching the DIABLO design matrix},
  year    = {2026},
  version = {v1.0.0},
  doi     = {ZENODO_DOI},
  url     = {https://doi.org/ZENODO_DOI}
}
```

## License

Code: MIT (see `LICENSE`). Text of the supplementary material: CC BY 4.0.

## AI assistance

Analysis code was written with the assistance of a large language model (Claude, Anthropic) from the author's specifications. The author reviewed the code, ran all analyses and checked the results against the raw output files.
