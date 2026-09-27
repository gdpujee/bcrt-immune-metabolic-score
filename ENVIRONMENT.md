# ENVIRONMENT

**Date**: 2026-09-20
**OS**: Darwin (Apple Git 2.39.5)
**Python**: 3.10.12
**Packages (VERIFIED via pip list / import)**:
- pandas 2.3.3 — VERIFIED
- numpy 2.2.6 — VERIFIED
- scikit-learn 1.7.2 — VERIFIED
- scipy 1.15.3 — VERIFIED
- statsmodels 0.15.0 — VERIFIED (PHReg available, tested import)
- matplotlib 3.10.9 — VERIFIED
- seaborn 0.13.2 — VERIFIED
- requests 2.33.1 — VERIFIED
- lifelines — NOT INSTALLED (pip OSError outside workspace) → NOT_USED, replaced by statsmodels
- R — NOT AVAILABLE on the original 2026-09-20 host; historical note superseded by the 2026-09-26 verification below.

**Network (VERIFIED 2026-09-20)**:
- NCBI GEO HTTPS 200 — VERIFIED
- GEO FTP matrix for GSE42568 (22,565,943 bytes, Last-Modified 2026-07-06) — VERIFIED header-level
- UCSC Xena https 200 — VERIFIED landing reachable

**Randomness**: seed 42 prespecified for all stochastic steps (LASSO CV splits, permutations)
**Genome/annotation**: to be frozen after download (GPL570 annotation + GENCODE/Ensembl for TCGA symbols)

**2026-09-26 PH verification environment (VERIFIED):** Rscript 4.6.1, survival 3.8.6, jsonlite 2.0.0. `Rscript scripts/32_ph_exact.R` generates `results/raw/ph_diagnostics_exact.json` and `logs/ph_exact.log`; this is the exact ranked-time Grambsch–Therneau test. Python statsmodels 0.15.0 remains the Cox/figure runtime; `scripts/28_ph_diagnostics.py` supplies an approximate residual-correlation screen and split-time sensitivity. The early/late split estimates were independently reproduced in R.
