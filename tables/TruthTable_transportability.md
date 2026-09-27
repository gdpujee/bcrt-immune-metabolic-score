<!-- GENERATED FILE — do not edit by hand. Regenerate with `python3 scripts/51_build_truth_table.py`; every value is read from results/raw and cross-checked against the manuscript, supplement, Table 3 and the response letter (TRUTH-001). -->

# Truth table — locked-score validation statistics, all cohorts

| Statistic | GSE20685 | SCAN-B | METABRIC | Source (results/raw/) |
| :--- | :--- | :--- | :--- | :--- |
| Cohort n / deaths | 327 / 83 | 3,273 / 336 | 1,980 / 1,143 | `corrective_summary, rnaseq_summary, metabric_summary` |
| Overall unadjusted HR per SD (95% CI) | 1.59 [1.29,1.98] | 1.44 [1.31,1.59] | 1.13 [1.07,1.20] | `corrective_summary, rnaseq_summary, metabric_summary` |
| Adjusted HR per SD, primary model (95% CI) | 1.60 [1.28,1.99] | 1.46 [1.33,1.61] | 1.08 [1.004,1.168] | `adjusted_per_sd, metabric_adjusted` |
| Missing data: primary / sensitivity | median-imputed n=327; complete case n=325, 83 deaths, HR 1.60 | median-imputed n=3,273; complete case n=2,963, 286 deaths, HR 1.41 | complete case n=1,815, 1,041 deaths; imputation sensitivity n=1,980, HR 1.08 | `adjusted_per_sd, missingdata_sensitivity, metabric_adjusted` |
| C-index (95% CI) | 0.656 [0.60,0.71] | 0.595 [0.56,0.63] | 0.573 [0.554,0.594] | `corrective_summary, rnaseq_summary, metabric_summary` |
| Calibration slope | 0.46 | 0.30 | 0.12 | `adjusted_per_sd` |
| Locked-cutoff High/Low | 48 / 279 | 156 / 3,117 | 1,909 / 71 | `cutoff_splits` |
| Binary HR (95% CI) | 1.95 [1.16,3.29] | 3.12 [2.25,4.32] | 1.00 | `corrective_summary, rnaseq_summary, metabric_summary` |
| cox.zph (Grambsch-Therneau) p, score | 0.0078 | 0.229 | 4.6e-20 | `ph_diagnostics_exact` |
| Split 5y early HR (deaths) | 1.85 [1.42,2.42] (51) | 1.45 [1.31,1.60] (311) | 1.54 [1.41,1.68] (427) | `ph_diagnostics, scanb_split5y` |
| Split 5y late HR (deaths) | 1.22 [0.84,1.75] (32) | 1.42 [0.99,2.03] (25) | 0.92 [0.85,0.99] (716) | `ph_diagnostics, scanb_split5y` |
| Adjusted-model cox.zph (score / global p) | 0.011 / 0.078 | 0.369 / 0.106 | 1.084e-18 / 5.65e-31 | `ph_diagnostics_exact, metabric_adjusted_zph` |
| Time-varying beta (risk x log t) | -0.376 | -0.090 | -0.293 | `ph_timevarying` |
| Time-varying LRT (chi2, p) | 6.53, 0.0106 | 1.42, 0.233 | 77.46, 1.36e-18 | `ph_timevarying` |
| HR(t) 1y -> 10y | 2.52 -> 1.06 | 1.56 -> 1.27 | 1.90 -> 0.97 | `ph_timevarying` |
