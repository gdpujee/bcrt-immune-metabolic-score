# Immune–metabolic transcriptional score in breast cancer: analysis code and derived data

Analysis code and derived data tables accompanying the manuscript *What transfers and what does not: cross-platform evaluation of a fixed immune–metabolic transcriptional score for breast cancer overall survival*,
submitted to *Breast Cancer Research and Treatment*.

**Archived software version:** v1.0.2. The version DOI is listed on its Zenodo record.
*(Concept DOI for all versions: [10.5281/zenodo.22994650](https://doi.org/10.5281/zenodo.22994650))*

## What the analysis does

A 14-gene immune–metabolic score was built in GSE42568 (104 tumours, 35 deaths). The recorded score was evaluated in GSE20685, but the model-lock record and initial GSE20685 analysis share a commit, so their order cannot be established; we treat GSE20685 as supportive same-platform evidence. Subsequent commits document application of the fixed parameters without refitting or recentering to SCAN-B and METABRIC. The primary cross-platform evidence therefore comes from these two later evaluations; the absolute cutoff and calibration are reported as separate failure dimensions.

| Cohort | Accession | Assay | Analysable | Deaths | Evidence role |
|---|---|---|---|---|---|
| GSE20685 | GEO | Affymetrix GPL570 | 327 tumours | 83 | Supportive same-platform evaluation; lock order unresolved |
| SCAN-B | GSE96058 | RNA-seq | 3,273 patients | 336 | Subsequent fixed-parameter cross-platform validation |
| METABRIC | cBioPortal `brca_metabric` | Illumina HT-12 | 1,980 patients | 1,143 | Subsequent fixed-parameter cross-platform validation |

The primary estimand is the continuous per-SD association from Cox regression. Same-direction continuous associations were observed in GSE20685 and in the two later fixed-parameter applications (per-SD HR 1.59 [95% CI 1.29–1.98], 1.44 [1.31–1.59] and 1.13 [1.07–1.20]); the absolute cutoff did not transport. The exploratory five-class PAM50 heterogeneity pattern seen in SCAN-B
was not reproduced in the broader METABRIC CLAUDIN_SUBTYPE analysis (global
interaction p=0.42). Effects were time-dependent in the two
microarray cohorts, so those estimates are averages over follow-up.

## Repository layout

| Path | Contents |
|---|---|
| `scripts/` | The whole pipeline, in execution order (see below) |
| `results/raw/` | Per-analysis result files (JSON / TSV) read by the manuscript builders |
| `results/derived/` | The locked model: `locked_model.json`, `locked_coefficients.tsv` |
| `tables/` | The published tables as TSV (`tables/archive/` holds superseded versions, kept for provenance) |
| `figures/` | Figure sources, PDF and PNG |
| `MODEL_LOCK_MANIFEST.md` | The frozen model audit trail: 14 genes, coefficients, normalization, cutoff, file hashes, commit timeline |

## What is not in this repository

The primary expression and clinical data, and the probe-mapped expression matrices
derived from them. All of it is public and is downloaded from the accessions above by
`scripts/02_parse_geo.py`, `scripts/10_parse_scanb.py` and the cBioPortal REST API; the
matrices are left out because a single one exceeds GitHub's 100 MB per-file limit.
The repository holds only what the analysis derives from them at a publishable size. A fresh end-to-end run therefore needs network access to retrieve source data; downstream stages can run offline after inputs have been downloaded and cached locally.

The manuscript, the journal submission package, and the project's internal audit
trail (review notes, decision log, evidence ledger) are not part of this release.

## Reproducing the analysis

### Environment

Python 3.10.12 with the pinned set in `requirements.txt`:

```sh
python3 -m pip install -r requirements.txt
```

`ENVIRONMENT.md` records the verified versions and the network checks. Two
proportional-hazards steps are written in R (`scripts/32_ph_exact.R`,
`scripts/35_ph_timevarying.R`) and need R ≥ 4.6 with `survival` and `jsonlite`; the
exact Grambsch–Therneau test is only available there. Every stochastic step uses
seed 42, fixed in advance.

### Order

Scripts are numbered in execution order and each carries a module docstring stating
its contract. The stages are:

| Stage | Scripts | What they do |
|---|---|---|
| Parse and QC | `02`–`03` | Download and parse the GEO series; QC and endpoints |
| Candidate pool | `04` | Differential expression, KEGG intersection → 297 candidates |
| Training (non-inferential) | `05` | Univariable screen, LASSO-Cox, coefficient refit, cutoff lock |
| Same-platform cohort evaluation | `06`, `07`, `09` | GSE20685, robustness, corrective passes |
| RNA-seq validation | `10`–`11`, `17` | SCAN-B parse and locked validation |
| Cross-platform validation | `27`, `33`, `38`, `41`, `42` | METABRIC, clinically adjusted Cox, one per-SD scale, multiplicity, cutoff splits |
| Proportional hazards | `28`, `32`, `35` | Residual screen, exact Grambsch–Therneau test, time-varying coefficient model |
| Subgroup and sensitivity | `12`–`15`, `18`–`20`, `29`, `31` | Clinical value, PAM50, ER-stratified, probe mapping, gene overlap, bootstrap stability |
| Figures and tables | `08`, `16`, `21`, `30` | Figure and table generation, sample manifest |

The manuscript builders (`22`–`25`, `34`, `36`, `39`, `44`, `46`) and the consistency
gate (`40`) render the submission package. They expect the manuscript sources and a
working tree that is not part of this release, so they are not runnable from a fresh
clone of the archive alone; the analysis stages above are.

### Determinism

Every builder is idempotent: re-running the chain over unchanged inputs reproduces
byte-identical artifacts, including the PDFs (document creation timestamps are
pinned). `MANIFEST.sha256` in the release archive lists a SHA-256 for every file, so
a download can be verified with `sha256sum -c MANIFEST.sha256`.

## Citation

Cite the archived release and the article. The concept DOI is above; the exact
version DOI for the release you used is on the Zenodo record page.
`CITATION.cff` carries machine-readable metadata if your reference manager reads it.

## Licence

MIT — see `LICENSE`. It covers the analysis code and the derived tables authored
here. The primary data are not redistributed and remain under the terms of their
providers (NCBI GEO, cBioPortal); third-party gene sets and reference lists used as
inputs remain under their own licences.

## Contact

Qiang Li — corresponding author, Foshan, Guangdong, China (unaffiliated /
independent researcher). Email: gdpujee@gmail.com, ORCID
[0009-0003-6404-335X](https://orcid.org/0009-0003-6404-335X).

Danhua He — Guangdong Provincial Hospital of Chinese Medicine, Guangzhou, Guangdong,
China. ORCID [0009-0007-0238-2494](https://orcid.org/0009-0007-0238-2494).
