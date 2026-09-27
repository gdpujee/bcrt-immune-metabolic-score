# QC Report (QC-001, 2026-09-20)

## GSE42568 (training): n_total=121, tumors w/ OS=104 (events=35, censored=69), normals=17
- OS years median 6.00; censored median 6.74; max 8.28
- Expression: 21755 genes x 121 samples; range [2.31,16.05]
## GSE20685 (validation): n_total=327, valid OS=327 (deaths=83, censored=244)
- OS years median 8.10; censored median 9.20
- Metastasis events 83; Expression 21755 genes x 327 samples
## GSE45827 (biology-only): n_total=155, excluded cell lines=14, bio=141 (tumor=130, normal=11)
- Batches [1, 2, 3, 4, 5, 8, 9, 10]; NO survival (by design); Expression 14522 genes x 155 samples (pre-filtered, 14.5k genes)

## Endpoint harmonization
- Primary OS: 42568 death/censor + days/365.25; 20685 event_death + follow_up_duration_years. Both all-cause OS from diagnosis/sample; units years for KM/Cox.
- Secondary (separate, never pooled): 42568 RFS, 20685 metastasis. Definitions differ -> reported separately.
- Inclusion frozen: 42568 train = 104 tumors w/ OS; 20685 validation = 327 w/ OS; 45827 biology = 141 (130+11) excl. 14 cell lines.
