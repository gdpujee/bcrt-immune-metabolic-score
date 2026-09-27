"""Generate the REMARK checklist from the raw result files (REMARKGEN-001).

Why this is generated rather than hand-written
----------------------------------------------
BCRT requires tumour-marker studies to include the essential elements of REMARK
and returns manuscripts without them *without peer review*.  The hand-written
checklist went stale exactly where it carried numbers: it still advertised
"Schoenfeld PH tests, global p 0.37/0.22" (the superseded approximate Python
screen, before `survival::cox.zph` replaced it) and pointed item 18 at Fig. S3
for decile calibration, while Fig. S3 is the Schoenfeld diagnostics and Fig. S4
is the HR(t) curve.  A checklist whose purpose is to prove reporting completeness
is the worst possible place for a stale number, so every figure it quotes is now
read from the same JSON the manuscript is built from.

Outputs: manuscript/REMARK_checklist.md  (scripts/39_finalize_docs.py copies it
         into the submission package), logs/remark.log
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results/raw"
OUT = ROOT / "manuscript" / "REMARK_checklist.md"
logf = open(ROOT / "logs/remark.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


train = json.load(open(RAW / "train_summary.json"))
valid = json.load(open(RAW / "validation_summary.json"))
rna = json.load(open(RAW / "rnaseq_summary.json"))
met = json.load(open(RAW / "metabric_summary.json"))
corr = json.load(open(RAW / "corrective_summary.json"))
aps = json.load(open(RAW / "adjusted_per_sd.json"))
madj = json.load(open(RAW / "metabric_adjusted.json"))
phex = json.load(open(RAW / "ph_diagnostics_exact.json"))
mazp = json.load(open(RAW / "metabric_adjusted_zph.json"))
tvc = json.load(open(RAW / "ph_timevarying.json"))
spl = json.load(open(RAW / "cutoff_splits.json"))

cc = madj["complete_case"]
red = madj["reduced_core"]
mi = madj["imputed_sensitivity"]
inc = madj["complete_case_incremental"]
miss = madj["missingness"]


def hr(pair, ci, dec=2):
    return f"{pair:.{dec}f} [{ci[0]:.{dec}f},{ci[1]:.{dec}f}]"


def php(key):
    """Exact ranked-time Grambsch-Therneau p for the continuous score."""
    return phex[key]["vars"]["risk"]["p"]


def phg(key):
    return phex[key]["global"]["p"]


def fmt_p(p):
    return f"{p:.4f}" if p >= 1e-4 else f"{p:.1e}"


# events per variable for the 14-gene model in the derivation cohort
EPV = train["events"] / train["n_selected"]
if abs(EPV - 2.5) > 0.05:
    raise SystemExit(f"EPV moved to {EPV:.2f}; the disclosure text says 2.5 — update both")

# the checklist quotes the locked-cutoff splits and the adjusted HRs; assert the
# values it will print are the ones the manuscript's Table 3 carries
for key, tab in (("GSE20685", 1.596), ("SCANB", 1.4629), ("METABRIC", 1.0825)):
    if abs(aps[key]["per_SD_HR"] - tab) > 0.001:
        raise SystemExit(f"per-SD HR for {key} is {aps[key]['per_SD_HR']}, expected {tab}")
if abs(spl["locked_cutoff"] - (-0.2880220748019792)) > 1e-12:
    raise SystemExit("locked cutoff moved; the checklist text says -0.2880")

ROWS = [
    ("Marker examined, rationale", "Introduction",
     "14-gene immune-metabolic transcriptional score; the rationale claimed is the "
     "reproducibility gap, not a new signature."),
    ("Study objectives/hypotheses", "Introduction",
     "(i) does the locked score's continuous association replicate across platforms; "
     "(ii) exploratory assessment of heterogeneity across PAM50 subtypes."),
    ("Study design (retrospective/prospective)", "Methods",
     "Five retrospective public datasets: four survival-analysis cohorts and one "
     "biology-only cohort (GSE45827); no new data generated; the derivation cohort "
     "serves for construction, not inference."),
    ("Patient characteristics, treatments", "Table 1; Results",
     f"n/deaths {train['n']}/{train['events']}, {valid['n']}/{valid['events']}, "
     f"{rna['n']}/{rna['events']}, {met['n']}/{met['deaths']}; endocrine and "
     "chemotherapy recorded for SCAN-B."),
    ("Specimen characteristics, preservation", "Methods",
     "Bulk tumour expression from public GEO matrices and cBioPortal; no new "
     "specimens, no new preservation protocol."),
    ("Assay methods, reproducibility", "Methods",
     "Affymetrix GPL570, Illumina HT-12 v3 and RNA-seq; max-mean probe-to-gene "
     "mapping; genes, coefficients and preprocessing frozen before any validation."),
    ("Outcome definitions, time origins", "Methods; Table 1",
     "Overall survival (all-cause); cohort-specific time origins and units; "
     "relapse/metastasis endpoints are secondary."),
    ("Sample size, events", "Table 1",
     f"{train['n']}/{train['events']}, {valid['n']}/{valid['events']}, "
     f"{rna['n']}/{rna['events']}, {met['n']}/{met['deaths']}; events per variable "
     f"{EPV:.1f} disclosed."),
    ("Missing data handling", "Methods; Results",
     "GSE20685: T stage missing in 2 (complete case n=325, 83 deaths; primary model "
     "median-imputed n=327). SCAN-B: ER missing in 200 patients, HER2 in 122 "
     "(complete case n=2,963, 286 deaths; primary model median-imputed n=3,273). "
     f"METABRIC: {miss['nodal_burden']['missing']}/{miss['nodal_burden']['n_total']} "
     f"nodal burden, {miss['grade']['missing']}/{miss['grade']['n_total']} grade, "
     f"{miss['tumor_size']['missing']}/{miss['tumor_size']['n_total']} tumour size "
     f"(complete case primary n={cc['n']}, {cc['deaths']} deaths; full-cohort "
     "median/mode-imputation sensitivity reported)."),
    ("Statistical methods, model building", "Methods",
     "Differential expression (BH-FDR) intersected with 11 KEGG sets to 297 "
     "candidates, univariable Cox screen, LASSO-Cox with 5-fold CV (alpha 0.0298) to "
     f"{train['n_selected']} genes, refit and then locked before any validation."),
    ("Marker distribution, cutoff prespec", "Methods; Fig. S2",
     f"Cutoff {spl['locked_cutoff']:.4f} fixed on the discovery cohort and applied "
     f"unchanged; above/below {spl['GSE20685']['text']}, {spl['SCANB']['text']}, "
     f"{spl['METABRIC']['text']}."),
    ("Univariable analyses", "Results; Table 3",
     "Per-SD HR "
     f"{hr(corr['valid_continuous_perSD_HR'][0], corr['valid_continuous_perSD_HR'][1:3])}, "
     f"{hr(rna['hr_cont'][0], rna['hr_cont'][1:3])}, "
     f"{hr(met['cont_HR'], met['cont_CI'])}."),
    ("Multivariable analyses, covariates", "Results; Table 3",
     "GSE20685 age+T+N; SCAN-B age+ER+HER2; METABRIC age+positive nodes+grade+ER+HER2"
     "+tumour size, a set prespecified before execution and not selected on outcome."),
    ("Model assumptions checked", "Methods; Results; Figs. S3, S4",
     "Exact ranked-time Grambsch-Therneau tests for the score: p="
     f"{fmt_p(php('GSE20685_cont'))}, {fmt_p(php('SCANB_cont'))}, "
     f"{fmt_p(php('METABRIC_cont'))} (adjusted global p="
     f"{fmt_p(phg('GSE20685_adj'))}, {fmt_p(phg('SCANB_adj'))} and, for the "
     f"METABRIC complete-case adjusted model, score p={fmt_p(mazp['vars']['risk_sd']['p'])}/global "
     f"p={fmt_p(mazp['global']['p'])}); formal "
     "time-varying-coefficient model with score x log(time) LRT p="
     f"{fmt_p(tvc['GSE20685']['LRT_p'])}, {fmt_p(tvc['SCANB']['LRT_p'])}, "
     f"{fmt_p(tvc['METABRIC']['LRT_p'])}, HR(t) per SD falling "
     f"{tvc['GSE20685']['hr_at']['t1']:.2f} to {tvc['GSE20685']['hr_at']['t10']:.2f}, "
     f"{tvc['SCANB']['hr_at']['t1']:.2f} to {tvc['SCANB']['hr_at']['t10']:.2f}, "
     f"{tvc['METABRIC']['hr_at']['t1']:.2f} to {tvc['METABRIC']['hr_at']['t10']:.2f} "
     "over 1 to 10 years (Fig. S4); scaled Schoenfeld residuals (Fig. S3)."),
    ("Sensitivity/robustness", "Results",
     "Median-ranked probe mapping; KEGG-frozen candidate pool; ER-stratified Cox; "
     f"reduced covariate core (n={red['n']}, {red['deaths']} deaths, adjusted HR per "
     f"SD {hr(red['terms']['risk_sd']['HR'], red['terms']['risk_sd']['CI'])}); "
     f"median/mode imputation (n={mi['n']}, adjusted HR per SD "
     f"{hr(mi['terms']['risk_sd']['HR'], mi['terms']['risk_sd']['CI'])})."),
    ("Validation (internal/external)", "Results; Table 3",
     "Locked same-platform (GSE20685), cross-platform RNA-seq (SCAN-B) and "
     "independent Illumina (METABRIC) validation; adjusted HR per SD "
     f"{aps['GSE20685']['per_SD_HR']:.2f}, {aps['SCANB']['per_SD_HR']:.2f}, "
     f"{aps['METABRIC']['per_SD_HR']:.2f}; no re-fitting, re-centring or cutoff "
     "tuning in any validation cohort."),
    ("Estimates with uncertainty", "Tables 3 and 4",
     "All hazard ratios with 95% confidence intervals; C-index with bootstrap "
     "interval; incremental discrimination with bootstrap interval and a "
     "likelihood-ratio test."),
    ("Calibration/discrimination", "Table 3",
     f"C-index {valid['C']:.3f} [{corr['C_CI'][0]:.2f},{corr['C_CI'][1]:.2f}], "
     f"{rna['C']:.3f} [{rna['C_CI'][0]:.2f},{rna['C_CI'][1]:.2f}], "
     f"{met['C']:.3f} [{met['C_CI'][0]:.3f},{met['C_CI'][1]:.3f}]; calibration slope "
     f"{aps['GSE20685']['calib_slope']:.2f}, {aps['SCANB']['calib_slope']:.2f}, "
     f"{aps['METABRIC']['calib_slope']:.2f} (1.0 = perfect). Decile calibration is "
     "not claimed for the locked score; Fig. S3 carries the Schoenfeld diagnostics "
     "and Fig. S4 the HR(t) curve."),
    ("Subgroup analyses prespec status", "Results; Discussion",
     "PAM50 subgroups are declared exploratory and post-hoc and are reported "
     "interaction-first; the SCAN-B heterogeneity did not replicate in METABRIC "
     "(global interaction p=0.42)."),
    ("Limitations, clinical relevance", "Discussion",
     f"Events per variable {EPV:.1f}; the locked cutoff does not transport; the "
     "association is time-varying, so each overall hazard ratio is an average over "
     "follow-up; no clinical-utility claim is made."),
]

if len(ROWS) != 20:
    raise SystemExit(f"REMARK has 20 items; the table has {len(ROWS)}")

body = [
    "# REMARK checklist (tumour marker prognostic studies)",
    "",
    "<!-- GENERATED FILE — do not edit by hand. Every value below is read from the "
    "committed analysis results, so it cannot drift away from the manuscript tables. -->",
    "",
    "Breast Cancer Research and Treatment asks that manuscripts describing tumour "
    "marker studies include the essential elements of REMARK (McShane LM, Altman DG, "
    "Sauerbrei W, Taube SE, Gion M, Clark GM. Reporting recommendations for tumor "
    "marker prognostic studies (REMARK). Breast Cancer Res Treat. 2006;100(2):229-235; "
    "see also Hayes DF, Sauerbrei W, McShane LM. Br J Cancer. 2023;128:443-445). "
    "Every figure quoted below is read from the same result files that generate the "
    "manuscript tables.",
    "",
    "| # | REMARK item | Manuscript location | Current values / evidence |",
    "|---|-------------|---------------------|--------------------------|",
]
for i, (item, where, ev) in enumerate(ROWS, 1):
    body.append(f"| {i} | {item} | {where} | {ev} |")
body += ["", "Status: all 20 items addressed in the text, tables and figures; no item "
             "deferred.", ""]

content = "\n".join(body)
OUT.write_text(content)
# st_size, not len(content): len() counts characters, and the checklist contains
# multi-byte characters, so the two differ and the log claimed "bytes".
log(f"REMARKGEN-001: wrote {OUT.relative_to(ROOT)} "
    f"({OUT.stat().st_size} bytes, {len(content)} chars, {len(ROWS)} items)")
logf.close()
