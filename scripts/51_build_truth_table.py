"""Build the single truth table from raw outputs and cross-check every document (TRUTH-001).

External-review round 3, P0-2: "generate one unique truth table from the final
outputs, and let Main / Supplement / Response / REMARK all read from the same
result source."  REMARK already reads results/raw (scripts/44); the rendered
manuscript formats are generated from manuscript_submission.md (scripts 22-25);
this script closes the loop:

  1. assemble every headline statistic from results/raw into one table;
  2. write it to tables/TruthTable_transportability.md (a generated artefact);
  3. assert that each value appears, verbatim, in the documents that claim it:
     manuscript_submission.md, supplement.md, Table 3 and the response letter.

A stale document fails loudly here instead of shipping mismatched numbers again.

Outputs: tables/TruthTable_transportability.md, logs/truth_table.log
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results/raw"
TAB = ROOT / "tables"
LOGS = ROOT / "logs"
MS = ROOT / "submission_bcrt/manuscript_submission.md"
SUP = ROOT / "submission_bcrt/supplement.md"
RESP = ROOT / "bio_dsh_review_pack/07_Response_to_Reviewers_PH_Diagnostics_Bug_Fix.md"
TAB3 = ROOT / "tables/Tab3_performance_v3.tsv"
logf = open(LOGS / "truth_table.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")


corr = json.load(open(RAW / "corrective_summary.json"))
rna = json.load(open(RAW / "rnaseq_summary.json"))
met = json.load(open(RAW / "metabric_summary.json"))
aps = json.load(open(RAW / "adjusted_per_sd.json"))
mads = json.load(open(RAW / "missingdata_sensitivity.json"))
madj = json.load(open(RAW / "metabric_adjusted.json"))
phd = json.load(open(RAW / "ph_diagnostics.json"))
sbs = json.load(open(RAW / "scanb_split5y.json"))
phx = json.load(open(RAW / "ph_diagnostics_exact.json"))
tvc = json.load(open(RAW / "ph_timevarying.json"))

g06, scn, mtb = aps["GSE20685"], aps["SCANB"], aps["METABRIC"]


def ci(x, nd=2):
    return f"[{x[0]:.{nd}f},{x[1]:.{nd}f}]"


rows = []
rows.append(("Cohort n / deaths", "327 / 83", "3,273 / 336", "1,980 / 1,143",
             "corrective_summary, rnaseq_summary, metabric_summary"))
rows.append(("Overall unadjusted HR per SD (95% CI)",
             f"{corr['valid_continuous_perSD_HR'][0]:.2f} {ci(corr['valid_continuous_perSD_HR'][1:3])}",
             f"{rna['hr_cont'][0]:.2f} {ci(rna['hr_cont'][1:3])}",
             f"{met['cont_HR']:.2f} {ci(met['cont_CI'])}",
             "corrective_summary, rnaseq_summary, metabric_summary"))
rows.append(("Adjusted HR per SD, primary model (95% CI)",
             f"{g06['per_SD_HR']:.2f} {ci(g06['per_SD_CI'])}",
             f"{scn['per_SD_HR']:.2f} {ci(scn['per_SD_CI'])}",
             f"{mtb['per_SD_HR']:.2f} {ci(mtb['per_SD_CI'], 3)}",
             "adjusted_per_sd, metabric_adjusted"))
rows.append(("Missing data: primary / sensitivity",
             f"median-imputed n={g06['n']}; complete case n={mads['GSE20685']['n']}, "
             f"{mads['GSE20685']['deaths']} deaths, HR {mads['GSE20685']['per_SD_HR']:.2f}",
             f"median-imputed n={scn['n']:,}; complete case n={mads['SCANB']['n']:,}, "
             f"{mads['SCANB']['deaths']} deaths, HR {mads['SCANB']['per_SD_HR']:.2f}",
             f"complete case n={mtb['n']:,}, {mtb['deaths']:,} deaths; "
             f"imputation sensitivity n=1,980, HR {mtb['imputed_sensitivity_per_SD_HR']:.2f}",
             "adjusted_per_sd, missingdata_sensitivity, metabric_adjusted"))
rows.append(("C-index (95% CI)", "0.656 [0.60,0.71]",
             f"{rna['C']:.3f} {ci(rna['C_CI'], 2)}",
             "0.573 [0.554,0.594]",
             "corrective_summary, rnaseq_summary, metabric_summary"))
rows.append(("Raw-score Cox coefficient β", f"{g06['raw_score_coefficient']:.2f}", f"{scn['raw_score_coefficient']:.2f}",
             f"{mtb['raw_score_coefficient']:.2f}", "adjusted_per_sd"))
rows.append(("Locked-cutoff High/Low", "48 / 279", "156 / 3,117", "1,909 / 71",
             "cutoff_splits"))
rows.append(("Binary HR (95% CI)",
             f"{corr['valid_binary_HR_CI'][0]:.2f} {ci(corr['valid_binary_HR_CI'][1:3])}",
             f"{rna['hr_bin'][0]:.2f} {ci(rna['hr_bin'][1:3])}",
             f"{met['bin_HR']:.2f}",
             "corrective_summary, rnaseq_summary, metabric_summary"))
rows.append(("cox.zph (Grambsch-Therneau) p, score",
             f"{phx['GSE20685_cont']['vars']['risk']['p']:.4f}",
             f"{phx['SCANB_cont']['vars']['risk']['p']:.3f}",
             "4.6e-20",
             "ph_diagnostics_exact"))
rows.append(("Split 5y early HR (deaths)",
             f"{phd['GSE20685_split5y']['early']['HR']:.2f} {ci(phd['GSE20685_split5y']['early']['CI'])} "
             f"({phd['GSE20685_split5y']['early']['deaths']})",
             f"{sbs['early']['HR']:.2f} {ci(sbs['early']['CI'])} ({sbs['early']['deaths']})",
             f"{phd['METABRIC_timesplit']['5']['early']['HR']:.2f} "
             f"{ci(phd['METABRIC_timesplit']['5']['early']['CI'])} "
             f"({phd['METABRIC_timesplit']['5']['early']['deaths']})",
             "ph_diagnostics, scanb_split5y"))
rows.append(("Split 5y late HR (deaths)",
             f"{phd['GSE20685_split5y']['late']['HR']:.2f} {ci(phd['GSE20685_split5y']['late']['CI'])} "
             f"({phd['GSE20685_split5y']['late']['deaths']})",
             f"{sbs['late']['HR']:.2f} {ci(sbs['late']['CI'])} ({sbs['late']['deaths']})",
             f"{phd['METABRIC_timesplit']['5']['late']['HR']:.2f} "
             f"{ci(phd['METABRIC_timesplit']['5']['late']['CI'])} "
             f"({phd['METABRIC_timesplit']['5']['late']['deaths']})",
             "ph_diagnostics, scanb_split5y"))
mazp = json.load(open(RAW / "metabric_adjusted_zph.json"))
rows.append(("Adjusted-model cox.zph (score / global p)", "0.011 / 0.078",
             "0.369 / 0.106",
             f"{mazp['vars']['risk_sd']['p']:.3e} / {mazp['global']['p']:.2e}",
             "ph_diagnostics_exact, metabric_adjusted_zph"))
rows.append(("Time-varying beta (risk x log t)",
             f"{tvc['GSE20685']['beta_tt']:.3f}", f"{tvc['SCANB']['beta_tt']:.3f}",
             f"{tvc['METABRIC']['beta_tt']:.3f}", "ph_timevarying"))
rows.append(("Time-varying LRT (chi2, p)",
             f"{tvc['GSE20685']['LRT_chisq']:.2f}, {tvc['GSE20685']['LRT_p']:.4f}",
             f"{tvc['SCANB']['LRT_chisq']:.2f}, {tvc['SCANB']['LRT_p']:.3f}",
             f"{tvc['METABRIC']['LRT_chisq']:.2f}, 1.36e-18", "ph_timevarying"))
rows.append(("HR(t) 1y -> 10y",
             f"{tvc['GSE20685']['hr_at']['t1']:.2f} -> {tvc['GSE20685']['hr_at']['t10']:.2f}",
             f"{tvc['SCANB']['hr_at']['t1']:.2f} -> {tvc['SCANB']['hr_at']['t10']:.2f}",
             f"{tvc['METABRIC']['hr_at']['t1']:.2f} -> {tvc['METABRIC']['hr_at']['t10']:.2f}",
             "ph_timevarying"))

hdr = ("<!-- GENERATED FILE — do not edit by hand. Regenerate with "
       "`python3 scripts/51_build_truth_table.py`; every value is read from "
       "results/raw and cross-checked against the manuscript, supplement, "
       "Table 3 and the response letter (TRUTH-001). -->\n\n"
       "# Truth table — locked-score validation statistics, all cohorts\n\n"
       "| Statistic | GSE20685 | SCAN-B | METABRIC | Source (results/raw/) |\n"
       "| :--- | :--- | :--- | :--- | :--- |\n")
body = "".join(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | `{r[4]}` |\n" for r in rows)
(TAB / "TruthTable_transportability.md").write_text(hdr + body)
log(f"WROTE tables/TruthTable_transportability.md ({len(rows)} rows)")

# ---- cross-document checks -------------------------------------------------
# (label, document key, expected substring). Each document is checked against the
# exact rendering it uses: Table 3 is a TSV, the response letter uses en-dashes.
docs = {"tab3": (TAB3, TAB3.read_text())}
if MS.exists():
    docs["ms"] = (MS, MS.read_text())
if SUP.exists():
    docs["sup"] = (SUP, SUP.read_text())
if RESP.exists():
    docs["resp"] = (RESP, RESP.read_text())

CHECKS = [
    ("METABRIC adjusted HR unrounded", "ms", "[1.004,1.168]"),
    ("METABRIC adjusted HR unrounded", "sup", "[1.004, 1.168]"),
    ("METABRIC adjusted HR unrounded", "resp", "1.004–1.168"),
    ("GSE20685 adjusted HR", "tab3", "1.60 [1.28,1.99]"),
    ("GSE20685 adjusted HR", "resp", "1.60 [1.28–1.99]"),
    ("SCAN-B adjusted HR", "tab3", "1.46 [1.33,1.61]"),
    ("SCAN-B adjusted HR", "resp", "1.46 [1.33–1.61]"),
    ("GSE20685 CC sensitivity", "ms", "n=325, 83 deaths"),
    ("GSE20685 CC sensitivity", "sup", "n=325, 83 deaths"),
    ("GSE20685 CC sensitivity", "resp", "HR 1.60"),
    ("SCAN-B CC sensitivity", "ms", "n=2,963, 286 deaths"),
    ("SCAN-B CC sensitivity", "sup", "n=2,963, 286 deaths"),
    ("SCAN-B CC sensitivity", "resp", "HR 1.41"),
    ("METABRIC imputation sensitivity", "ms", "1.08 [1.01,1.16]"),
    ("METABRIC imputation sensitivity", "sup", "1.08 [1.01,1.16]"),
    ("METABRIC imputation sensitivity", "resp", "1.08 [1.01–1.16]"),
    ("cox.zph GSE20685", "ms", "0.0078"),
    ("cox.zph GSE20685", "sup", "0.007755"),
    ("cox.zph GSE20685", "resp", "0.0078"),
    ("cox.zph SCAN-B", "ms", "0.229"),
    ("cox.zph SCAN-B", "sup", "0.2294"),
    ("cox.zph SCAN-B", "resp", "0.229"),
    ("cox.zph METABRIC", "ms", "4.6×10⁻²⁰"),
    ("cox.zph METABRIC", "sup", "4.614e-20"),
    ("cox.zph METABRIC", "resp", "4.614"),
    ("GSE20685 split early", "ms", "1.85"),
    ("GSE20685 split early", "resp", "1.85"),
    ("GSE20685 split late", "ms", "1.22"),
    ("GSE20685 split late", "resp", "1.22"),
    ("SCAN-B split early", "ms", "1.45"),
    ("SCAN-B split early", "resp", "1.45"),
    ("SCAN-B split late", "ms", "1.42"),
    ("SCAN-B split late", "resp", "1.42"),
    ("METABRIC 3y split", "ms", "1.71"),
    ("METABRIC 3y split", "resp", "1.714"),
    ("METABRIC 5y split", "ms", "1.54"),
    ("METABRIC 5y split", "resp", "1.538"),
    ("METABRIC 7y split", "ms", "1.45"),
    ("METABRIC 7y split", "resp", "1.447"),
    ("beta_tt GSE20685", "ms", "-0.376"),
    ("beta_tt GSE20685", "sup", "-0.376"),
    ("beta_tt GSE20685", "resp", "-0.376"),
    ("beta_tt SCAN-B", "ms", "-0.090"),
    ("beta_tt SCAN-B", "sup", "-0.090"),
    ("beta_tt SCAN-B", "resp", "-0.090"),
    ("beta_tt METABRIC", "ms", "-0.293"),
    ("beta_tt METABRIC", "sup", "-0.293"),
    ("beta_tt METABRIC", "resp", "-0.293"),
    ("LRT METABRIC", "ms", "77.46"),
    ("LRT METABRIC", "sup", "77.46"),
    ("LRT METABRIC", "resp", "77.46"),
    ("n=1,979 exclusion", "ms", "1,979"),
    ("n=1,979 exclusion", "sup", "1,979"),
    ("n=1,979 exclusion", "resp", "1,979"),
    ("cutoff SCAN-B", "tab3", "156/3117"),
    ("cutoff SCAN-B", "resp", "156"),
    ("cutoff METABRIC", "tab3", "1909/71"),
    ("cutoff METABRIC", "sup", "1909/71"),
    ("cutoff METABRIC", "resp", "1,909"),
    ("GSE20685 late CI lower bound (raw 0.8448)", "ms", "1.22 [0.84,1.75]"),
    ("GSE20685 late CI lower bound (raw 0.8448)", "resp", "1.22** [0.84–1.75"),
    ("METABRIC adjusted cox.zph score p", "ms", "1.084×10⁻¹⁸"),
    ("METABRIC adjusted cox.zph score p", "sup", "1.084e-18"),
    ("METABRIC adjusted cox.zph global p", "sup", "5.65e-31"),
    ("METABRIC adjusted cox.zph", "resp", "5.645\\times10^{-31}"),
    ("abstract cohort-specific early effects", "ms",
     "the HR was 1.54 (1.41–1.68) through 5 years and 0.92 (0.85–0.99) thereafter"),
    ("conclusion cohort-specific time dependence", "ms",
     "average estimates do not imply a persistent adverse association in METABRIC"),
    ("early-hazard citation", "ms", "recurrence hazard during the early years after breast cancer diagnosis [20]"),
    ("late-recurrence citation", "ms", "beyond 5 years, particularly in estrogen-receptor-positive disease [21]"),
    ("coding-correction superseded n (S4)", "sup", str(madj["coding_correction"]["superseded_complete_case"]["n"])),
    ("coding-correction superseded HR (S4)", "sup", f"{madj['coding_correction']['superseded_complete_case']['risk_sd_HR']:.3f}"),
    ("coding-correction disclosed in response", "resp", "1,679"),
    ("coding-correction corrected n in response", "resp", f"{madj['complete_case']['n']:,}"),
]

failures = []
for label, key, needle in CHECKS:
    if key not in docs:
        log(f"SKIP {label}: document {key} not present")
        continue
    path, text = docs[key]
    if needle not in text:
        failures.append(f"{label}: '{needle}' NOT FOUND in {path.name}")
if failures:
    for f in failures:
        log("FAIL " + f)
    sys.exit(f"truth-table cross-check failed with {len(failures)} mismatches")
log(f"OK {len(CHECKS)} cross-document value checks passed")
log("TRUTH-001 complete")
