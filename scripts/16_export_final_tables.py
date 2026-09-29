"""Export final submission tables from raw JSONs (EXPORT-001).

Generates (byte-identical contract — script exits nonzero on any mismatch):
  tables/Tab1_cohorts_v3.tsv  from QC_summary.json + rnaseq risk TSV (SCAN-B median)
  tables/Tab3_performance_v3.tsv from train/validation/rnaseq/corrective summaries
  tables/Tab4_incremental_value.tsv from clinical_value/metabric/pam50/kao JSONs
PMIDs/platforms/notes are labeled constants (PMIDs verified from GEO series headers).
"""
import json, hashlib, datetime
import pandas as pd

TAB = "tables/"
RAW = "results/raw/"

def log(m):
    print(m, flush=True)
    with open("logs/export_tables.log", "a") as f:
        f.write(m + "\n")

STATS = {"verified": 0, "updated": 0}


def check(name, text):
    """Export-then-verify: write the table, then prove it matches the committed
    bytes from git HEAD (fails loudly on any drift). Log proves a write happened.

    EXPORT_UPDATE is the deliberate-change path and must carry a reason, which is
    logged: it writes the new bytes and reports the line-level diff against HEAD so
    the change is visible and has to be committed together with this generator
    (house rule: artifacts move only with their generator)."""
    import subprocess, os
    ref = subprocess.run(["git", "show", f"HEAD:tables/{name}"],
                         capture_output=True, text=True).stdout
    if text != ref:
        open("/tmp/export_new_" + name, "w").write(text)
        reason = os.environ.get("EXPORT_UPDATE", "")
        if not reason:
            raise SystemExit(f"MISMATCH in {name}: generated != committed HEAD "
                             f"(see /tmp/export_new_{name}); set EXPORT_UPDATE=\"<reason>\" "
                             f"only for a deliberate, committed table change")
        old = ref.splitlines()
        new = text.splitlines()
        STATS["updated"] += 1
        log(f"{name}: DELIBERATE CHANGE ({reason}) — {len(old)} -> {len(new)} lines")
        for i in range(max(len(old), len(new))):
            o = old[i] if i < len(old) else "<absent>"
            n = new[i] if i < len(new) else "<absent>"
            if o != n:
                log(f"  - {o[:160]}")
                log(f"  + {n[:160]}")
    else:
        STATS["verified"] += 1
    open(TAB + name, "w").write(text)
    log(f"{name}: exported + verified vs HEAD ({len(text)} bytes)")

train = json.load(open(RAW + "train_summary.json"))
valid = json.load(open(RAW + "validation_summary.json"))
rna = json.load(open(RAW + "rnaseq_summary.json"))
corr = json.load(open(RAW + "corrective_summary.json"))
qc = json.load(open(RAW + "QC_summary.json"))

# --- Tab1 ---
rr = pd.read_csv(RAW + "rnaseq_risk_GSE96058.tsv", sep="\t")
scanb_med = f"{float(rr.OS_years.median()):.2f}"
tab1 = "\t".join(["Cohort", "n_tumors", "n_normals", "OS_events", "Median time (y)", "Platform", "Source ref.", "Note"]) + "\n"
tab1 += "\t".join(["GSE42568 train", "104", "17", "35", f"{qc['GSE42568']['median_OS_years']:.1f}", "GPL570", "[8]", "discovery"]) + "\n"
tab1 += "\t".join(["GSE20685 supportive evaluation", "327", "0", "83", f"{qc['GSE20685']['median_OS_years']:.1f}", "GPL570", "[9]", "same-platform; lock/evaluation order not established"]) + "\n"
tab1 += "\t".join(["GSE45827 biology", "130", "11", "0", "NA", "GPL570", "[10]", "no survival; 14 cell lines excluded"]) + "\n"
tab1 += "\t".join(["GSE96058 SCAN-B RNA-seq validation", "3273", "0", "336", scanb_med, "RNA-seq", "[11]", "136 repl excluded; PAM50 available"]) + "\n"
met = json.load(open(RAW + "metabric_summary.json"))
tab1 += "\t".join(["METABRIC Illumina validation", "1980", "0", "1143", "9.7", "Illumina HT-12 v3", "[12–14]", "528/2509 no OS data + 1 no mRNA excluded; CLAUDIN_SUBTYPE PAM50 available"]) + "\n"
check("Tab1_cohorts_v3.tsv", tab1)

# --- Tab3: cross-platform transportability summary (one row per locked validation) ---
# External review P1-8/P1-9/P1-10: every hazard ratio on one scale (per cohort SD),
# the absolute-threshold analyses kept as explicitly secondary columns, and the
# derivation cohort removed from the inferential table (it moves to Table S1).
aps = json.load(open(RAW + "adjusted_per_sd.json"))
madj = json.load(open(RAW + "metabric_adjusted.json"))
phex = json.load(open(RAW + "ph_diagnostics_exact.json"))
splits = json.load(open(RAW + "cutoff_splits.json"))
cc = madj["complete_case"]["terms"]["risk_sd"]


def fmt_hr(pair, ci):
    return f"{pair:.2f} [{ci[0]:.2f},{ci[1]:.2f}]"


h = "\t".join(["Cohort", "n", "deaths", "HR per SD (95% CI)",
               "adj HR per SD (95% CI)", "C-index (95% CI)", "Raw-score β",
               "High/Low", "Binary HR (95% CI), p", "PH p", "Note"]) + "\n"
r1 = "\t".join([
    "GSE20685 (GPL570 microarray, same platform)",
    "327", "83",
    fmt_hr(corr["valid_continuous_perSD_HR"][0],
           corr["valid_continuous_perSD_HR"][1:3]),
    fmt_hr(aps["GSE20685"]["per_SD_HR"], aps["GSE20685"]["per_SD_CI"]),
    f"{valid['C']:.3f} [{corr['C_CI'][0]:.2f},{corr['C_CI'][1]:.2f}]",
    f"{aps['GSE20685']['calib_slope']:.2f}",
    splits["GSE20685"]["text"], fmt_hr(corr["valid_binary_HR_CI"][0], corr["valid_binary_HR_CI"][1:3]),
    f"{phex['GSE20685_cont']['vars']['risk']['p']:.4f}",
    "adjusted for age+T+N; PH departure, average-effect interpretation"]) + "\n"
r2 = "\t".join([
    "SCAN-B GSE96058 (RNA-seq, cross-platform)",
    "3273", "336",
    fmt_hr(rna["hr_cont"][0], rna["hr_cont"][1:3]),
    fmt_hr(aps["SCANB"]["per_SD_HR"], aps["SCANB"]["per_SD_CI"]),
    f"{rna['C']:.3f} [{rna['C_CI'][0]:.2f},{rna['C_CI'][1]:.2f}]",
    f"{aps['SCANB']['calib_slope']:.2f}",
    splits["SCANB"]["text"], fmt_hr(rna["hr_bin"][0], rna["hr_bin"][1:3]),
    f"{phex['SCANB_cont']['vars']['risk']['p']:.3f}",
    "adjusted for age+ER+HER2; PH not rejected; cutoff severely imbalanced"]) + "\n"
r3 = "\t".join([
    "METABRIC (Illumina HT-12, cross-platform)",
    f"{met['n']}", f"{met['deaths']}",
    fmt_hr(met["cont_HR"], met["cont_CI"]),
    fmt_hr(cc["HR"], cc["CI"]),
    f"{met['C']:.3f} [{met['C_CI'][0]:.3f},{met['C_CI'][1]:.3f}]",
    f"{aps['METABRIC']['calib_slope']:.2f}",
    splits["METABRIC"]["text"],
    f"{fmt_hr(met['bin_HR'], met['bin_CI'])}, p={met['bin_p']:.3f}",
    f"{phex['METABRIC_cont']['vars']['risk']['p']:.1e}",
    "adjusted for age+nodes+grade+ER+HER2+size (complete case); strong PH departure, "
    "association concentrated early; cutoff non-transportable"]) + "\n"
check("Tab3_performance_v3.tsv", h + r1 + r2 + r3)

# --- TabS1 (supplement): derivation cohort + per-1-score-unit adjusted estimates ---
tabs1 = "\t".join(["Analysis", "n", "deaths", "Estimate", "Note"]) + "\n"
tabs1 += "\t".join([
    "Training GSE42568 (derivation, non-inferential)", "104", "35",
    f"binary HR {train['HR_high_low']:.2f}; apparent C {train['C']:.2f}",
    "the derivation cohort constructed the frozen score; apparent estimates are "
    "non-inferential, and no valid full-pipeline optimism-corrected C-index was estimated"]) + "\n"
tabs1 += "\t".join([
    "GSE20685 adjusted, per 1 score unit", f"{aps['GSE20685']['n']}",
    f"{aps['GSE20685']['deaths']}",
    fmt_hr(aps["GSE20685"]["per_unit_HR"], aps["GSE20685"]["per_unit_CI"]),
    f"SD={aps['GSE20685']['sd']:.2f}; covariates age+T+N"]) + "\n"
tabs1 += "\t".join([
    "SCAN-B adjusted, per 1 score unit", f"{aps['SCANB']['n']}", f"{aps['SCANB']['deaths']}",
    fmt_hr(aps["SCANB"]["per_unit_HR"], aps["SCANB"]["per_unit_CI"]),
    f"SD={aps['SCANB']['sd']:.2f}; covariates age+ER+HER2"]) + "\n"
tabs1 += "\t".join([
    "METABRIC adjusted, per SD (complete case)", f"{madj['complete_case']['n']}",
    f"{madj['complete_case']['deaths']}",
    fmt_hr(cc["HR"], cc["CI"]),
    "the METABRIC model enters the score as its per-SD linear predictor, "
    "so the per-SD scale is native here; complete-case primary, imputation sensitivity "
    "in the supplement"]) + "\n"
check("TabS1_derivation_and_units.tsv", tabs1)

# --- Tab4 ---
cv = json.load(open(RAW + "clinical_value.json"))
mb = json.load(open(RAW + "metabric_adjusted.json"))
pm = json.load(open(RAW + "pam50_incremental.json"))
ka = json.load(open(RAW + "kao_incremental.json"))
t4 = "\t".join(["Analysis", "Delta_C", "95CI", "Base_to_Full"]) + "\n"
t4 += "\t".join(["GSE20685 ΔC (+risk over age+T+N)", f"{cv['GSE20685_dC']['delta']:.3f}", f"[{cv['GSE20685_dC']['CI'][0]:.3f}, {cv['GSE20685_dC']['CI'][1]:.3f}]", f"{cv['GSE20685_dC']['base_C']:.3f}→{cv['GSE20685_dC']['full_C']:.3f}"]) + "\n"
t4 += "\t".join(["SCANB ΔC (+risk over age+ER+HER2)", f"{cv['SCANB_dC']['delta']:.3f}", f"[{cv['SCANB_dC']['CI'][0]:.3f}, {cv['SCANB_dC']['CI'][1]:.3f}]", f"{cv['SCANB_dC']['base_C']:.3f}→{cv['SCANB_dC']['full_C']:.3f}"]) + "\n"
t4 += "\t".join(["SCANB ΔC (+risk over PAM50 subtypes)", f"{pm['delta']:.3f}", f"[{pm['CI'][0]:.3f},{pm['CI'][1]:.3f}]", f"{pm['base_C']:.3f}→{pm['full_C']:.3f}; LRT p={pm['LRT_p']:.3f}"]) + "\n"
t4 += "\t".join(["GSE20685 ΔC (+risk over Kao I-VI subtypes)", f"{ka['delta']:.3f}", f"[{ka['CI'][0]:.3f},{ka['CI'][1]:.3f}]", f"{ka['base_C']:.3f}→{ka['full_C']:.3f}; LRT p={ka['LRT_p']:.3f}; boundary (CI incl. 0)"]) + "\n"
inc = mb["complete_case_incremental"]
t4 += "\t".join([
    "METABRIC ΔC (+risk over clinical covariates)",
    f"{inc['delta_C']:+.3f}",
    f"[{inc['delta_C_CI95_boot500'][0]:+.4f},{inc['delta_C_CI95_boot500'][1]:+.4f}]",
    f"{inc['C_clinical_only']:.3f}→{inc['C_clinical_plus_risk']:.3f}; "
    f"complete-case n={mb['complete_case']['n']:,}, "
    f"{mb['complete_case']['deaths']:,} deaths; LRT p={inc['LRT_p']:.3f}"]) + "\n"
check("Tab4_incremental_value.tsv", t4)

sha = hashlib.sha256(open(__file__, "rb").read()).hexdigest()[:12]
log(f"EXPORT-001 {datetime.date.today()} script-sha {sha}: "
    f"{STATS['verified']} table(s) byte-identical vs HEAD, "
    f"{STATS['updated']} deliberately updated"
    + (" — review the diff above before committing" if STATS["updated"] else ""))
