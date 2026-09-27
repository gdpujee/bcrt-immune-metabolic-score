"""Put every adjusted model on one scale: adjusted HR per SD of the locked score.

External review P1-8 ("unify the statistical unit of HR") and P1-10 (one
cross-platform transportability table) both require the adjusted and unadjusted
hazard ratios to be reported in the same unit.  The manuscript's existing
adjusted models were reported per 1 score unit (GSE20685 SD=1.01, SCAN-B
SD=1.21), which makes the three cohorts non-comparable.

This script refits each cohort's prespecified adjusted model exactly as
documented (GSE20685: age + T + N; SCAN-B: age + ER + HER2; METABRIC: the
prespecified six-covariate set), with the same median/mode imputation used in the
published models, and reports the adjusted HR (95% CI) for a one-SD increase in
the locked score so that it is directly comparable with the primary unadjusted
per-SD estimate.  Reproduction of the previously published per-unit point
estimates is asserted, so this is a change of scale, not of model.

Outputs: results/raw/adjusted_per_sd.json, logs/adjusted_per_sd.log
"""
import json
import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results/raw"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
logf = open(LOGS / "adjusted_per_sd.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


def impute(df, cols, binary=()):
    d = df.copy()
    for c in cols:
        v = pd.to_numeric(d[c], errors="coerce")
        fill = v.mode()[0] if c in binary and not v.mode().empty else v.median()
        d[c] = v.fillna(fill)
    return d


def fit_adj(df, covs, label, sd):
    """Adjusted Cox with the locked score on the per-1-unit scale, then rescale."""
    cols = ["risk_unit"] + covs
    X = np.column_stack([df[c].values.astype(float) for c in cols])
    r = PHReg(df.time.values, X, df.event.values.astype(int)).fit(disp=0)
    b, se = r.params[0], r.bse[0]
    per_unit = float(np.exp(b))
    per_unit_ci = [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))]
    # one SD of the score on the raw (per-unit) scale == multiplying the linear
    # predictor by sd, so log(HR) scales linearly and the CI follows
    k = float(sd)
    out = {"n": int(len(df)), "deaths": int(df.event.sum()), "sd": round(k, 4),
           "per_unit_HR": round(per_unit, 4),
           "per_unit_CI": [round(per_unit_ci[0], 4), round(per_unit_ci[1], 4)],
           "per_unit_p": float(r.pvalues[0]),
           "per_SD_HR": round(float(np.exp(b * k)), 4),
           "per_SD_CI": [round(float(np.exp((b - 1.96 * se) * k)), 4),
                         round(float(np.exp((b + 1.96 * se) * k)), 4)],
           "per_SD_p": float(r.pvalues[0])}
    log(f"{label}: n={out['n']} deaths={out['deaths']} SD={k:.4f} "
        f"per-unit HR={out['per_unit_HR']:.3f} -> per-SD HR={out['per_SD_HR']:.3f} "
        f"[{out['per_SD_CI'][0]:.3f},{out['per_SD_CI'][1]:.3f}] p={out['per_SD_p']:.3g}")
    return out


out = {}

# ---- GSE20685: age + T + N (as published) ----
g = pd.read_csv(META / "GSE20685_clinical_curated.tsv", sep="\t")
gr = pd.read_csv(RAW / "validation_risk_GSE20685.tsv", sep="\t")
g = g.merge(gr[["GSM", "risk"]], on="GSM")
g = g.dropna(subset=["follow_up_duration_years", "event_death", "risk"])
g = g.rename(columns={"follow_up_duration_years": "time", "event_death": "event",
                      "risk": "risk_unit"})
g = impute(g, ["risk_unit", "age_at_diagnosis", "t_stage", "n_stage"])
out["GSE20685"] = fit_adj(g, ["age_at_diagnosis", "t_stage", "n_stage"],
                          "GSE20685 (age+T+N)", sd=g.risk_unit.std())
out["GSE20685"]["covariates"] = "age_at_diagnosis + t_stage + n_stage (median-imputed)"
assert abs(out["GSE20685"]["per_unit_HR"] - 1.5923) < 5e-3, out["GSE20685"]

# ---- SCAN-B: age + ER + HER2 (as published) ----
# `scripts/11_rnaseq_validate.py` fits the published adjusted model on the aligned
# analysis cohort (the 3273 non-replicate samples with OS), coerces all four columns
# with pd.to_numeric and then imputes EVERY column with the cohort median — including
# the two binary ones.  That is reproduced verbatim here, because the point of this
# script is to change the scale, not the model, and the assertion below only means
# something if the fit is the same fit.
rr = pd.read_csv(RAW / "rnaseq_risk_GSE96058.tsv", sep="\t")
cl = pd.read_csv(META / "GSE96058_clinical_raw.tsv", sep="\t")
# er_status / her2_status are recorded as 1/0, NOT as the words positive/negative;
# an earlier revision mapped the strings and silently produced an all-NaN covariate.
for col in ("er_status", "her2_status"):
    vals = set(pd.to_numeric(cl[col], errors="coerce").dropna().unique())
    if not vals <= {0.0, 1.0}:
        raise SystemExit(f"SCAN-B {col} has unexpected values {sorted(vals)[:6]} — "
                         f"the coding assumption behind the published model changed")
merged = cl.merge(rr[["sample", "risk", "OS_years", "OS_event"]],
                  left_on="title", right_on="sample", how="inner")
sc = pd.DataFrame({
    "time": pd.to_numeric(merged.OS_years, errors="coerce"),
    "event": pd.to_numeric(merged.OS_event, errors="coerce"),
    "risk_unit": pd.to_numeric(merged.risk, errors="coerce"),
    "age_at_diagnosis": pd.to_numeric(merged.age_at_diagnosis, errors="coerce"),
    "er_status": pd.to_numeric(merged.er_status, errors="coerce"),
    "her2_status": pd.to_numeric(merged.her2_status, errors="coerce"),
})
sc = sc.dropna(subset=["time", "event", "risk_unit"])
MODEL_COLS = ["risk_unit", "age_at_diagnosis", "er_status", "her2_status"]
sc[MODEL_COLS] = sc[MODEL_COLS].fillna(sc[MODEL_COLS].median())
log(f"SCAN-B cohort: n={len(sc)} deaths={int(sc.event.sum())} "
    f"(expected n=3273, deaths=336)")
if (len(sc), int(sc.event.sum())) != (3273, 336):
    raise SystemExit("SCAN-B analysis cohort does not match the published one")
out["SCANB"] = fit_adj(sc, ["age_at_diagnosis", "er_status", "her2_status"],
                       "SCAN-B (age+ER+HER2)", sd=sc.risk_unit.std())
out["SCANB"]["covariates"] = "age_at_diagnosis + er_status + her2_status (median-imputed)"
out["SCANB"]["imputation"] = ("cohort-median on all four columns, as in "
                              "scripts/11_rnaseq_validate.py")
assert abs(out["SCANB"]["per_unit_HR"] - 1.3677) < 5e-3, out["SCANB"]

# ---- METABRIC: the PRIMARY (complete-case) estimate, matching Table 3 ----
# The transportability table's adjusted column must be the primary estimate, not the
# imputation sensitivity; the two happen to round to the same value here, which is
# exactly why the distinction has to be made explicit rather than left to luck.
adj = json.load(open(RAW / "metabric_adjusted.json"))
ccp = adj["complete_case"]
mi = adj["imputed_sensitivity"]["terms"]["risk_sd"]
ccr = ccp["terms"]["risk_sd"]
out["METABRIC"] = {
    "n": ccp["n"], "deaths": ccp["deaths"],
    "sd": None, "per_unit_HR": None, "per_unit_CI": None, "per_unit_p": None,
    "per_SD_HR": round(ccr["HR"], 4),
    "per_SD_CI": [round(ccr["CI"][0], 4), round(ccr["CI"][1], 4)],
    "per_SD_p": ccr["p"],
    "primary_model": "complete-case, prespecified six-covariate set",
    "imputed_sensitivity_per_SD_HR": round(mi["HR"], 4),
    "imputed_sensitivity_per_SD_CI": [round(mi["CI"][0], 4), round(mi["CI"][1], 4)],
    "covariates": "age + positive nodes + grade + ER + HER2 + tumor size "
                  "(T category; prespecified, complete case)"}
log(f"METABRIC (complete case, per-SD): n={out['METABRIC']['n']} "
    f"deaths={out['METABRIC']['deaths']} HR={out['METABRIC']['per_SD_HR']:.3f} "
    f"[{out['METABRIC']['per_SD_CI'][0]:.3f},{out['METABRIC']['per_SD_CI'][1]:.3f}] "
    f"p={out['METABRIC']['per_SD_p']:.3g}")
log(f"METABRIC imputation sensitivity (per-SD): HR "
    f"{out['METABRIC']['imputed_sensitivity_per_SD_HR']:.3f} "
    f"[{out['METABRIC']['imputed_sensitivity_per_SD_CI'][0]:.3f},"
    f"{out['METABRIC']['imputed_sensitivity_per_SD_CI'][1]:.3f}]")

out["note"] = ("All adjusted hazard ratios are expressed per one SD of the locked "
               "linear predictor within each cohort (the primary analysis unit), so "
               "the adjusted and unadjusted columns of the transportability table are "
               "directly comparable. Per-1-score-unit values are retained in the "
               "supplementary table. The score itself is never refit.")

# ---- METABRIC calibration slope (Cox coefficient of the raw score; 1.0 = perfect) ----
mc = pd.read_csv(META / "METABRIC_clinical_dl.tsv", sep="\t").drop_duplicates("patientId")
mr = pd.read_csv(RAW / "metabric_risk.tsv", sep="\t")
mc = mc.merge(mr, on="patientId")
mt = pd.to_numeric(mc.OS_MONTHS, errors="coerce") / 12
me = mc.OS_STATUS.astype(str).str.contains("DECEASED", na=False).astype(int)
keep = mt.notna() & (mt > 0)
rs = PHReg(mt[keep].values, mc.risk.values[keep].reshape(-1, 1), me[keep].values).fit(disp=0)
out["METABRIC"]["calib_slope"] = round(float(rs.params[0]), 4)
out["METABRIC"]["calib_slope_SE"] = round(float(rs.bse[0]), 4)
log(f"METABRIC calibration slope (Cox coef of raw score, 1.0=perfect): "
    f"{out['METABRIC']['calib_slope']:.4f} SE {out['METABRIC']['calib_slope_SE']:.4f}")
out["GSE20685"]["calib_slope"] = 0.4647
out["GSE20685"]["calib_slope_SE"] = 0.1093
out["SCANB"]["calib_slope"] = 0.3029
out["SCANB"]["calib_slope_SE"] = 0.0406
out["calib_slope_source"] = ("Cox coefficient of the locked score refit in each "
                             "validation cohort (1.0 = perfect calibration); GSE20685 "
                             "and SCAN-B values from results/raw/corrective_summary.json "
                             "and rnaseq_summary.json, METABRIC computed here.")

json.dump(out, open(RAW / "adjusted_per_sd.json", "w"), indent=2)
log("WROTE results/raw/adjusted_per_sd.json (ADJSD-001)")
logf.close()
