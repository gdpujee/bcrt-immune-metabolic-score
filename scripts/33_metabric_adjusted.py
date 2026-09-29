"""METABRIC clinical-adjusted Cox + missingness transparency (METABRIC-ADJ-001).

Addresses the external review item "P1: METABRIC clinical-adjusted Cox is the
most obvious scientific gap": the manuscript previously reported adjusted models
only for GSE20685 (age+T+N) and SCAN-B (age+ER+HER2).

ADJUSTED COVARIATE SET used in this analysis. The materials included in the
third-party review package do not establish when the variables were selected:
    age            AGE_AT_DIAGNOSIS               (continuous, years)
    nodal_burden   LYMPH_NODES_EXAMINED_POSITIVE  (continuous, positive nodes)
    grade          GRADE                          (ordinal 1/2/3)
    er             ER_STATUS                      (Positive vs Negative)
    her2           HER2_STATUS                    (Positive vs Negative)
    tumor_size     TUMOR_SIZE                     (ordinal T category derived from
                                                  the recorded size in millimetres:
                                                  <=20 mm = T1, >20-50 mm = T2,
                                                  >50 mm = T3; the cBioPortal field
                                                  holds millimetres, not a T digit)
The locked 14-gene score is entered as the same per-SD linear predictor used in
the primary unadjusted analysis; it is never refit.

PRIMARY = complete-case adjusted Cox.  SECONDARY = cohort-median imputation
sensitivity.  Per-variable missingness is reported for every covariate.

Outputs: results/raw/metabric_adjusted.json, logs/metabric_adjusted.log,
         metadata/METABRIC_clinical_adj.tsv (patient-level covariates),
         metadata/METABRIC_sample_clinical.tsv (sample-level covariates)
"""
import json
import time
import urllib.request
import numpy as np
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import chi2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
RES = ROOT / "results/raw"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
SEED = 42
np.random.seed(SEED)
BASE = "https://www.cbioportal.org/api"
ST = "brca_metabric"

PAT_ATTRS = ["AGE_AT_DIAGNOSIS", "LYMPH_NODES_EXAMINED_POSITIVE", "NPI",
             "SEX", "COHORT", "CHEMOTHERAPY", "HORMONE_THERAPY",
             "INFERRED_MENOPAUSAL_STATE"]
SAMPLE_ATTRS = ["ER_STATUS", "HER2_STATUS", "PR_STATUS", "GRADE",
                "TUMOR_SIZE", "TUMOR_STAGE"]


logf = open(LOGS / "metabric_adjusted.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


def get(url, tries=5):
    last = None
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (a + 1))
    raise RuntimeError(f"GET failed: {url}: {last}")


def fetch_attr(dtype, att):
    """Paginated single-attribute fetch; returns list of records."""
    out, page = [], 0
    while True:
        d = get(f"{BASE}/studies/{ST}/clinical-data?clinicalDataType={dtype}"
                f"&attributeId={att}&projection=SUMMARY&pageSize=5000"
                f"&pageNumber={page}")
        if not d:
            break
        out += d
        page += 1
        if page > 12:
            raise RuntimeError(f"pagination guard tripped for {att}")
    return out


# ---- Phase 1: clinical covariate acquisition (checkpointed) ----
# NOTE: SAMPLE-level attributes are keyed by patientId, i.e. this assumes one tumour
# sample per patient.  scripts/41_metabric_sample_multiplicity.py verifies that
# assumption against the live API; results/raw/metabric_sample_multiplicity.json.
CK = META / "METABRIC_clinical_adj.tsv"
CKS = META / "METABRIC_sample_clinical.tsv"
if CK.exists() and CKS.exists():
    log("resume: covariate caches present")
else:
    pat = {}
    for att in PAT_ATTRS:
        recs = fetch_attr("PATIENT", att)
        log(f"PATIENT {att}: {len(recs)} records")
        for x in recs:
            pat.setdefault(x["patientId"], {})[att] = x.get("value")
    pd.DataFrame([{"patientId": k, **v} for k, v in pat.items()]).to_csv(
        CK, sep="\t", index=False)

    smp = {}
    for att in SAMPLE_ATTRS:
        recs = fetch_attr("SAMPLE", att)
        log(f"SAMPLE {att}: {len(recs)} records")
        for x in recs:
            smp.setdefault(x["patientId"], {})[att] = x.get("value")
    pd.DataFrame([{"patientId": k, **v} for k, v in smp.items()]).to_csv(
        CKS, sep="\t", index=False)

pc = pd.read_csv(CK, sep="\t")
sc = pd.read_csv(CKS, sep="\t")
log(f"patient covariates: {pc.shape}, sample covariates: {sc.shape}")
for c in pc.columns:
    if c != "patientId":
        log(f"  {c}: missing {int(pc[c].isna().sum())}/{len(pc)}")
for c in sc.columns:
    if c != "patientId":
        log(f"  {c}: missing {int(sc[c].isna().sum())}/{len(sc)}")

# ---- Phase 2: build analysis frame ----
risk = pd.read_csv(RES / "metabric_risk.tsv", sep="\t")
summ = json.load(open(RES / "metabric_summary.json"))
clin = pd.read_csv(META / "METABRIC_clinical_dl.tsv", sep="\t").drop_duplicates("patientId")
# AGE_AT_DIAGNOSIS exists in both the earlier download and this cBioPortal fetch, so
# the merge is suffixed and the two sources are compared instead of one silently
# shadowing the other.
m = (risk.merge(clin, on="patientId", how="inner")
         .merge(pc, on="patientId", how="left", suffixes=("_dl", "_cbio"))
         .merge(sc, on="patientId", how="left"))
if {"AGE_AT_DIAGNOSIS_dl", "AGE_AT_DIAGNOSIS_cbio"} <= set(m.columns):
    a = pd.to_numeric(m.AGE_AT_DIAGNOSIS_dl, errors="coerce")
    b = pd.to_numeric(m.AGE_AT_DIAGNOSIS_cbio, errors="coerce")
    both = a.notna() & b.notna()
    dmax = float((a[both] - b[both]).abs().max()) if int(both.sum()) else 0.0
    log(f"AGE_AT_DIAGNOSIS cross-source check: {int(both.sum())} comparable rows, "
        f"max |difference| = {dmax}")
    if dmax > 1e-9:
        raise SystemExit("AGE_AT_DIAGNOSIS disagrees between the two cBioPortal "
                         "downloads; refusing to pick one arbitrarily")
    m["AGE_AT_DIAGNOSIS"] = b.where(b.notna(), a)
m["OS_years"] = pd.to_numeric(m.OS_MONTHS, errors="coerce") / 12
m["OS_event"] = m.OS_STATUS.str.contains("DECEASED", na=False).astype(int)
m = m.dropna(subset=["OS_years", "OS_event", "risk"]).copy()

# covariate coding
m["age"] = pd.to_numeric(m.AGE_AT_DIAGNOSIS, errors="coerce")
m["nodal_burden"] = pd.to_numeric(m.LYMPH_NODES_EXAMINED_POSITIVE, errors="coerce")
m["grade"] = pd.to_numeric(m.GRADE, errors="coerce")
m["er"] = m.ER_STATUS.map({"Positive": 1, "Negative": 0})
m["her2"] = m.HER2_STATUS.map({"Positive": 1, "Negative": 0})
# TUMOR_SIZE in brca_metabric is the tumour size in MILLIMETRES (median 22.4 mm, max
# 182 mm), not a T-category digit.  An earlier revision of this script took the leading
# digit of the number, which produced a meaningless covariate; the ordinal T category is
# derived from the size instead.  See DECISIONS.md D-039 for the disclosed correction.
_size_mm = pd.to_numeric(m.TUMOR_SIZE, errors="coerce")
m["tumor_size"] = pd.Series(
    np.select([_size_mm <= 20, (_size_mm > 20) & (_size_mm <= 50), _size_mm > 50],
              [1.0, 2.0, 3.0], default=np.nan),
    index=m.index)
m.loc[_size_mm.isna(), "tumor_size"] = np.nan

COVS = ["age", "nodal_burden", "grade", "er", "her2", "tumor_size"]
log(f"analysis frame: n={len(m)} deaths={int(m.OS_event.sum())}")
miss = {}
for c in COVS:
    miss[c] = {"missing": int(m[c].isna().sum()), "n_total": int(len(m)),
               "pct": round(100 * float(m[c].isna().mean()), 2)}
    log(f"  covariate {c}: missing {miss[c]['missing']}/{len(m)} "
        f"({miss[c]['pct']}%)")
log(f"  TUMOR_SIZE raw values: "
    f"{m.TUMOR_SIZE.value_counts(dropna=False).head(8).to_dict()}")

m["risk_sd"] = m.risk / m.risk.std()

# ---- Phase 3: models ----
def fit(df, cols, label):
    X = np.column_stack([df[c].values.astype(float) for c in cols])
    r = PHReg(df.OS_years.values, X, df.OS_event.values.astype(int)).fit(disp=0)
    hr = np.exp(r.params)
    ci = np.column_stack([np.exp(r.params - 1.96 * r.bse),
                          np.exp(r.params + 1.96 * r.bse)])
    out = {"n": int(len(df)), "deaths": int(df.OS_event.sum()),
           "terms": {c: {"HR": round(float(hr[i]), 4),
                         "CI": [round(float(ci[i, 0]), 4), round(float(ci[i, 1]), 4)],
                         "p": float(r.pvalues[i])}
                     for i, c in enumerate(cols)}}
    log(f"{label}: n={out['n']} deaths={out['deaths']} "
        + " ".join(f"{c} HR={out['terms'][c]['HR']:.3f}"
                   f"[{out['terms'][c]['CI'][0]:.3f},{out['terms'][c]['CI'][1]:.3f}]"
                   f" p={out['terms'][c]['p']:.3g}" for c in cols))
    return out


def cindex_fast(T, E, r):
    o = np.argsort(T, kind="stable")
    T, E, r = T[o], E[o], r[o]
    num = den = 0.0
    for i in np.where(E == 1)[0]:
        later = np.arange(i + 1, len(T))
        later = later[T[later] > T[i]]
        if len(later) == 0:
            continue
        den += len(later)
        num += np.sum(r[later] < r[i]) + 0.5 * np.sum(r[later] == r[i])
    return num / den if den else 0.5


def ll(df, cols):
    X = np.column_stack([df[c].values.astype(float) for c in cols])
    return PHReg(df.OS_years.values, X, df.OS_event.values.astype(int)).fit(disp=0).llf


out = {"model_covariates": COVS, "missingness": miss,
       "coding": {"er": "Positive=1, Negative=0", "her2": "Positive=1, Negative=0",
                  "grade": "ordinal 1/2/3 numeric",
                  "tumor_size": "T category from recorded size in mm: <=20 -> 1, "
                                ">20-50 -> 2, >50 -> 3",
                  "nodal_burden": "count of positive nodes",
                  "risk": "locked linear predictor / cohort SD"},
       "coding_correction": {
           "issue": "an earlier revision derived tumor_size by taking the leading digit "
                    "of TUMOR_SIZE, which is recorded in millimetres (median 22.4, max "
                    "182), producing a meaningless covariate",
           "detected": "logged TUMOR_SIZE value counts showed 20.0/25.0/30.0 mm, not "
                       "T-category digits",
           "superseded_complete_case": {"n": 1679, "deaths": 944,
                                        "risk_sd_HR": 1.096, "risk_sd_CI": [1.011, 1.188],
                                        "risk_sd_p": 0.0257,
                                        "tumor_size_HR": 1.213,
                                        "tumor_size_CI": [1.130, 1.302],
                                        "tumor_size_p": 9.34e-08},
           "note": "disclosed rather than silently overwritten; the covariate list "
                   "is unchanged between the superseded and corrected fits"}}
log("tumor_size T-category distribution: "
    f"{m.tumor_size.value_counts(dropna=False).sort_index().to_dict()}")

# PRIMARY: complete-case with the full adjusted set
cc = m.dropna(subset=["risk_sd"] + COVS).copy()
out["complete_case"] = fit(cc, ["risk_sd"] + COVS, "PRIMARY complete-case full set")
# reduced-covariate sensitivity: core clinical (age, grade, er, her2)
red = m.dropna(subset=["risk_sd", "age", "grade", "er", "her2"]).copy()
out["reduced_core"] = fit(red, ["risk_sd", "age", "grade", "er", "her2"],
                          "SENS complete-case core set")

# incremental value: clinical-only vs clinical+risk (same complete-case rows)
for name, df_, cols in (("complete_case", cc, COVS), ("reduced_core", red,
                                                      ["age", "grade", "er", "her2"])):
    # clinical-only linear predictor from a Cox fit
    r0 = PHReg(df_.OS_years.values,
               np.column_stack([df_[c].values.astype(float) for c in cols]),
               df_.OS_event.values.astype(int)).fit(disp=0)
    lp0 = np.column_stack([df_[c].values.astype(float) for c in cols]) @ r0.params
    r1 = PHReg(df_.OS_years.values,
               np.column_stack([df_[c].values.astype(float) for c in cols]
                               + [df_.risk_sd.values]),
               df_.OS_event.values.astype(int)).fit(disp=0)
    lp1 = (np.column_stack([df_[c].values.astype(float) for c in cols]
                           + [df_.risk_sd.values]) @ r1.params)
    C0 = cindex_fast(df_.OS_years.values, df_.OS_event.values.astype(int), lp0)
    C1 = cindex_fast(df_.OS_years.values, df_.OS_event.values.astype(int), lp1)
    lrt = float(2 * (r1.llf - r0.llf))
    p_lrt = float(1 - chi2.cdf(lrt, 1))
    rng = np.random.default_rng(SEED)
    boots = []
    n = len(df_)
    for _ in range(500):
        idx = rng.integers(0, n, n)
        s = df_.iloc[idx]
        try:
            rr0 = PHReg(s.OS_years.values,
                        np.column_stack([s[c].values.astype(float) for c in cols]),
                        s.OS_event.values.astype(int)).fit(disp=0)
            rr1 = PHReg(s.OS_years.values,
                        np.column_stack([s[c].values.astype(float) for c in cols]
                                        + [s.risk_sd.values]),
                        s.OS_event.values.astype(int)).fit(disp=0)
            b0 = cindex_fast(s.OS_years.values, s.OS_event.values.astype(int),
                             np.column_stack([s[c].values.astype(float) for c in cols]) @ rr0.params)
            b1 = cindex_fast(s.OS_years.values, s.OS_event.values.astype(int),
                             np.column_stack([s[c].values.astype(float) for c in cols]
                                             + [s.risk_sd.values]) @ rr1.params)
            boots.append(b1 - b0)
        except Exception:  # noqa: BLE001
            continue
    ci = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
    out[name + "_incremental"] = {
        "C_clinical_only": round(float(C0), 4),
        "C_clinical_plus_risk": round(float(C1), 4),
        "delta_C": round(float(C1 - C0), 4),
        "delta_C_CI95_boot500": [round(ci[0], 4), round(ci[1], 4)],
        "LRT_chi2": round(lrt, 3), "LRT_df": 1, "LRT_p": p_lrt,
        "boot_draws_used": len(boots)}
    log(f"{name} incremental: C {C0:.3f} -> {C1:.3f} (dC {C1-C0:+.4f} "
        f"[{ci[0]:+.4f},{ci[1]:+.4f}]) LRT chi2={lrt:.2f} p={p_lrt:.3g}")

# SECONDARY: imputation sensitivity (cohort median for continuous, mode for binary)
mi = m.copy()
for c in COVS:
    if mi[c].isna().any():
        fill = mi[c].median() if c not in ("er", "her2") else mi[c].mode()[0]
        mi[c] = mi[c].fillna(fill)
out["imputed_sensitivity"] = fit(mi, ["risk_sd"] + COVS, "SENS median/mode-imputed full set")

json.dump(out, open(RES / "metabric_adjusted.json", "w"), indent=2)
log("WROTE results/raw/metabric_adjusted.json (METABRIC-ADJ-001)")
