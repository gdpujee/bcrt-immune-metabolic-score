"""Independent METABRIC validation of the locked 14-gene score (METABRIC-001).

Second PAM50-heterogeneity cohort (cBioPortal brca_metabric, Illumina HT-12).
Phase 1 (download, checkpoint-resumed): per-patient clinical
(OS_MONTHS/OS_STATUS/CLAUDIN_SUBTYPE/AGE_AT_DIAGNOSIS/ER_STATUS/HER2_STATUS/GRADE)
+ batched mRNA for the 14 locked genes -> data/raw/METABRIC_* (gitignored).
Phase 2 (analysis): LOCKED genes/coefs/training means/SDs applied once
(no refit, no cutoff tuning); primary continuous per-SD Cox; binary locked
cutoff secondary (skew disclosed); PAM50-stratified + global interaction.
Outputs: results/raw/metabric_summary.json, logs/metabric.log
"""
import json
import time
import urllib.request
import pandas as pd
import numpy as np
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import chi2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
RES = ROOT / "results/raw"
DER = ROOT / "results/derived"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
SEED = 42
np.random.seed(SEED)
BASE = "https://www.cbioportal.org/api"
ST = "brca_metabric"
GENES = {"FBP1": 2203, "ALDH2": 217, "PRPS2": 5634, "TKT": 7086,
         "IDNK": 414328, "UQCRHL": 440567, "GRK6": 2870, "IKBKB": 3551,
         "RAP1B": 5908, "TNFRSF19": 55504, "CXCL14": 9547, "DLG1": 1739,
         "TNFRSF12A": 51330, "PLCG1": 5335}


def log(m):
    print(m, flush=True)
    with open(LOGS / "metabric.log", "a") as f:
        f.write(m + "\n")


def get(url, tries=5):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if a == tries - 1:
                raise
            time.sleep(2 * (a + 1))


def post(url, payload, tries=5):
    import urllib.error
    for a in range(tries):
        try:
            req = urllib.request.Request(
                url, data=json.dumps(payload).encode(),
                headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if a == tries - 1:
                raise
            time.sleep(3 * (a + 1))


# ---- Phase 1a: patient IDs ----
pids, pg = [], 0
while True:
    d = get(f"{BASE}/studies/{ST}/patients?pageSize=500&pageNumber={pg}&projection=ID")
    if not d:
        break
    pids += [x["patientId"] for x in d]
    pg += 1
log(f"METABRIC patients: {len(pids)}")

# ---- Phase 1b: clinical (checkpointed) ----
WANT = ["OS_MONTHS", "OS_STATUS", "CLAUDIN_SUBTYPE", "AGE_AT_DIAGNOSIS",
        "ER_STATUS", "HER2_STATUS", "GRADE"]
ck = META / "METABRIC_clinical_dl.tsv"
done = set()
if ck.exists():
    done = set(pd.read_csv(ck, sep="\t").patientId.tolist())
    log(f"resume: {len(done)} cached")
rows = []
for i, pid in enumerate(pids):
    if pid in done:
        continue
    d = get(f"{BASE}/studies/{ST}/patients/{pid}/clinical-data?projection=SUMMARY&pageSize=100")
    rec = {"patientId": pid}
    for x in d:
        if x.get("clinicalAttributeId") in WANT:
            rec[x["clinicalAttributeId"]] = x.get("value")
    rows.append(rec)
    if len(rows) % 100 == 0:
        pd.DataFrame(rows).to_csv(ck, sep="\t", index=False,
                                  mode="a" if ck.exists() else "w",
                                  header=not ck.exists())
        log(f"clinical {i + 1}/{len(pids)}")
        rows = []
if rows:
    pd.DataFrame(rows).to_csv(ck, sep="\t", index=False,
                              mode="a" if ck.exists() else "w",
                              header=not ck.exists())
clin = pd.read_csv(ck, sep="\t").drop_duplicates("patientId")
log(f"clinical rows: {len(clin)}")
log("CLAUDIN values: " + str(clin.CLAUDIN_SUBTYPE.value_counts(dropna=False).to_dict()))

# ---- Phase 1c: mRNA batch ----
mrna = post(f"{BASE}/molecular-profiles/brca_metabric_mrna/molecular-data/fetch",
            {"entrezGeneIds": list(GENES.values()),
             "sampleListId": "brca_metabric_mrna"})
me = pd.DataFrame(mrna)
inv = {v: k for k, v in GENES.items()}
me["symbol"] = me.entrezGeneId.map(inv)
X = me.pivot(index="symbol", columns="sampleId", values="value")
X.to_csv(RAW / "METABRIC_mrna_14g.tsv", sep="\t")
log(f"mRNA matrix: {X.shape[0]}/14 genes x {X.shape[1]} samples")
open(RAW / "METABRIC_fetch.json", "w").write(json.dumps(
    {"date": "2026-09-26", "study": ST, "profile": "brca_metabric_mrna",
     "genes": GENES}, indent=2))

# ---- Phase 2: locked evaluation ----
locked = json.load(open(DER / "locked_model.json"))
genes, coefs = locked["genes"], np.array(locked["coefs"])
means = pd.Series(locked["scaling"]["means"])
sds = pd.Series(locked["scaling"]["sds"])
cutoff = locked["cutoff"]
m = clin.merge(X.T, left_on="patientId", right_index=True, how="inner")
m["OS_years"] = pd.to_numeric(m.OS_MONTHS, errors="coerce") / 12
m["OS_event"] = m.OS_STATUS.str.contains("DECEASED", na=False).astype(int)
m = m.dropna(subset=["OS_years", "OS_event"] + genes).copy()
Z = (m[genes] - means) / sds
m["risk"] = Z.values @ coefs
T, E = m.OS_years.values, m.OS_event.values.astype(int)
sd = m.risk.std()
rc = PHReg(T, (m.risk.values / sd).reshape(-1, 1), E).fit(disp=0)
hr = float(np.exp(rc.params[0]))
ci = [float(np.exp(rc.params[0] - 1.96 * rc.bse[0])),
      float(np.exp(rc.params[0] + 1.96 * rc.bse[0]))]
hi = (m.risk.values > cutoff).astype(int)
rb = PHReg(T, hi.reshape(-1, 1), E).fit(disp=0)
bin_hr = float(np.exp(rb.params[0]))
bin_ci = [float(np.exp(rb.params[0] - 1.96 * rb.bse[0])),
          float(np.exp(rb.params[0] + 1.96 * rb.bse[0]))]
bin_p = float(rb.pvalues[0])
log(f"PRIMARY continuous per-SD: n={len(m)} deaths={int(E.sum())} "
    f"HR={hr:.3f} [{ci[0]:.3f},{ci[1]:.3f}] p={float(rc.pvalues[0]):.3g}")
log(f"SECONDARY locked cutoff: High={int(hi.sum())}/Low={int((1-hi).sum())} "
    f"HR={bin_hr:.3f} [{bin_ci[0]:.3f},{bin_ci[1]:.3f}] p={bin_p:.3g}")
out = {"n": int(len(m)), "deaths": int(E.sum()), "cont_HR": hr,
       "cont_CI": ci, "cont_p": float(rc.pvalues[0]),
       "bin_HR": bin_hr, "bin_CI": bin_ci, "bin_p": bin_p, "cutoff": cutoff,
       "pam50": {}, "interaction": None}
for st, s in m.groupby("CLAUDIN_SUBTYPE"):
    try:
        r = PHReg(s.OS_years.values, (s.risk.values / s.risk.std()).reshape(-1, 1),
                  s.OS_event.values.astype(int)).fit(disp=0)
        out["pam50"][str(st)] = {"n": int(len(s)),
                                 "deaths": int(s.OS_event.sum()),
                                 "HR": round(float(np.exp(r.params[0])), 3),
                                 "p": float(r.pvalues[0])}
        log(f"PAM50 {st}: n={len(s)} deaths={int(s.OS_event.sum())} "
            f"HR={float(np.exp(r.params[0])):.3f} p={float(r.pvalues[0]):.3g}")
    except Exception as e:
        out["pam50"][str(st)] = {"n": int(len(s)),
                                 "deaths": int(s.OS_event.sum()),
                                 "note": f"fit failed: {e}"}
        log(f"PAM50 {st}: fit failed ({e}) — N/events only")
D = pd.get_dummies(m.CLAUDIN_SUBTYPE, prefix="pam")
ref = sorted(D.columns)[1]
Xb = np.column_stack([(m.risk.values - m.risk.mean()) / m.risk.std(),
                      D.drop(columns=[ref]).values.astype(float)])
Xi = np.column_stack([(m.risk.values - m.risk.mean()) / m.risk.std() * D[c].values
                      for c in sorted(D.columns) if c != ref])
Xf = np.column_stack([Xb, Xi])
r0 = PHReg(T, Xb, E).fit(disp=0)
r1 = PHReg(T, Xf, E).fit(disp=0)
lrt = float(2 * (r1.llf - r0.llf))
df = Xf.shape[1] - Xb.shape[1]
pint = float(1 - chi2.cdf(lrt, df))
out["interaction"] = {"LR_chi2": lrt, "df": df, "p_nominal": pint, "ref": ref}
log(f"Global interaction: LR={lrt:.3f} df={df} p={pint:.4f}")
def cindex_fast(T, E, r):
    o = np.argsort(T, kind="stable")
    T, E, r = T[o], E[o], r[o]
    ev = np.where(E == 1)[0]
    num = den = 0.0
    for i in ev:
        later = np.arange(i + 1, len(T))
        later = later[T[later] > T[i]]
        if len(later) == 0:
            continue
        den += len(later)
        num += np.sum(r[later] < r[i]) + 0.5 * np.sum(r[later] == r[i])
    return num / den if den else 0.5


json.dump(out, open(RES / "metabric_summary.json", "w"), indent=2)
m[["patientId", "risk"]].to_csv(RES / "metabric_risk.tsv", sep="\t", index=False)
C0 = cindex_fast(T, E, m.risk.values)
Cb = []
for b in range(200):
    idx = np.random.default_rng(SEED + b).integers(0, len(T), len(T))
    Cb.append(cindex_fast(T[idx], E[idx], m.risk.values[idx]))
clo, chi = float(np.quantile(Cb, [0.025, 0.975])[0]), float(np.quantile(Cb, [0.025, 0.975])[1])
out["C"], out["C_CI"] = round(float(C0), 3), [round(clo, 3), round(chi, 3)]
json.dump(out, open(RES / "metabric_summary.json", "w"), indent=2)
log(f"C-index: {C0:.3f} [{clo:.3f},{chi:.3f}] (B=200)")
log("WROTE metabric_summary.json (METABRIC-001)")
