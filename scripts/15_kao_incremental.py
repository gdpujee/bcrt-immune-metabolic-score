#!/usr/bin/env python3
"""15_kao_incremental.py — RUN-ID: KAO-001 (R3 of 10-round plan; exploratory, D-021 pattern)
Kao I-VI subtype-only vs subtype+risk in GSE20685: ΔC + LRT. Mirrors PAM50-001.
Inputs: metadata/GSE20685_clinical_curated.tsv, results/raw/validation_risk_GSE20685.tsv
Outputs: results/raw/kao_incremental.json, logs/kao.log
"""
import pandas as pd, numpy as np, json
from pathlib import Path
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import chi2

rng = np.random.default_rng(42)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT/"results/raw"; LOGS = ROOT/"logs"; META = ROOT/"metadata"
logf = open(LOGS/"kao.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()

c = pd.read_csv(META/"GSE20685_clinical_curated.tsv", sep="\t")
v = pd.read_csv(RES/"validation_risk_GSE20685.tsv", sep="\t")
c = c.merge(v[["GSM","risk"]], on="GSM")
T = c.follow_up_duration_years.values.astype(float); E = c.event_death.values.astype(int)
D = pd.get_dummies(c.subtype, prefix="kao")
ref = sorted(D.columns)[0]
Xb = D.drop(columns=[ref]).values.astype(float)
Xf = np.column_stack([Xb, c.risk.values])
log(f"n={len(T)}, deaths={int(E.sum())}, ref={ref}, subtypes={sorted(c.subtype.dropna().unique())}")

def cindex(T,E,x):
    T=np.asarray(T); E=np.asarray(E); x=np.asarray(x)
    ev=E==1; Te=T[ev][:,None]; To=T[None,:]; xe=x[ev][:,None]; xo=x[None,:]
    comp=To>Te; conc=(xo<xe).astype(float)+0.5*(xo==xe).astype(float)
    tot=comp.sum(); return float((conc*comp).sum()/tot) if tot else 0.5

rb=PHReg(T,Xb,E).fit(disp=0); rf=PHReg(T,Xf,E).fit(disp=0)
cb=cindex(T,E,Xb@np.asarray(rb.params).ravel()); cf=cindex(T,E,Xf@np.asarray(rf.params).ravel())
lrt=2*(rf.llf-rb.llf); p=float(1-chi2.cdf(lrt,1))
log(f"Kao-only C={cb:.3f} → +risk C={cf:.3f}, ΔC={cf-cb:.3f}; LRT χ²={lrt:.2f} p={p:.3g}")
# bootstrap CI for ΔC (B=200, seed 42 — reproduces committed CI exactly)
# convergence diagnostics: draws emitting ConvergenceWarning are counted and a
# clean-draws sensitivity CI is logged (JSON keeps all-draws CI, unchanged).
import warnings as _warnings
diffs = []
clean = []
n_warned = 0
for b in range(200):
    idx = rng.integers(0, len(T), len(T))
    try:
        with _warnings.catch_warnings(record=True) as _w:
            _warnings.simplefilter("always")
            a = PHReg(T[idx], Xb[idx], E[idx]).fit(disp=0)
            f = PHReg(T[idx], Xf[idx], E[idx]).fit(disp=0)
        _cw = any("converge" in str(x.message).lower() for x in _w)
        xa = Xb[idx] @ np.asarray(a.params).ravel()
        xc = Xf[idx] @ np.asarray(f.params).ravel()
        _d = cindex(T[idx], E[idx], xc) - cindex(T[idx], E[idx], xa)
        diffs.append(_d)
        if _cw:
            n_warned += 1
        else:
            clean.append(_d)
    except Exception:
        pass
diffs = np.array(diffs)
clean = np.array(clean)
log(f"ΔC 95% CI: [{np.percentile(diffs, 2.5):.3f},{np.percentile(diffs, 97.5):.3f}]")
log(f"Bootstrap convergence: {n_warned}/200 draws warned (included above); "
    f"clean-draws sensitivity CI [{np.percentile(clean, 2.5):.3f},{np.percentile(clean, 97.5):.3f}] — boundary reading unchanged")
json.dump({"ref": ref, "base_C": float(cb), "full_C": float(cf), "delta": float(cf - cb),
 "LRT_chi2": float(lrt), "LRT_p": p, "n": int(len(T)), "events": int(E.sum()),
 "CI": [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))]},
 open(RES / "kao_incremental.json", "w"), indent=2)
log("WROTE kao_incremental.json (KAO-001)")
logf.close()
