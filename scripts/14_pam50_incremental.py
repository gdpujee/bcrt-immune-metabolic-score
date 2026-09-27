#!/usr/bin/env python3
"""14_pam50_incremental.py — RUN-ID: PAM50-001 (amendment D-021, exploratory)
PAM50-subtype-only vs PAM50+risk: ΔC with bootstrap CI + likelihood-ratio test, SCAN-B.
PAM50 as dummy variables (LumA reference). No refitting of locked score.
Inputs: results/raw/rnaseq_risk_GSE96058.tsv
Outputs: results/raw/pam50_incremental.json, logs/pam50.log
Seed 42.
"""
import pandas as pd, numpy as np, json
from pathlib import Path
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT/"results/raw"; LOGS = ROOT/"logs"
logf = open(LOGS/"pam50.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()
rng = np.random.default_rng(42)

r = pd.read_csv(RES/"rnaseq_risk_GSE96058.tsv", sep="\t")
T = r.OS_years.values.astype(float); E = r.OS_event.values.astype(int)
D = pd.get_dummies(r.pam50_subtype, prefix="pam", drop_first=False)
# reference LumA: drop it
Xb = D.drop(columns=["pam_LumA"]).values.astype(float)
Xf = np.column_stack([Xb, r.risk.values])
log(f"n={len(T)}, deaths={int(E.sum())}, base cols={list(D.drop(columns=['pam_LumA']).columns)}")

def cindex(T,E,x):
    # vectorized Harrell's C (ties 0.5)
    T=np.asarray(T); E=np.asarray(E); x=np.asarray(x)
    ev = E==1
    Te=T[ev][:,None]; To=T[None,:]; xe=x[ev][:,None]; xo=x[None,:]
    comp = To > Te
    conc = (xo < xe).astype(float) + 0.5*(xo == xe).astype(float)
    tot = comp.sum()
    return float((conc*comp).sum()/tot) if tot else 0.5

rb = PHReg(T,Xb,E).fit(disp=0); rf = PHReg(T,Xf,E).fit(disp=0)
cb = cindex(T,E,Xb@np.asarray(rb.params).ravel()); cf = cindex(T,E,Xf@np.asarray(rf.params).ravel())
# LRT: full vs reduced (1 df)
lrt = 2*(rf.llf - rb.llf); p_lrt = float(1-chi2.cdf(lrt,1))
log(f"PAM50-only C={cb:.3f} → +risk C={cf:.3f}, ΔC={cf-cb:.3f}; LRT χ²={lrt:.2f} p={p_lrt:.3g}")
# bootstrap CI for ΔC
diffs=[]
for b in range(500):
    idx = rng.integers(0,len(T),len(T))
    try:
        a=PHReg(T[idx],Xb[idx],E[idx]).fit(disp=0); c=PHReg(T[idx],Xf[idx],E[idx]).fit(disp=0)
        xa=Xb[idx]@np.asarray(a.params).ravel(); xc=Xf[idx]@np.asarray(c.params).ravel()
        diffs.append(cindex(T[idx],E[idx],xc)-cindex(T[idx],E[idx],xa))
    except Exception: pass
diffs=np.array(diffs)
log(f"ΔC 95% CI: [{np.percentile(diffs,2.5):.3f},{np.percentile(diffs,97.5):.3f}]")
json.dump({"base_C":float(cb),"full_C":float(cf),"delta":float(cf-cb),
 "CI":[float(np.percentile(diffs,2.5)),float(np.percentile(diffs,97.5))],
 "LRT_chi2":float(lrt),"LRT_p":p_lrt,"n":int(len(T)),"events":int(E.sum())},
 open(RES/"pam50_incremental.json","w"),indent=2)
log("WROTE pam50_incremental.json (PAM50-001)")
logf.close()
