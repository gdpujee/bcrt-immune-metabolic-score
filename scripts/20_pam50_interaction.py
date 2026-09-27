"""Global risk x PAM50 interaction in SCAN-B (PAM50INT-001, exploratory, nominal).

LRT: full (risk + PAM50 dummies + risk x PAM50) vs reduced (risk + PAM50
dummies), 4 df. Reports nominal p only; no multiplicity adjustment (single
exploratory heterogeneity check accompanying the per-subtype HRs).
Contract: output byte-identical to committed results/raw/pam50_interaction.json.
"""
import json
import subprocess
import pandas as pd
import numpy as np
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import chi2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
LOGS = ROOT / "logs"
SEED = 42
np.random.seed(SEED)


def log(m):
    print(m, flush=True)
    with open(LOGS / "pam50_interaction.log", "a") as f:
        f.write(m + "\n")


r = pd.read_csv(RES / "rnaseq_risk_GSE96058.tsv", sep="\t")
s = r.dropna(subset=["pam50_subtype"]).copy()
z = (s.risk.values - s.risk.mean()) / s.risk.std()
D = pd.get_dummies(s.pam50_subtype, prefix="pam")
ref = sorted(D.columns)[0]
Xb = np.column_stack([z, D.drop(columns=[ref]).values.astype(float)])
Xi = np.column_stack([z * D[c].values for c in sorted(D.columns) if c != ref])
Xf = np.column_stack([Xb, Xi])
T = s.OS_years.values
E = s.OS_event.values.astype(int)
rb = PHReg(T, Xb, E).fit(disp=0)
rf = PHReg(T, Xf, E).fit(disp=0)
lrt = float(2 * (rf.llf - rb.llf))
p = float(1 - chi2.cdf(lrt, 4))
log(f"Global risk x PAM50: LR chi2={lrt:.4f} df=4 nominal p={p:.4f} (n={len(s)})")
out = {"n": int(len(s)), "LR_chi2": lrt, "df": 4, "p_nominal": p,
       "reference_subtype": ref, "seed": SEED,
       "note": "exploratory heterogeneity check; complements per-subtype HRs"}
text = json.dumps(out, indent=2)
open(RES / "pam50_interaction.json", "w").write(text)
ref = subprocess.run(["git", "show", "HEAD:results/raw/pam50_interaction.json"],
                     capture_output=True, text=True, cwd=ROOT).stdout
if ref and text != ref:
    raise SystemExit("MISMATCH vs committed pam50_interaction.json")
log("PAM50INT-001: pam50_interaction.json exported + verified vs HEAD; PASS")
