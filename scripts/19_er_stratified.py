"""ER-stratified sensitivity in SCAN-B (ER-001, exploratory).

Per-SD Cox within ER strata from the locked risk table. ER status known for
3073/3273 patients (200 missing, excluded with disclosure — not silently).
Contract: output byte-identical to committed results/raw/er_stratified.json,
else exit nonzero. Archiving rule: HR round(x,3), p full precision.
"""
import json
import subprocess
import pandas as pd
import numpy as np
from statsmodels.duration.hazard_regression import PHReg
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
SEED = 42
np.random.seed(SEED)


def log(m):
    print(m, flush=True)
    with open(LOGS / "er_stratified.log", "a") as f:
        f.write(m + "\n")


rs = pd.read_csv(RES / "rnaseq_risk_GSE96058.tsv", sep="\t")
cu = pd.read_csv(META / "GSE96058_clinical_raw.tsv", sep="\t")
cu = cu[~cu.title.str.contains("repl", na=False)].copy()
cu["er"] = pd.to_numeric(cu.er_status, errors="coerce")
m = cu.merge(rs[["sample", "risk"]], left_on="title", right_on="sample")
known = m[m.er.isin([0, 1])].copy()
log(f"ER known {len(known)}/3273; missing {int(m.er.isna().sum())} excluded (disclosed)")
out = {}
for er, key in [(0, "ER0"), (1, "ER1")]:
    s = known[known.er == er]
    sd = s.risk.std()
    r = PHReg(s.overall_survival_days.values / 365.25,
              (s.risk.values / sd).reshape(-1, 1),
              s.overall_survival_event.values).fit(disp=0)
    hr = float(np.exp(r.params[0]))
    p = float(r.pvalues[0])
    out[key] = {"n": int(len(s)), "deaths": int(s.overall_survival_event.sum()),
                "HR": round(hr, 3), "p": p}
    log(f"{key}: n={len(s)} deaths={int(s.overall_survival_event.sum())} "
        f"HR={round(hr,3)} p={p:.4g}")
text = json.dumps(out, indent=2)
ref = subprocess.run(["git", "show", f"HEAD:results/raw/er_stratified.json"],
                     capture_output=True, text=True, cwd=ROOT).stdout
if text != ref:
    open("/tmp/er_new.json", "w").write(text)
    raise SystemExit("MISMATCH vs committed er_stratified.json (see /tmp/er_new.json)")
open(RES / "er_stratified.json", "w").write(text)
log("ER-001: er_stratified.json exported + verified vs HEAD; PASS")
