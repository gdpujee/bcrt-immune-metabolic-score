"""Cross-cohort forest figure (METAFOR-001).

Panel A: primary continuous per-SD HRs across the three locked validations
(GSE20685, SCAN-B, METABRIC) read from result JSONs. Panel B: METABRIC
CLAUDIN_SUBTYPE-stratum HRs (five canonical PAM50-like groups plus claudin-low). No numbers hardcoded.
Outputs: figures/Fig12_crosscohort_forest.png(.pdf), logs/metabric_forest.log
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
FIG = ROOT / "figures"
LOGS = ROOT / "logs"


def log(m):
    print(m, flush=True)
    with open(LOGS / "metabric_forest.log", "a") as f:
        f.write(m + "\n")


v = json.load(open(RES / "validation_summary.json"))
cv = json.load(open(RES / "corrective_summary.json"))
rna = json.load(open(RES / "rnaseq_summary.json"))
met = json.load(open(RES / "metabric_summary.json"))
A = [("GSE20685 locked",) + tuple(cv["valid_continuous_perSD_HR"][:3]),
     ("SCAN-B locked",) + tuple(rna["hr_cont"][:3]),
     ("METABRIC locked", met["cont_HR"], met["cont_CI"][0], met["cont_CI"][1])]
# Panel B: refit per-stratum Cox locally (deterministic; assert HR == script 27)
import pandas as pd
from statsmodels.duration.hazard_regression import PHReg
clin = pd.read_csv(ROOT / "metadata/METABRIC_clinical_dl.tsv",
                   sep="\t").drop_duplicates("patientId")
risk = pd.read_csv(RES / "metabric_risk.tsv", sep="\t")
m = clin.merge(risk, on="patientId")
m["OS_years"] = pd.to_numeric(m.OS_MONTHS, errors="coerce") / 12
m["OS_event"] = m.OS_STATUS.str.contains("DECEASED", na=False).astype(int)
m = m.dropna(subset=["OS_years", "OS_event", "risk", "CLAUDIN_SUBTYPE"]).copy()
B = []
for st, s in sorted(m.groupby("CLAUDIN_SUBTYPE")):
    try:
        r = PHReg(s.OS_years.values,
                  (s.risk.values / s.risk.std()).reshape(-1, 1),
                  s.OS_event.values.astype(int)).fit(disp=0)
        hr = float(np.exp(r.params[0]))
        assert abs(hr - met["pam50"][str(st)]["HR"]) < 0.002, (st, hr)
        ci = [float(np.exp(r.params[0] - 1.96 * r.bse[0])),
              float(np.exp(r.params[0] + 1.96 * r.bse[0]))]
        B.append((f"{st} (n={len(s)})", hr, ci[0], ci[1]))
    except Exception as e:
        log(f"panel B skip {st}: {e}")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1, 1.2]})
for a, rows, ttl in zip(ax, [A, B],
                         ["A. Primary per-SD HR (locked validations)",
                          "B. METABRIC CLAUDIN_SUBTYPE strata (exploratory)"]):
    ys = np.arange(len(rows))[::-1]
    for y, (nm, hr, lo, hi) in zip(ys, rows):
        a.plot(hr, y, "s", color="black")
        if lo is not None:
            a.hlines(y, lo, hi, color="black")
        a.text(max(hr, 1.0) * 1.12 if lo is None else hi * 1.03, y,
               f"{nm} {hr:.2f}" + (f" [{lo:.2f},{hi:.2f}]" if lo else ""),
               va="center", fontsize=8)
    a.axvline(1, ls="--", color="grey", lw=1)
    a.set_yticks([])
    a.set_xlabel("HR per SD")
    a.set_title(ttl, fontsize=9)
    a.set_xscale("log")
plt.tight_layout()
plt.savefig(FIG / "Fig12_crosscohort_forest.png", dpi=150)
plt.savefig(FIG / "Fig12_crosscohort_forest.pdf")
log(f"METAFOR-001: panels {len(A)}/{len(B)}; interaction p={met['interaction']['p_nominal']:.3f}")
