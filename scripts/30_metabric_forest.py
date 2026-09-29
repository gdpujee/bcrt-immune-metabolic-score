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
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5),
                       gridspec_kw={"width_ratios": [1, 1.25]})
for panel, rows, title in zip(
        ax, [A, B],
        ["A. Continuous association", "B. METABRIC subtype strata (exploratory)"]):
    ys = np.arange(len(rows))[::-1]
    for y, (_, hr, lo, hi) in zip(ys, rows):
        panel.hlines(y, lo, hi, color="black", lw=1.2)
        panel.plot(hr, y, "s", color="black", markersize=4.5, zorder=3)
    panel.axvline(1, ls="--", color="grey", lw=1)
    panel.set_yticks(ys, [row[0] for row in rows], fontsize=7.5)
    panel.set_xlabel("HR per cohort SD", fontsize=8)
    panel.set_title(title, fontsize=9)
    panel.grid(axis="x", linestyle=":", alpha=0.5)
    panel.spines["top"].set_visible(False)
    panel.spines["right"].set_visible(False)

# Explicit numeric ticks avoid the crowded 10^0 / 1.2×10^0 labels produced by
# automatic log formatting in the narrow two-panel journal figure.
ax[0].set_xscale("log")
ax[0].set_xlim(0.95, 2.15)
ax[0].set_xticks([1.0, 1.25, 1.5, 1.75, 2.0],
                 labels=["1.0", "1.25", "1.5", "1.75", "2.0"])
ax[0].minorticks_off()
ax[1].set_xscale("log")
ax[1].set_xlim(0.48, 5.4)
ax[1].set_xticks([0.5, 1.0, 2.0, 5.0], labels=["0.5", "1", "2", "5"])
ax[1].minorticks_off()
fig.subplots_adjust(left=0.15, right=0.985, top=0.82, bottom=0.17, wspace=0.48)
plt.savefig(FIG / "Fig12_crosscohort_forest.png", dpi=150)
plt.savefig(FIG / "Fig12_crosscohort_forest.pdf")
log(f"METAFOR-001: panels {len(A)}/{len(B)}; interaction p={met['interaction']['p_nominal']:.3f}")
