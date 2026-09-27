"""Schoenfeld residual PH diagnostics (PHZPH-001).

Approximate Schoenfeld rank-correlation screen for the locked-score Cox models:
per-variable z = corr * sqrt(events-1), with a correlation-based global check.
Exact Grambsch-Therneau tests are run independently by 32_ph_exact.R.
Models: GSE20685/SCAN-B continuous and covariate-adjusted, METABRIC continuous.
Figure: scaled residuals vs time for all three continuous scores.
Outputs: results/raw/ph_diagnostics.json, figures/PH_schoenfeld.png(.pdf),
logs/ph_diagnostics.log
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import chi2, norm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
DER = ROOT / "results/derived"
META = ROOT / "metadata"
FIG = ROOT / "figures"
LOGS = ROOT / "logs"
SEED = 42
np.random.seed(SEED)


def log(m):
    print(m, flush=True)
    with open(LOGS / "ph_diagnostics.log", "a") as f:
        f.write(m + "\n")


def schoenfeld(T, E, X, beta):
    """Return (raw residuals d x p, aligned ranked event times, d)."""
    X = np.asarray(X, float)
    xb = X @ np.asarray(beta, float).ravel()
    w = np.exp(xb - xb.max())
    ev = np.where(E == 1)[0]
    Tres = []
    for t in np.unique(T[ev]):
        idx = np.where((T == t) & (E == 1))[0]
        R = np.where(T >= t)[0]
        wr = w[R]
        Ex = (wr[:, None] * X[R]).sum(0) / wr.sum()
        for i in idx:
            Tres.append(X[i] - Ex)
    R = np.array(Tres)  # raw Schoenfeld residuals; GT scaling applied by caller
    # ranks MUST follow the same time-sorted order as R (ties averaged);
    # pairing time-ordered residuals with original-order ranks is invalid.
    from scipy.stats import rankdata
    g = rankdata(np.sort(T[ev])).astype(float)
    assert len(g) == len(R)
    return R, g, len(ev)


def zph(T, E, X, names):
    fit = PHReg(np.asarray(T, float), np.asarray(X, float),
                np.asarray(E, int)).fit(disp=0)
    S_raw, g, d = schoenfeld(T, E, X, fit.params)
    # Grambsch-Therneau scaling: S = R @ Cov(beta-hat)
    S = S_raw @ np.asarray(fit.cov_params())
    gc = g - g.mean()
    out = {"n": int(len(T)), "events": int(d), "vars": {}}
    zs = []
    for k, nm in enumerate(names):
        s = S[:, k]
        r = float(np.corrcoef(s, g)[0, 1]) if np.std(s) > 0 else 0.0
        z = r * np.sqrt(max(d - 1, 1))
        p = float(2 * norm.sf(abs(z)))
        p2 = float(chi2.sf(z * z, 1))
        assert abs(p - p2) < 1e-9, "normal/chi2 p mismatch"
        out["vars"][nm] = {"z": round(z, 3), "p": p}
        zs.append(z)
        log(f"{nm}: z={z:.3f} p={p:.4g} (d={d})")
    if len(names) > 1:
        R = np.corrcoef(S.T)
        try:
            gchi = float(np.array(zs) @ np.linalg.inv(R + 1e-9 * np.eye(len(zs))) @ np.array(zs))
            gp = float(chi2.sf(gchi, len(zs)))
        except Exception:
            gchi, gp = float("nan"), float("nan")
        out["global"] = {"chi2": round(gchi, 3), "df": len(names), "p": gp}
        log(f"GLOBAL: chi2={gchi:.3f} df={len(names)} p={gp:.4g}")
    return out, fit, (S, g)


res = {}
# GSE20685 continuous
c = pd.read_csv(META / "GSE20685_clinical_curated.tsv", sep="\t")
v = pd.read_csv(RES / "validation_risk_GSE20685.tsv", sep="\t")
c = c.merge(v[["GSM", "risk"]], on="GSM")
T = c.follow_up_duration_years.values.astype(float)
E = c.event_death.values.astype(int)
sd = c.risk.std()
log("== GSE20685 continuous ==")
o, _, sg = zph(T, E, (c.risk.values / sd).reshape(-1, 1), ["risk"])
res["GSE20685_cont"] = o
# GSE20685 adjusted
for cc in ["age_at_diagnosis", "t_stage", "n_stage"]:
    c[cc] = pd.to_numeric(c[cc], errors="coerce")
Xa = c[["risk", "age_at_diagnosis", "t_stage", "n_stage"]].fillna(
    c[["risk", "age_at_diagnosis", "t_stage", "n_stage"]].median()).values
log("== GSE20685 adjusted ==")
o, _, _ = zph(T, E, Xa, ["risk", "age", "T", "N"])
res["GSE20685_adj"] = o
# GSE20685 5y split (administrative censoring at 5y on full cohort)
e1T = np.minimum(T, 5)
e1E = ((T <= 5) & (E == 1)).astype(int)
e1r = (c.risk.values / sd).reshape(-1, 1)
rr = PHReg(e1T, e1r, e1E).fit(disp=0)
res["GSE20685_split5y"] = {
    "early": {"n": int(len(e1T)), "deaths": int(e1E.sum()),
              "HR": round(float(np.exp(rr.params[0])), 3),
              "CI": [round(float(np.exp(rr.params[0] - 1.96 * rr.bse[0])), 4),
                     round(float(np.exp(rr.params[0] + 1.96 * rr.bse[0])), 4)],
              "p": float(rr.pvalues[0])}}
late_mask = T > 5
rr2 = PHReg(T[late_mask] - 5, ((c.risk.values / sd)[late_mask]).reshape(-1, 1),
            E[late_mask]).fit(disp=0)
res["GSE20685_split5y"]["late"] = {
    "n": int(late_mask.sum()), "deaths": int(E[late_mask].sum()),
    "HR": round(float(np.exp(rr2.params[0])), 3),
    "CI": [round(float(np.exp(rr2.params[0] - 1.96 * rr2.bse[0])), 4),
           round(float(np.exp(rr2.params[0] + 1.96 * rr2.bse[0])), 4)],
    "p": float(rr2.pvalues[0])}
log(f"GSE20685 split5y early HR={res['GSE20685_split5y']['early']['HR']} "
    f"late HR={res['GSE20685_split5y']['late']['HR']}")
# SCAN-B continuous
r = pd.read_csv(RES / "rnaseq_risk_GSE96058.tsv", sep="\t")
T2, E2 = r.OS_years.values, r.OS_event.values.astype(int)
sd2 = r.risk.std()
log("== SCAN-B continuous ==")
o, _, sg2 = zph(T2, E2, (r.risk.values / sd2).reshape(-1, 1), ["risk"])
res["SCANB_cont"] = o
# SCAN-B adjusted
cu = pd.read_csv(META / "GSE96058_clinical_raw.tsv", sep="\t")
cu = cu[~cu.title.str.contains("repl", na=False)].copy()
m = cu.merge(r[["sample", "risk"]], left_on="title", right_on="sample")
for cc in ["age_at_diagnosis", "er_status", "her2_status"]:
    m[cc] = pd.to_numeric(m[cc], errors="coerce")
Xb = m[["risk", "age_at_diagnosis", "er_status", "her2_status"]].fillna(
    m[["risk", "age_at_diagnosis", "er_status", "her2_status"]].median()).values
T3 = (m.overall_survival_days.values / 365.25)
E3 = m.overall_survival_event.values.astype(int)
log("== SCAN-B adjusted ==")
o, _, _ = zph(T3, E3, Xb, ["risk", "age", "ER", "HER2"])
res["SCANB_adj"] = o
# METABRIC continuous (PHZPH-002: time-varying hazard assessment)
mc = pd.read_csv(META / "METABRIC_clinical_dl.tsv", sep="\t").drop_duplicates("patientId")
mr = pd.read_csv(RES / "metabric_risk.tsv", sep="\t")
mm = mc.merge(mr, on="patientId")
mm["Ty"] = pd.to_numeric(mm.OS_MONTHS, errors="coerce") / 12
mm["Ee"] = mm.OS_STATUS.str.contains("DECEASED", na=False).astype(int)
mm = mm.dropna(subset=["Ty", "Ee", "risk"]).copy()
mm = mm[np.isfinite(mm["Ty"].values)]
T4, E4 = mm["Ty"].values, mm["Ee"].values.astype(int)
sd4 = float(mm.risk.std())
log("== METABRIC continuous ==")
o, _, sg4 = zph(T4, E4, (mm.risk.values / sd4).reshape(-1, 1), ["risk"])
res["METABRIC_cont"] = o
splits = {}
for cut in [3, 5, 7]:
    # early: administrative censoring at cut on the FULL cohort
    # (time=min(T,cut), event=E&(T<=cut)); subsetting to T<=cut drops
    # at-risk person-time and biases the estimate.
    e1 = mm[["Ty", "Ee", "risk"]].copy()
    e1["Ty"] = np.minimum(e1["Ty"].values, cut)
    e1["Ee"] = ((mm["Ty"].values <= cut) & (mm["Ee"].values == 1)).astype(int)
    # late: conditional on survival past cut (delayed entry at cut)
    e2 = mm[mm.Ty > cut].copy()
    e2["Ty"] = e2["Ty"] - cut
    row = {}
    for tag, d in [("early", e1), ("late", e2)]:
        rr = PHReg(d.Ty.values, (d.risk.values / sd4).reshape(-1, 1),
                   d.Ee.values.astype(int)).fit(disp=0)
        row[tag] = {"n": int(len(d)), "deaths": int(d.Ee.sum()),
                    "HR": round(float(np.exp(rr.params[0])), 3),
                    "CI": [round(float(np.exp(rr.params[0] - 1.96 * rr.bse[0])), 4),
                           round(float(np.exp(rr.params[0] + 1.96 * rr.bse[0])), 4)],
                    "p": float(rr.pvalues[0])}
        log(f"METABRIC split{cut}y {tag}: n={len(d)} d={int(d.Ee.sum())} "
            f"HR={float(np.exp(rr.params[0])):.3f} p={float(rr.pvalues[0]):.3g}")
    splits[str(cut)] = row
res["METABRIC_timesplit"] = splits
json.dump(res, open(RES / "ph_diagnostics.json", "w"), indent=2)

# figure: scaled residuals vs ranked time
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for a, (S, g), ttl in zip(ax, [sg, sg2, sg4], ["GSE20685 continuous", "SCAN-B continuous", "METABRIC continuous"]):
    a.scatter(g, S[:, 0], s=4, alpha=0.4)
    q = np.quantile(g, np.linspace(0, 1, 11))
    bx = [(q[k] + q[k + 1]) / 2 for k in range(10)]
    by = [S[:, 0][(g >= q[k]) & (g <= q[k + 1])].mean() for k in range(10)]
    a.plot(bx, by, color="red", lw=2)
    a.axhline(0, ls="--", lw=1)
    a.set_xlabel("ranked event time")
    a.set_ylabel("scaled Schoenfeld residual (risk)")
    a.set_title(ttl, fontsize=9)
plt.tight_layout()
plt.savefig(FIG / "PH_schoenfeld.png", dpi=150)
plt.savefig(FIG / "PH_schoenfeld.pdf")
log("WROTE ph_diagnostics.json + PH_schoenfeld (PHZPH-001)")
