#!/usr/bin/env python3
"""06_validation.py — RUN-ID: VALID-001
Locked application of TRAIN-001 model to GSE20685 (n=327) + biology checks in GSE45827. NO retraining, NO cutoff tuning.
Inputs: results/derived/locked_model.json, data/processed/*, metadata/*_curated
Outputs: results/raw/validation_*.json/tsv, figures/valid_KM.png, tables/*
"""
import pandas as pd, numpy as np, json
from pathlib import Path
from statsmodels.duration.hazard_regression import PHReg
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT/"data/processed"; META = ROOT/"metadata"; RES = ROOT/"results/raw"; DER = ROOT/"results/derived"; FIG = ROOT/"figures"; TAB = ROOT/"tables"
for d in [RES, FIG, TAB]: d.mkdir(parents=True, exist_ok=True)
logf = open(ROOT/"logs/validation.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()

locked = json.load(open(DER/"locked_model.json"))
genes = locked["genes"]; coefs = np.array(locked["coefs"]); means = locked["scaling"]["means"]; sds = locked["scaling"]["sds"]; cutoff = locked["cutoff"]
log(f"Locked model v={locked['version']}: {len(genes)} genes, cutoff={cutoff:.4f}")

def cindex(T,E,risk):
    n=len(T); conc=tot=0
    for i in range(n):
        if E[i]!=1: continue
        for j in range(n):
            if T[j]>T[i]:
                tot+=1
                if risk[j]<risk[i]: conc+=1
                elif risk[j]==risk[i]: conc+=0.5
    return conc/tot if tot else 0.5

def km(T,E):
    order=np.argsort(T); T=T[order]; E=E[order]
    uniq=np.sort(np.unique(T[E==1])); surv=[]; s=1.0
    for t in uniq:
        at=(T>=t).sum(); d=int(((T==t)&(E==1)).sum())
        s=s*(1-d/at) if at>0 else s; surv.append(s)
    return uniq, np.array(surv)

def logrank(Ta,Ea,Tb,Eb):
    from scipy.stats import chi2
    times=np.sort(np.unique(np.concatenate([Ta[Ea==1],Tb[Eb==1]])))
    O1=E1=V=0.0
    for t in times:
        n1=(Ta>=t).sum(); n2=(Tb>=t).sum(); n=n1+n2
        d1=((Ta==t)&(Ea==1)).sum(); d2=((Tb==t)&(Eb==1)).sum(); d=d1+d2
        if n>1 and d>0:
            E1+=d*n1/n; O1+=d1; V+=(n1*n2*d*(n-d))/(n*n*(n-1))
    chi=(O1-E1)**2/V if V>0 else 0
    return float(chi), float(1-chi2.cdf(chi,1))

# --- GSE20685 locked application ---
c206 = pd.read_csv(META/"GSE20685_clinical_curated.tsv", sep="\t")
e206 = pd.read_csv(PROC/"GSE20685_expr_gene.tsv", sep="\t", index_col=0)
missing = [g for g in genes if g not in e206.index]
log(f"Missing genes in validation: {missing} (must be [] for locked validation)")
assert len(missing)==0, f"Missing {missing} — cannot claim external validation"
T = c206["OS_years"].values.astype(float); E = c206["OS_event"].values.astype(int)
X = e206.loc[genes, c206["GSM"]].T
# apply TRAINING scaling (frozen)
Z = (X - pd.Series(means)) / pd.Series(sds)
risk = Z.values @ coefs
group = np.where(risk > cutoff, "High", "Low")
log(f"Validation risk: High={(group=='High').sum()}, Low={(group=='Low').sum()}")
hi=(group=="High"); lo=(group=="Low")
th,sh = km(T[hi],E[hi]); tl,sl = km(T[lo],E[lo])
chi2v, plog = logrank(T[hi],E[hi],T[lo],E[lo])
res_g = PHReg(T, hi.astype(int), E).fit(disp=0)
HR = float(np.exp(res_g.params[0])); p_hr = float(res_g.pvalues[0])
# multivariable adjusted (age + T/N/M + subtype as proxies; subtype categorical -> dummies)
c206["age"] = pd.to_numeric(c206["age_at_diagnosis"], errors="coerce")
for c in ["t_stage","n_stage","m_stage"]:
    c206[c] = pd.to_numeric(c206[c], errors="coerce")
adj = pd.DataFrame({"risk": risk, "age": c206["age"], "t": c206["t_stage"], "n": c206["n_stage"]})
adj = adj.fillna(adj.median())
res_adj = PHReg(T, adj.values, E).fit(disp=0)
log(f"Adjusted Cox (risk+age+T+N): risk HR={np.exp(res_adj.params[0]):.2f} p={res_adj.pvalues[0]:.3g}")
C = cindex(T,E,risk)
# Descriptive known-status cumulative/dynamic classification AUC: events by t0
# are cases, follow-up beyond t0 are controls, and censored observations at or
# before t0 without an event are excluded. This is not IPCW-corrected.
def td_auc(t0):
    y=((T<=t0)&(E==1)).astype(int); mask=((T<=t0)&(E==1))|(T>t0)
    if y[mask].sum()==0 or y[mask].sum()==mask.sum(): return float("nan"), int(mask.sum())
    return float(roc_auc_score(y[mask], risk[mask])), int(mask.sum())
aucs={f"{t}y": td_auc(t) for t in [1,3,5]}
log(f"VALIDATION: log-rank p={plog:.4g}, HR={HR:.2f} p={p_hr:.4g}, C={C:.3f}, AUCs={aucs}")
# save
pd.DataFrame({"GSM": c206["GSM"], "risk": risk, "group": group, "OS_years": T, "OS_event": E}).to_csv(RES/"validation_risk_GSE20685.tsv", sep="\t", index=False)
with open(RES/"validation_summary.json","w") as f:
    json.dump({"n":int(len(T)),"events":int(E.sum()),"HR":HR,"p_hr":p_hr,"logrank_p":plog,"C":C,
     "aucs":{k:[float(v[0]) if v[0]==v[0] else None,v[1]] for k,v in aucs.items()},
     "adj_risk_HR":float(np.exp(res_adj.params[0])),"adj_risk_p":float(res_adj.pvalues[0])},f,indent=2)
plt.figure(figsize=(6,5))
plt.step(np.concatenate([[0],th]), np.concatenate([[1],sh]), where="post", label=f"High (n={hi.sum()})")
plt.step(np.concatenate([[0],tl]), np.concatenate([[1],sl]), where="post", label=f"Low (n={lo.sum()})")
plt.ylim(0,1.02); plt.xlabel("Years"); plt.ylabel("OS probability")
plt.legend(); plt.tight_layout(); plt.savefig(FIG/"valid_KM.png", dpi=150); plt.savefig(FIG/"valid_KM.pdf")

# --- GSE45827 biology: risk distribution by subtype + tumor vs normal ---
c458 = pd.read_csv(META/"GSE45827_clinical_curated.tsv", sep="\t")
e458 = pd.read_csv(PROC/"GSE45827_expr_gene.tsv", sep="\t", index_col=0)
avail = [g for g in genes if g in e458.index]
log(f"45827 genes available {len(avail)}/{len(genes)} (pre-filtered matrix)")
X458 = e458.loc[avail, c458["GSM"]].T
Z458 = (X458 - pd.Series({g:means[g] for g in avail})) / pd.Series({g:sds[g] for g in avail})
risk458 = Z458.values @ np.array([coefs[genes.index(g)] for g in avail])
c458["risk_14g"] = risk458  # full 14-gene locked score applied descriptively (no survival endpoint — exploratory, NOT locked validation)
# tumor vs normal + subtype box stats
bio = c458[~c458["is_cell_line"]].copy()
from scipy.stats import mannwhitneyu
t_r = bio[bio.diagnosis=="Breast cancer"]["risk_14g"].values
n_r = bio[bio.diagnosis.str.contains("None", na=False)]["risk_14g"].values
_, p_tn = mannwhitneyu(t_r, n_r, alternative="two-sided")
log(f"45827 tumor(n={len(t_r)}) vs normal(n={len(n_r)}): median {np.median(t_r):.3f} vs {np.median(n_r):.3f}, MW p={p_tn:.3g} (EXPLORATORY full 14-gene score, no survival endpoint)")
bio[["GSM","diagnosis","tumor_subtype","batch","risk_14g"]].to_csv(RES/"biology_risk_GSE45827.tsv", sep="\t", index=False)
# subtype medians
log("Subtype medians:\n"+str(bio.groupby("tumor_subtype")["risk_14g"].median()))
log("WROTE validation + biology outputs")
logf.close()
