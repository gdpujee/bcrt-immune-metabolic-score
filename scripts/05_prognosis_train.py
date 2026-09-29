#!/usr/bin/env python3
"""05_prognosis_train.py — RUN-ID: TRAIN-001
Training ONLY on GSE42568 tumors (n=104, 35 events). Univariable Cox screen + LASSO-Cox (5-fold CV, seed 42) + multivariable refit + locked RiskScore.
Outputs: results/raw/train_*.json/tsv, results/derived/locked_model.json, figures/train_KM.png, logs/train.log
Leakage: all selection/scaling/cutoff fit on training only. Validation untouched.
"""
import pandas as pd, numpy as np, json
from pathlib import Path
from statsmodels.duration.hazard_regression import PHReg
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT/"data/processed"; META = ROOT/"metadata"; RES = ROOT/"results/raw"; DER = ROOT/"results/derived"; FIG = ROOT/"figures"; LOGS = ROOT/"logs"
for d in [RES, DER, FIG]: d.mkdir(parents=True, exist_ok=True)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SEED = 42; np.random.seed(SEED)
logf = open(LOGS/"train.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()

# load
clin = pd.read_csv(META/"GSE42568_clinical_curated.tsv", sep="\t")
expr = pd.read_csv(PROC/"GSE42568_expr_gene.tsv", sep="\t", index_col=0)
pool = pd.read_csv(RES/"candidate_pool.tsv", sep="\t")["symbol"].tolist()
train_gsm = clin[(clin.tissue=="breast cancer") & clin.OS_time_days.notna() & clin.OS_event.notna()]["GSM"].tolist()
log(f"Training tumors: {len(train_gsm)} (events={(clin.set_index('GSM').loc[train_gsm,'OS_event']==1).sum()})")
X = expr[train_gsm].T  # samples x genes
X = X[[g for g in pool if g in X.columns]]
log(f"Candidate pool available in training: {X.shape[1]} / {len(pool)}")
T = clin.set_index("GSM").loc[train_gsm, "OS_time_days"].values.astype(float)/365.25
E = clin.set_index("GSM").loc[train_gsm, "OS_event"].values.astype(int)
# z-score per gene using training stats (freeze)
means = X.mean(axis=0); sds = X.std(axis=0, ddof=0).replace(0, 1.0)
Z = (X - means)/sds

# 1) Full-derivation-cohort univariable Cox screen p<0.01. This screen is run once
# before the folds below; the 5-fold CV tunes alpha within the screened set and is
# not a nested estimate of full-pipeline performance.
rows=[]
for g in Z.columns:
    try:
        res = PHReg(T, Z[[g]].values, E).fit(disp=0)
        p = float(res.pvalues[0]); hr = float(np.exp(res.params[0]))
    except Exception as e:
        p, hr = 1.0, 1.0
    rows.append((g,p,hr))
uni = pd.DataFrame(rows, columns=["symbol","p_uni","HR_uni"]).sort_values("p_uni")
uni.to_csv(RES/"train_univariable_cox.tsv", sep="\t", index=False)
sig = uni[uni.p_uni<0.01]
log(f"Univariable p<0.01: {len(sig)} / {len(uni)}")
if len(sig)==0:
    raise SystemExit("FAILED: no genes pass univariable screen — STOP per skill (negative result preserved, no threshold shopping)")
Zs = Z[sig.symbol.tolist()]

# 2) LASSO-Cox via elastic_net L1_wt=1, 5-fold CV over alpha grid, maximize C-index
def cindex(T, E, risk):
    # Harrell's C: comparable pairs where earlier event has higher risk
    n=len(T); conc=0; tot=0
    for i in range(n):
        if E[i]!=1: continue
        for j in range(n):
            if T[j] > T[i]:
                tot+=1
                if risk[j] < risk[i]: conc+=1
                elif risk[j]==risk[i]: conc+=0.5
    return conc/tot if tot else 0.5

alphas = np.logspace(-3, 1, 20)  # 0.001..10
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
cv_scores = []
for a in alphas:
    cs=[]
    for tr, va in skf.split(Zs.values, E):
        try:
            fit = PHReg(T[tr], Zs.values[tr], E[tr]).fit_regularized(method="elastic_net", alpha=a, L1_wt=1.0, maxiter=200)
            pvec = np.asarray(fit.params).ravel()
            nz = np.abs(pvec) > 1e-6
            if nz.sum()==0:
                cs.append(0.5); continue
            risk_va = Zs.values[va][:, nz] @ pvec[nz]
            cs.append(cindex(T[va], E[va], np.asarray(risk_va).ravel()))
        except Exception as e:
            cs.append(0.5)
    cv_scores.append(float(np.mean(cs)))
    log(f"alpha={a:.4f} CV C={np.mean(cs):.3f}")
best_alpha = alphas[int(np.argmax(cv_scores))]
log(f"Best alpha={best_alpha:.4f} (CV C={max(cv_scores):.3f})")
# refit on full training with best alpha
final_reg = PHReg(T, Zs.values, E).fit_regularized(method="elastic_net", alpha=best_alpha, L1_wt=1.0, maxiter=500)
params = np.asarray(final_reg.params).ravel()
sel = sig.symbol.tolist()
selected = [(g, float(p)) for g,p in zip(sel, params) if abs(p) > 1e-6]
log(f"LASSO selected {len(selected)} / {len(sel)}: {selected}")
if len(selected)==0:
    raise SystemExit("FAILED: LASSO selected 0 genes — negative result, stop (no shopping)")
sel_genes = [g for g,_ in selected]

# 3) multivariable refit (unpenalized) on selected
Xm = Zs[sel_genes].values
res_mv = PHReg(T, Xm, E).fit(disp=0)
coefs = np.asarray(res_mv.params).ravel()
ses = np.asarray(res_mv.bse).ravel()
hrs = np.exp(coefs); lo = np.exp(coefs-1.96*ses); hi = np.exp(coefs+1.96*ses)
mv = pd.DataFrame({"symbol": sel_genes, "coef": coefs, "SE": ses, "HR": hrs, "HR_lo95": lo, "HR_hi95": hi, "p": np.asarray(res_mv.pvalues).ravel()})
mv.to_csv(RES/"train_multivariable_cox.tsv", sep="\t", index=False)
log("Multivariable:\n"+mv.to_string(index=False))

# risk score + median cutoff (locked)
risk = Xm @ coefs
cutoff = float(np.median(risk))
group = np.where(risk > cutoff, "High", "Low")
log(f"Risk median cutoff (locked) = {cutoff:.4f}; High={(group=='High').sum()}, Low={(group=='Low').sum()}")

# KM + log-rank (hand-rolled)
def km(T, E, times):
    order = np.argsort(T); T=T[order]; E=E[order]
    uniq = np.sort(np.unique(T[E==1]))
    surv=[]; n_at=[]; s=1.0; n=len(T)
    idx=0
    for t in uniq:
        at_risk = (T>=t).sum(); d=(T[E==1]==t).sum() if (T==t).any() else ((T==t)&(E==1)).sum()
        # simpler: d = events at t
        d = int(((T==t)&(E==1)).sum())
        s = s*(1-d/at_risk) if at_risk>0 else s
        surv.append(s); n_at.append(at_risk)
    return np.array(uniq), np.array(surv), np.array(n_at)

def logrank(Ta,Ea,Tb,Eb):
    from scipy.stats import chi2
    times = np.sort(np.unique(np.concatenate([Ta[Ea==1], Tb[Eb==1]])))
    O1=E1=V=0.0
    for t in times:
        n1=(Ta>=t).sum(); n2=(Tb>=t).sum(); n=n1+n2
        d1=((Ta==t)&(Ea==1)).sum(); d2=((Tb==t)&(Eb==1)).sum(); d=d1+d2
        if n>1 and d>0:
            e1=d*n1/n
            O1+=d1; E1+=e1
            V+= (n1*n2*d*(n-d))/(n*n*(n-1)) if n>1 else 0
    chi = (O1-E1)**2/V if V>0 else 0
    p = float(1-chi2.cdf(chi,1))
    return float(chi), p

hi = (group=="High"); lo = (group=="Low")
th,sh,nh = km(T[hi],E[hi],T); tl,sl,nl = km(T[lo],E[lo],T)
chi2v, plog = logrank(T[hi],E[hi],T[lo],E[lo])
# HR high vs low (univariable Cox on group)
gnum = hi.astype(int)
res_g = PHReg(T, gnum, E).fit(disp=0)
HR_gl = float(np.exp(res_g.params[0])); p_gl = float(res_g.pvalues[0])
# C-index overall
C = cindex(T,E,risk)
# Descriptive known-status cumulative/dynamic classification AUC at 1/3/5y.
# Cases are observed events by t; controls are patients followed beyond t;
# patients censored on or before t without an event are excluded. This is not an
# IPCW-corrected survival ROC estimator.
def td_auc(t0):
    y = ((T<=t0)&(E==1)).astype(int)
    # controls: T>t0; exclude censored before t0
    mask = ((T<=t0)&(E==1)) | (T>t0)
    if y[mask].sum()==0 or y[mask].sum()==mask.sum(): return float("nan"), 0
    return float(roc_auc_score(y[mask], risk[mask])), int(mask.sum())
aucs = {f"{t}y": td_auc(t) for t in [1,3,5]}
log(f"log-rank chi2={chi2v:.2f} p={plog:.4g}; HR High-vs-Low={HR_gl:.2f} p={p_gl:.4g}; C={C:.3f}; AUCs={aucs}")

# save locked model
locked = {"version":"v1.0 TRAIN-001 2026-09-20","seed":SEED,"genes":sel_genes,"coefs":[float(x) for x in coefs],
 "scaling":{"means":{g:float(means[g]) for g in sel_genes},"sds":{g:float(sds[g]) for g in sel_genes}},
 "cutoff":cutoff,"alpha":float(best_alpha),"n_train":int(len(T)),"n_events":int(E.sum()),
 "C_train":float(C),"HR_high_low":float(HR_gl),"logrank_p":float(plog)}
with open(DER/"locked_model.json","w") as f: json.dump(locked,f,indent=2)
mv.to_csv(DER/"locked_coefficients.tsv",sep="\t",index=False)
# KM plot
plt.figure(figsize=(6,5))
plt.step(np.concatenate([[0],th]), np.concatenate([[1],sh]), where="post", label=f"High (n={hi.sum()})")
plt.step(np.concatenate([[0],tl]), np.concatenate([[1],sl]), where="post", label=f"Low (n={lo.sum()})")
plt.ylim(0,1.02); plt.xlabel("Years"); plt.ylabel("OS probability")
plt.legend(); plt.tight_layout(); plt.savefig(FIG/"train_KM.png",dpi=150); plt.savefig(FIG/"train_KM.pdf")
# summary json
with open(RES/"train_summary.json","w") as f:
    json.dump({"n":int(len(T)),"events":int(E.sum()),"n_uni_sig":int(len(sig)),"n_selected":int(len(sel_genes)),
      "genes":sel_genes,"coefs":[float(x) for x in coefs],"cutoff":cutoff,"C":float(C),
     "HR_high_low":float(HR_gl),"logrank_p":float(plog),"aucs":{k:[float(v[0]) if v[0]==v[0] else None,v[1]] for k,v in aucs.items()}},f,indent=2)
log("WROTE locked_model.json, train_KM.png, train_summary.json")
logf.close()
