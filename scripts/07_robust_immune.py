#!/usr/bin/env python3
"""07_robust_immune.py — RUN-ID: ROBUST-001
Robustness (training adjusted Cox, alternative scaling, subtype stratification) + immune z-mean scores + figure/table source data.
"""
import pandas as pd, numpy as np, json
from pathlib import Path
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROC=ROOT/"data/processed"; META=ROOT/"metadata"; RES=ROOT/"results/raw"; DER=ROOT/"results/derived"; FIG=ROOT/"figures"; TAB=ROOT/"tables"
for d in [RES,FIG,TAB]: d.mkdir(parents=True, exist_ok=True)
logf=open(ROOT/"logs/robust.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()

locked=json.load(open(DER/"locked_model.json")); genes=locked["genes"]; coefs=np.array(locked["coefs"])
kegg=json.load(open(RES/"KEGG_genesets.json")); gs=kegg["genesets"]
# --- training adjusted ---
c425=pd.read_csv(META/"GSE42568_clinical_curated.tsv",sep="\t"); e425=pd.read_csv(PROC/"GSE42568_expr_gene.tsv",sep="\t",index_col=0)
tr=c425[(c425.tissue=="breast cancer")&c425.OS_time_days.notna()].copy()
T=tr.OS_time_days.values/365.25; E=tr.OS_event.values.astype(int)
X=e425.loc[genes,tr.GSM].T; Z=(X-pd.Series(locked["scaling"]["means"]))/pd.Series(locked["scaling"]["sds"])
risk=Z.values@coefs
tr["risk"]=risk
for c in ["age","grade","size","lymph_node_status","er_status"]:
    tr[c]=pd.to_numeric(tr[c],errors="coerce")
adj=tr[["risk","age","grade","size"]].fillna(tr[["risk","age","grade","size"]].median())
ra=PHReg(T,adj.values,E).fit(disp=0)
log(f"Train adjusted (risk+age+grade+size): risk HR={np.exp(ra.params[0]):.2f} p={ra.pvalues[0]:.3g} (n=104, events=35; EPV caution: 35/4=8.75 <10 -> reported as exploratory)")
# alternative scaling: median/IQR
med=X.median(); iqr=(X.quantile(0.75)-X.quantile(0.25)).replace(0,1)
Z2=(X-med)/iqr; risk2=Z2.values@coefs
def cindex(T,E,r):
    n=len(T);cc=tt=0
    for i in range(n):
        if E[i]!=1: continue
        for j in range(n):
            if T[j]>T[i]:
                tt+=1
                if r[j]<r[i]: cc+=1
                elif r[j]==r[i]: cc+=0.5
    return cc/tt if tt else 0.5
log(f"Alternative median/IQR scaling C={cindex(T,E,risk2):.3f} vs z-score C={cindex(T,E,risk):.3f}")
# --- validation subtype stratification ---
c206=pd.read_csv(META/"GSE20685_clinical_curated.tsv",sep="\t")
vr=pd.read_csv(RES/"validation_risk_GSE20685.tsv",sep="\t")
c206=c206.merge(vr[["GSM","risk","group"]],on="GSM")
# subtype column: characteristics 'subtype' (type I-VI)
for st in sorted(c206["subtype"].dropna().unique()):
    s=c206[c206.subtype==st]; hr=np.nan; p=np.nan
    if s.OS_event.sum()>=5 and len(s)>=20:
        try:
            r=PHReg(s.OS_years.values,s.HR_dummy.values if False else (s.group=="High").astype(int).values,s.OS_event.values).fit(disp=0)
            hr=float(np.exp(r.params[0])); p=float(r.pvalues[0])
        except Exception as e: pass
    log(f"20685 {st}: n={len(s)}, deaths={int(s.OS_event.sum())}, HR High-vs-Low={hr} p={p}")
# --- immune z-mean scores (outcome-blind gene sets, computed per cohort) ---
def zmean_scores(expr_gsm, gene_sets):
    # expr: genes x samples; z-score per gene within cohort, then mean over set members present
    mu=expr_gsm.mean(axis=1); sd=expr_gsm.std(axis=1,ddof=0).replace(0,1)
    Zc=(expr_gsm.T-mu).T/sd if False else expr_gsm.sub(mu,axis=0).div(sd,axis=0)
    out={}
    for name,members in gene_sets.items():
        m=[g for g in members if g in expr_gsm.index]
        if len(m)>=5:
            out[name]=Zc.loc[m].mean(axis=0)
    return pd.DataFrame(out)

for gse in ["GSE42568","GSE20685"]:
    expr=pd.read_csv(PROC/f"{gse}_expr_gene.tsv",sep="\t",index_col=0)
    # restrict to tumors with risk
    if gse=="GSE42568":
        clin=pd.read_csv(META/"GSE42568_clinical_curated.tsv",sep="\t"); keep=clin[(clin.tissue=="breast cancer")&clin.OS_time_days.notna()]["GSM"].tolist()
        risk_s=pd.read_csv(RES/"validation_risk_GSE20685.tsv",sep="\t") if False else None
        risks=pd.Series(risk,index=keep)
        groups=pd.Series(np.where(risk>locked["cutoff"],"High","Low"),index=keep)
    else:
        keep=c206["GSM"].tolist(); risks=pd.Series(c206["risk"].values,index=keep); groups=pd.Series(c206["group"].values,index=keep)
    scores=zmean_scores(expr[keep], {k:v for k,v in gs.items()})
    scores["risk"]=risks; scores["group"]=groups.values
    scores.to_csv(RES/f"immune_scores_{gse}.tsv",sep="\t")
    # compare High vs Low per pathway
    rows=[]
    for pw in gs.keys():
        if pw not in scores.columns: continue
        h=scores[scores.group=="High"][pw].values; l=scores[scores.group=="Low"][pw].values
        try:
            _, p = mannwhitneyu(h, l, alternative="two-sided")
            if np.isnan(p): p = 1.0
        except: p = 1.0
        rows.append((pw,float(np.median(h)-np.median(l)),float(p),len(h),len(l)))
    comp=pd.DataFrame(rows,columns=["pathway","median_diff_High_minus_Low","p_MW","n_High","n_Low"])
    from statsmodels.stats.multitest import multipletests
    _,comp["p_adj_BH"],_,_=multipletests(comp["p_MW"],method="fdr_bh")
    comp.to_csv(RES/f"immune_compare_{gse}.tsv",sep="\t",index=False)
    log(f"{gse} immune compare (High vs Low, BH):\n"+comp.to_string(index=False))

# --- checkpoint genes box data ---
CHECKPOINTS=["PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2","IDO1","CD8A","CD8B","GZMA","GZMB","PRF1"]
for gse in ["GSE42568","GSE20685"]:
    expr=pd.read_csv(PROC/f"{gse}_expr_gene.tsv",sep="\t",index_col=0)
    avail=[g for g in CHECKPOINTS if g in expr.index]
    log(f"{gse} checkpoints available: {avail}")
    if gse=="GSE42568":
        keep=tr.GSM.tolist(); grp=np.where(risk>locked["cutoff"],"High","Low")
    else:
        keep=c206.GSM.tolist(); grp=c206.group.values
    dat=expr.loc[avail,keep].T; dat["group"]=grp
    dat.to_csv(RES/f"checkpoint_expr_{gse}.tsv",sep="\t",index=False)
    for g in avail:
        h=dat[dat.group=="High"][g].values; l=dat[dat.group=="Low"][g].values
        try:
            _, p = mannwhitneyu(h, l, alternative="two-sided")
            if np.isnan(p): p = 1.0
        except: p = 1.0
        log(f"{gse} {g}: High med {np.median(h):.2f} vs Low {np.median(l):.2f}, p={p:.3g}")

# --- figures: risk distribution + heatmap data ---
fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].hist(risk,bins=20,alpha=0.7); ax[0].axvline(locked["cutoff"],color="r",ls="--"); ax[0].set_xlabel("RiskScore")
vr2=pd.read_csv(RES/"validation_risk_GSE20685.tsv",sep="\t")
ax[1].hist(vr2.risk.values,bins=20,alpha=0.7); ax[1].axvline(locked["cutoff"],color="r",ls="--"); ax[1].set_xlabel("RiskScore")
plt.tight_layout(); plt.savefig(FIG/"risk_distributions.png",dpi=150); plt.savefig(FIG/"risk_distributions.pdf")
# heatmap source: top 14 genes z-scored, ordered by risk (train)
Zh=(X-pd.Series(locked["scaling"]["means"]))/pd.Series(locked["scaling"]["sds"])
order=np.argsort(risk)
Zh.iloc[order].to_csv(RES/"heatmap_train_14g_z.tsv",sep="\t")
log("WROTE robust + immune outputs, risk_distributions.png")
logf.close()
