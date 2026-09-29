#!/usr/bin/env python3
"""09_corrective.py — RUN-ID: CORRECT-001
Addresses Reviewers A+B computable issues: bootstrap CIs, calibration, continuous primary,
validation-median sensitivity, BH checkpoints, pathway-minus-risk-genes, loss tables, 5y truncation.
"""
import pandas as pd, numpy as np, json
from pathlib import Path
from statsmodels.duration.hazard_regression import PHReg
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
RES=ROOT/"results/raw"; DER=ROOT/"results/derived"; LOGS=ROOT/"logs"
logf=open(LOGS/"corrective.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()
rng=np.random.default_rng(42)

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

locked=json.load(open(DER/"locked_model.json")); genes=locked["genes"]; coefs=np.array(locked["coefs"])
# --- validation continuous primary + CIs ---
vr=pd.read_csv(RES/"validation_risk_GSE20685.tsv",sep="\t")
T=vr.OS_years.values; E=vr.OS_event.values; risk=vr.risk.values
# continuous Cox (per 1-SD risk)
sd=risk.std()
res_c=PHReg(T,(risk/sd).reshape(-1,1),E).fit(disp=0)
hr_c=float(np.exp(res_c.params[0])); lo_c=float(np.exp(res_c.params[0]-1.96*res_c.bse[0])); hi_c=float(np.exp(res_c.params[0]+1.96*res_c.bse[0]))
log(f"VALID continuous per-SD: HR={hr_c:.2f} 95%CI [{lo_c:.2f},{hi_c:.2f}] p={float(res_c.pvalues[0]):.3g}")
# binary HR CI (already has HR, add CI)
hi=(vr.group=="High").astype(int).values
res_b=PHReg(T,hi,E).fit(disp=0)
hr_b=float(np.exp(res_b.params[0])); lo_b=float(np.exp(res_b.params[0]-1.96*res_b.bse[0])); hi_b=float(np.exp(res_b.params[0]+1.96*res_b.bse[0]))
log(f"VALID binary High-vs-Low: HR={hr_b:.2f} 95%CI [{lo_b:.2f},{hi_b:.2f}] p={float(res_b.pvalues[0]):.3g}")
# Bootstrap CIs for Harrell's C and the descriptive known-status cumulative/
# dynamic AUCs below (not IPCW-corrected; censored before/at t are excluded).
# 1000 resamples, seed 42.
B=1000; Cb=[]; A1=[];A3=[];A5=[]
n=len(T)
for b in range(B):
    idx=rng.integers(0,n,n)
    Cb.append(cindex(T[idx],E[idx],risk[idx]))
    for arr,t0 in [(A1,1),(A3,3),(A5,5)]:
        y=((T[idx]<=t0)&(E[idx]==1)).astype(int); mask=((T[idx]<=t0)&(E[idx]==1))|(T[idx]>t0)
        try: arr.append(roc_auc_score(y[mask],risk[idx][mask]))
        except: arr.append(np.nan)
def ci(a):
    a=np.array(a); a=a[~np.isnan(a)]
    return float(np.percentile(a,2.5)), float(np.percentile(a,97.5))
log(f"C 95%CI bootstrap: [{ci(Cb)[0]:.3f},{ci(Cb)[1]:.3f}] (point {cindex(T,E,risk):.3f})")
log(f"AUC1y 95%CI: [{ci(A1)[0]:.3f},{ci(A1)[1]:.3f}]")
log(f"AUC3y 95%CI: [{ci(A3)[0]:.3f},{ci(A3)[1]:.3f}]")
log(f"AUC5y 95%CI: [{ci(A5)[0]:.3f},{ci(A5)[1]:.3f}]")
# Raw locked-score Cox coefficient and descriptive 5y KM by score quartile.
res_slope=PHReg(T,risk.reshape(-1,1),E).fit(disp=0)
log(f"Raw-score Cox coefficient (validation cohort): {float(res_slope.params[0]):.3f} SE {float(res_slope.bse[0]):.3f}")
# 5y KM observed by risk quartile
vr["q"]=pd.qcut(vr.risk,4,labels=["Q1","Q2","Q3","Q4"])
for q in ["Q1","Q2","Q3","Q4"]:
    s=vr[vr.q==q]; km5=((s.OS_years>5).sum()+( ((s.OS_years<=5)&(s.OS_event==1)).sum()==0 ))/len(s) if len(s) else np.nan
    # proper KM at 5y:
    tt=s.OS_years.values; ee=s.OS_event.values; order=np.argsort(tt); tt=tt[order]; ee=ee[order]
    surv=1.0
    for t in np.sort(np.unique(tt[ee==1])):
        if t<=5:
            at=(tt>=t).sum(); d=int(((tt==t)&(ee==1)).sum()); surv*=(1-d/at) if at else 1
    log(f"Calib {q}: n={len(s)}, deaths={int(s.OS_event.sum())}, KM@5y={surv:.3f}, mean risk={s.risk.mean():.3f}")
# 5y truncation sensitivity (common horizon): censor all at 5y
T5=np.minimum(T,5); E5=((T<=5)&(E==1)).astype(int)
res5=PHReg(T5,hi,E5).fit(disp=0)
log(f"5y-truncated binary HR={float(np.exp(res5.params[0])):.2f} p={float(res5.pvalues[0]):.3g} (events={int(E5.sum())})")
# validation-median exploratory
med=float(vr.risk.median())
g2=np.where(risk>med,"High","Low")
from scipy.stats import chi2
def logrank(Ta,Ea,Tb,Eb):
    times=np.sort(np.unique(np.concatenate([Ta[Ea==1],Tb[Eb==1]])))
    O1=E1=V=0.0
    for t in times:
        n1=(Ta>=t).sum(); n2=(Tb>=t).sum(); n=n1+n2
        d1=((Ta==t)&(Ea==1)).sum(); d2=((Tb==t)&(Eb==1)).sum(); d=d1+d2
        if n>1 and d>0: E1+=d*n1/n; O1+=d1; V+=(n1*n2*d*(n-d))/(n*n*(n-1))
    chi=(O1-E1)**2/V if V>0 else 0
    return float(chi), float(1-chi2.cdf(chi,1))
h2=(g2=="High")
_,p2=logrank(T[h2],E[h2],T[~h2],E[~h2])
res2=PHReg(T,h2.astype(int),E).fit(disp=0)
log(f"Exploratory validation-median cutoff {med:.3f}: High={(h2).sum()}/Low={(~h2).sum()}, HR={float(np.exp(res2.params[0])):.2f} p={float(res2.pvalues[0]):.3g}, logrank p={p2:.3g} (NOT locked; exploratory)")

# --- BH checkpoints ---
chk=pd.read_csv(RES/"checkpoint_expr_GSE20685.tsv",sep="\t")
ps=[]; names=[]
for g in [c for c in chk.columns if c!="group"]:
    h=chk[chk.group=="High"][g].values; l=chk[chk.group=="Low"][g].values
    try:
        _, p = mannwhitneyu(h, l, alternative="two-sided")
        if np.isnan(p): p = 1.0
    except: p = 1.0
    ps.append(p); names.append(g)
_,padj,_,_=multipletests(ps,method="fdr_bh")
for g,p,pa in zip(names,ps,padj): log(f"Checkpoint {g}: raw p={p:.3g}, BH adj={pa:.3g}")
pd.DataFrame({"gene":names,"p_raw":ps,"p_adj_BH":padj}).to_csv(RES/"checkpoint_BH_GSE20685.tsv",sep="\t",index=False)

# --- pathway minus 14 risk genes ---
kegg=json.load(open(RES/"KEGG_genesets.json"))["genesets"]
riskset=set(genes)
for name,members in kegg.items():
    overlap=set(members)&riskset
    if overlap: log(f"KEGG {name}: contains risk genes {sorted(overlap)} (will be excluded in sensitivity)")
# recompute z-mean excluding risk genes (validation)
e206=pd.read_csv(ROOT/"data/processed/GSE20685_expr_gene.tsv",sep="\t",index_col=0)
keep=vr.GSM.tolist(); groups=vr.group.values
mu=e206[keep].mean(axis=1); sd=e206[keep].std(axis=1,ddof=0).replace(0,1)
Zc=e206[keep].sub(mu,axis=0).div(sd,axis=0)
rows=[]
for name,members in kegg.items():
    m=[g for g in members if g in e206.index and g not in riskset]
    if len(m)<5: rows.append((name,0,len(m),1.0,1.0)); continue
    sc=Zc.loc[m].mean(axis=0).values
    h=sc[groups=="High"]; l=sc[groups=="Low"]
    try:
        _, p = mannwhitneyu(h, l, alternative="two-sided")
        if np.isnan(p): p = 1.0
    except: p = 1.0
    rows.append((name,float(np.median(h)-np.median(l)),len(m),float(p)))
comp=pd.DataFrame(rows,columns=["pathway","median_diff","n_genes_excl_risk","p_MW"])
_,comp["p_adj_BH"],_,_=multipletests(comp.p_MW,method="fdr_bh")
log("Pathway EXCLUDING 14 risk genes (validation):\n"+comp[["pathway","median_diff","n_genes_excl_risk","p_MW","p_adj_BH"]].to_string(index=False))
comp.to_csv(RES/"immune_compare_excl_risk_GSE20685.tsv",sep="\t",index=False)

# --- loss tables (A2) ---
pool=pd.read_csv(RES/"candidate_pool.tsv",sep="\t")
inter=set(pd.read_csv(ROOT/"data/processed/GSE45827_expr_gene.tsv",sep="\t",usecols=["symbol"])["symbol"])
log(f"Loss by 45827 pre-filter: pool 297 -> in-intersection {(pool.symbol.isin(inter)).sum()} (lost {(~pool.symbol.isin(inter)).sum()}); 14 risk genes in intersection: {sum(g in inter for g in genes)}/14")
# KEGG overlap audit (A5): metabolic vs immune disjoint by KEGG design?
met=set().union(*[set(kegg[k]) for k in ["Glycolysis_Gluconeogenesis","Fatty_acid_degradation","Fatty_acid_biosynthesis","Oxidative_phosphorylation","Pentose_phosphate","Glutathione_metabolism"]])
imm=set().union(*[set(kegg[k]) for k in ["Cytokine_cytokine_receptor","Chemokine_signaling","Antigen_processing_presentation","NK_cytotoxicity","T_cell_receptor"]])
log(f"KEGG met∩imm = {len(met&imm)} (0 expected: KEGG metabolic vs immune pathway gene sets are disjoint by curation; 0 overlap is NOT an artefact)")

# --- save corrective summary ---
with open(RES/"corrective_summary.json","w") as f:
    json.dump({"valid_continuous_perSD_HR":[hr_c,lo_c,hi_c,float(res_c.pvalues[0])],
     "valid_binary_HR_CI":[hr_b,lo_b,hi_b,float(res_b.pvalues[0])],
     "C_CI":list(ci(Cb)),"AUC1_CI":list(ci(A1)),"AUC3_CI":list(ci(A3)),"AUC5_CI":list(ci(A5)),
     "raw_score_coefficient":[float(res_slope.params[0]),float(res_slope.bse[0])],
     "trunc5y_HR":[float(np.exp(res5.params[0])),float(res5.pvalues[0])],
     "valid_median_exploratory":[med,float(np.exp(res2.params[0])),float(res2.pvalues[0]),float(p2)]},f,indent=2)
log("WROTE corrective_summary.json, checkpoint_BH, immune_compare_excl_risk")
logf.close()
