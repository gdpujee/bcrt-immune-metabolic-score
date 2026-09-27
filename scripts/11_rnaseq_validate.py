#!/usr/bin/env python3
"""11_rnaseq_validate.py — RUN-ID: RNASEQ-001
Locked cross-platform application of TRAIN-001 14-gene model to GSE96058 SCAN-B RNA-seq (3273 unique, repl excluded).
No retraining, no cutoff tuning. Primary: continuous per-SD Cox. Secondary: locked microarray cutoff.
Inputs: data/raw/GSE96058_gene_expression.csv.gz, metadata/GSE96058_clinical_raw.tsv, results/derived/locked_model.json
Outputs: results/raw/rnaseq_*.json/tsv, figures/RNASEQ_KM.png, logs/rnaseq.log
"""
import gzip, json
from pathlib import Path
import pandas as pd, numpy as np
from statsmodels.duration.hazard_regression import PHReg
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu, chi2
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT/"results/raw"; DER = ROOT/"results/derived"; FIG = ROOT/"figures"; LOGS = ROOT/"logs"; META = ROOT/"metadata"
logf = open(LOGS/"rnaseq.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()
rng = np.random.default_rng(42)

locked = json.load(open(DER/"locked_model.json"))
genes = locked["genes"]; coefs = np.array(locked["coefs"])
means = pd.Series(locked["scaling"]["means"]); sds = pd.Series(locked["scaling"]["sds"]); cutoff = locked["cutoff"]
log(f"Locked {len(genes)} genes, cutoff {cutoff:.4f} (microarray scale; applied to RNA-seq with frozen params)")

# load 14 genes, exclude repl cols
with gzip.open(ROOT/"data/raw/GSE96058_gene_expression.csv.gz","rt") as f:
    header = [c.strip('"') for c in f.readline().rstrip().split(",")]
    keep_idx = [0] + [i for i,c in enumerate(header) if i>0 and "repl" not in c]
    keep_cols = [header[i] for i in keep_idx]
    log(f"Columns total {len(header)-1}, unique {len(keep_cols)-1} (excluded {len(header)-len(keep_cols)} repl)")
    rows = {}
    for line in f:
        g = line.split(",",1)[0].strip().strip('"')
        if g in genes:
            parts = line.rstrip().split(",")
            rows[g] = [float(parts[i]) for i in keep_idx[1:]]
expr = pd.DataFrame(rows, index=[c for c in keep_cols[1:]]).T  # genes x samples(F-titles)
log(f"Expr {expr.shape}, genes missing: {[g for g in genes if g not in expr.index]}")
assert all(g in expr.index for g in genes), "missing locked genes — cannot claim validation"
# clinical map via title (non-repl)
clin = pd.read_csv(META/"GSE96058_clinical_raw.tsv", sep="\t")
clin["is_repl"] = clin.title.str.contains("repl", na=False)
cu = clin[~clin.is_repl].copy()
cu["OS_years"] = pd.to_numeric(cu.overall_survival_days, errors="coerce")/365.25
cu["OS_event"] = pd.to_numeric(cu.overall_survival_event, errors="coerce")
cu = cu[cu.OS_years.notna() & cu.OS_event.notna()].copy()
log(f"Unique patients with OS: {len(cu)} (deaths={int((cu.OS_event==1).sum())})")
# align
common = [t for t in expr.columns if t in set(cu.title)]
log(f"Aligned samples: {len(common)}")
cu = cu.set_index("title").loc[common]
X = expr[common].T.loc[:, genes]
Z = (X - means)/sds
risk = Z.values @ coefs
cu["risk"] = risk
T = cu.OS_years.values.astype(float); E = cu.OS_event.values.astype(int)

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

def logrank(Ta,Ea,Tb,Eb):
    times=np.sort(np.unique(np.concatenate([Ta[Ea==1],Tb[Eb==1]])))
    O1=E1=V=0.0
    for t in times:
        n1=(Ta>=t).sum(); n2=(Tb>=t).sum(); n=n1+n2
        d1=((Ta==t)&(Ea==1)).sum(); d2=((Tb==t)&(Eb==1)).sum(); d=d1+d2
        if n>1 and d>0: E1+=d*n1/n; O1+=d1; V+=(n1*n2*d*(n-d))/(n*n*(n-1))
    ch=(O1-E1)**2/V if V>0 else 0
    return float(ch), float(1-chi2.cdf(ch,1))

# primary continuous
sd = risk.std()
rc = PHReg(T,(risk/sd).reshape(-1,1),E).fit(disp=0)
hr_c=float(np.exp(rc.params[0])); ci_c=[float(np.exp(rc.params[0]-1.96*rc.bse[0])),float(np.exp(rc.params[0]+1.96*rc.bse[0]))]
log(f"PRIMARY continuous per-SD: HR={hr_c:.2f} [{ci_c[0]:.2f},{ci_c[1]:.2f}] p={float(rc.pvalues[0]):.3g}")
# secondary locked binary
grp = np.where(risk>cutoff,"High","Low")
hi=(grp=="High")
rb=PHReg(T,hi.astype(int),E).fit(disp=0)
hr_b=float(np.exp(rb.params[0])); ci_b=[float(np.exp(rb.params[0]-1.96*rb.bse[0])),float(np.exp(rb.params[0]+1.96*rb.bse[0]))]
_,plog=logrank(T[hi],E[hi],T[~hi],E[~hi])
C=cindex(T,E,risk)
log(f"SECONDARY locked cutoff {cutoff:.3f}: High={hi.sum()}/Low={(~hi).sum()}, HR={hr_b:.2f} [{ci_b[0]:.2f},{ci_b[1]:.2f}] p={float(rb.pvalues[0]):.3g}, logrank p={plog:.3g}, C={C:.3f}")
# bootstrap CIs
B=500; Cb=[]; A3=[];A5=[]
n=len(T)
for b in range(B):
    idx=rng.integers(0,n,n)
    Cb.append(cindex(T[idx],E[idx],risk[idx]))
    for arr,t0 in [(A3,3),(A5,5)]:
        y=((T[idx]<=t0)&(E[idx]==1)).astype(int); mask=((T[idx]<=t0)&(E[idx]==1))|(T[idx]>t0)
        try: arr.append(roc_auc_score(y[mask],risk[idx][mask]))
        except: arr.append(np.nan)
def ci(a):
    a=np.array(a); a=a[~np.isnan(a)]; return [float(np.percentile(a,2.5)),float(np.percentile(a,97.5))]
log(f"C 95% [{ci(Cb)[0]:.3f},{ci(Cb)[1]:.3f}]; AUC3y [{ci(A3)[0]:.3f},{ci(A3)[1]:.3f}]; AUC5y [{ci(A5)[0]:.3f},{ci(A5)[1]:.3f}]")
# calibration slope + quartiles
rs=PHReg(T,risk.reshape(-1,1),E).fit(disp=0)
log(f"Calib slope {float(rs.params[0]):.3f} SE {float(rs.bse[0]):.3f}")
cu["q"]=pd.qcut(cu.risk,4,labels=["Q1","Q2","Q3","Q4"])
for q in ["Q1","Q2","Q3","Q4"]:
    s=cu[cu.q==q]; tt=s.OS_years.values; ee=s.OS_event.values.astype(int); order=np.argsort(tt); tt=tt[order]; ee=ee[order]
    surv=1.0
    for t in np.sort(np.unique(tt[ee==1])):
        if t<=5:
            at=(tt>=t).sum(); d=int(((tt==t)&(ee==1)).sum()); surv*=(1-d/at) if at else 1
    log(f"{q}: n={len(s)} deaths={int(s.OS_event.sum())} KM@5y={surv:.3f} mean_risk={s.risk.mean():.2f}")
# adjusted (age + ER + HER2 + size proxies available in SCAN-B)
for c in ["age_at_diagnosis","er_status","her2_status","tumor_size"]:
    cu[c]=pd.to_numeric(cu[c],errors="coerce")
adj=cu[["risk","age_at_diagnosis","er_status","her2_status"]].fillna(cu[["risk","age_at_diagnosis","er_status","her2_status"]].median())
ra=PHReg(T,adj.values,E).fit(disp=0)
log(f"Adjusted (risk+age+ER+HER2): risk HR={float(np.exp(ra.params[0])):.2f} p={float(ra.pvalues[0]):.3g}")
# PAM50 stratification (powered: LumA 1657, LumB 729, Her2/Basal available)
for st in sorted(cu.pam50_subtype.dropna().unique()):
    s=cu[cu.pam50_subtype==st]
    if len(s)>=50 and s.OS_event.sum()>=10:
        try:
            r=PHReg(s.OS_years.values,(s.risk.values/s.risk.std()).reshape(-1,1),s.OS_event.values).fit(disp=0)
            log(f"PAM50 {st}: n={len(s)} deaths={int(s.OS_event.sum())} per-SD HR={float(np.exp(r.params[0])):.2f} p={float(r.pvalues[0]):.3g}")
        except Exception as e: log(f"PAM50 {st}: fit failed {e}")
    else:
        log(f"PAM50 {st}: n={len(s)} deaths={int(s.OS_event.sum())} — underpowered, N/events only")
# instrument batch
for inst in cu.instrument_model.dropna().unique():
    s=cu[cu.instrument_model==inst]
    log(f"Instrument {inst}: n={len(s)} deaths={int(s.OS_event.sum())} mean_risk={s.risk.mean():.2f}")
# 5y truncation
T5=np.minimum(T,5); E5=((T<=5)&(E==1)).astype(int)
r5=PHReg(T5,(risk/sd).reshape(-1,1),E5).fit(disp=0)
log(f"5y-truncated continuous HR={float(np.exp(r5.params[0])):.2f} p={float(r5.pvalues[0]):.3g} (events={int(E5.sum())})")

# save
cu.reset_index()[["title","GSM","risk","OS_years","OS_event","pam50_subtype","instrument_model"]].rename(columns={"title":"sample","GSM":"GSM"}).to_csv(RES/"rnaseq_risk_GSE96058.tsv",sep="\t",index=False)
json.dump({"n":int(len(T)),"events":int(E.sum()),"hr_cont":[hr_c]+ci_c+[float(rc.pvalues[0])],
 "hr_bin":[hr_b]+ci_b+[float(rb.pvalues[0])],"logrank_p":plog,"C":C,"C_CI":ci(Cb),
 "AUC3_CI":ci(A3),"AUC5_CI":ci(A5),"slope":[float(rs.params[0]),float(rs.bse[0])],
 "adj_hr":float(np.exp(ra.params[0])),"adj_p":float(ra.pvalues[0])},open(RES/"rnaseq_summary.json","w"),indent=2)
# KM plot
def km(T,E):
    o=np.argsort(T); T=T[o];E=E[o]; u=np.sort(np.unique(T[E==1])); s=[];v=1.0
    for t in u:
        at=(T>=t).sum(); d=int(((T==t)&(E==1)).sum()); v=v*(1-d/at) if at else v; s.append(v)
    return u,np.array(s)
th,sh=km(T[hi],E[hi]); tl,sl=km(T[~hi],E[~hi])
plt.figure(figsize=(6,5))
plt.step(np.concatenate([[0],th]),np.concatenate([[1],sh]),where="post",label=f"High (n={hi.sum()})")
plt.step(np.concatenate([[0],tl]),np.concatenate([[1],sl]),where="post",label=f"Low (n={(~hi).sum()})")
plt.ylim(0,1.02);plt.xlabel("Years");plt.ylabel("OS probability")
plt.legend();plt.tight_layout();plt.savefig(FIG/"RNASEQ_KM.png",dpi=150); plt.savefig(FIG/"RNASEQ_KM.pdf")
log("WROTE rnaseq outputs + RNASEQ_KM.png")
logf.close()
