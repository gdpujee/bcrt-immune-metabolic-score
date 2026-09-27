#!/usr/bin/env python3
"""12_clinical_value.py — RUN-ID: CLINIC-001
Incremental clinical value (ΔC clinical vs +risk), nomogram construction data, calibration plot data,
DCA (5y, known-status simplification disclosed), treatment-interaction screens. Exploratory unless noted.
Inputs: metadata curated TSVs, results/raw/{validation_risk,rnaseq_risk}.tsv
Outputs: results/raw/clinical_value.json, results/raw/nomogram_points.tsv, results/raw/calibration_*.tsv,
results/raw/dca_*.tsv, results/raw/treatment_interaction.tsv, figures/Nomogram.png, figures/Calibration_SCANB.png, figures/DCA_SCANB.png
Seed 42. No model refitting on test beyond evaluation Coxs (nomogram = presentation of SCAN-B multivariable fit).
"""
import pandas as pd, numpy as np, json
from pathlib import Path
from statsmodels.duration.hazard_regression import PHReg
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT/"results/raw"; FIG = ROOT/"figures"; LOGS = ROOT/"logs"; META = ROOT/"metadata"
logf = open(LOGS/"clinical_value.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()
rng = np.random.default_rng(42)

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

def dC_CI(T,E,X_base,X_full,B=500):
    # ΔC with bootstrap CI (same resamples)
    n=len(T); obs=None; diffs=[]
    rb=PHReg(T,X_base,E).fit(disp=0); rf=PHReg(T,X_full,E).fit(disp=0)
    cb=cindex(T,E,X_base@np.asarray(rb.params).ravel()); cf=cindex(T,E,X_full@np.asarray(rf.params).ravel())
    obs=cf-cb
    for b in range(B):
        idx=rng.integers(0,n,n)
        try:
            rb2=PHReg(T[idx],X_base[idx],E[idx]).fit(disp=0); rf2=PHReg(T[idx],X_full[idx],E[idx]).fit(disp=0)
            cb2=cindex(T[idx],E[idx],X_base[idx]@np.asarray(rb2.params).ravel())
            cf2=cindex(T[idx],E[idx],X_full[idx]@np.asarray(rf2.params).ravel())
            diffs.append(cf2-cb2)
        except Exception: pass
    diffs=np.array(diffs)
    return obs, float(np.percentile(diffs,2.5)), float(np.percentile(diffs,97.5)), float(cb), float(cf)

out={}
# ---------- GSE20685 ΔC ----------
c206=pd.read_csv(META/"GSE20685_clinical_curated.tsv",sep="\t")
vr=pd.read_csv(RES/"validation_risk_GSE20685.tsv",sep="\t")
c206=c206.merge(vr[["GSM","risk"]],on="GSM")
for cc in ["age_at_diagnosis","t_stage","n_stage"]: c206[cc]=pd.to_numeric(c206[cc],errors="coerce")
T=c206.follow_up_duration_years.values.astype(float); E=c206.event_death.values.astype(int)
Xb=c206[["age_at_diagnosis","t_stage","n_stage"]].fillna(c206[["age_at_diagnosis","t_stage","n_stage"]].median()).values
Xf=np.column_stack([Xb, c206.risk.values])
obs,lo,hi,cb,cf=dC_CI(T,E,Xb,Xf)
log(f"20685 ΔC (+risk over age+T+N): {obs:.3f} [{lo:.3f},{hi:.3f}] (base {cb:.3f} → full {cf:.3f})")
out["GSE20685_dC"]={"delta":obs,"CI":[lo,hi],"base_C":cb,"full_C":cf}
# ---------- SCAN-B ΔC ----------
rs=pd.read_csv(RES/"rnaseq_risk_GSE96058.tsv",sep="\t")
cu=pd.read_csv(META/"GSE96058_clinical_raw.tsv",sep="\t")
cu["is_repl"]=cu.title.str.contains("repl",na=False); cu=cu[~cu.is_repl].copy()
for cc in ["age_at_diagnosis","er_status","her2_status"]: cu[cc]=pd.to_numeric(cu[cc],errors="coerce")
cu=cu.merge(rs[["sample","risk"]],left_on="title",right_on="sample")
T2=cu.overall_survival_days.values.astype(float)/365.25; E2=cu.overall_survival_event.values.astype(int)
Xb2=cu[["age_at_diagnosis","er_status","her2_status"]].fillna(cu[["age_at_diagnosis","er_status","her2_status"]].median()).values
Xf2=np.column_stack([Xb2, cu.risk.values])
obs2,lo2,hi2,cb2,cf2=dC_CI(T2,E2,Xb2,Xf2)
log(f"SCANB ΔC (+risk over age+ER+HER2): {obs2:.3f} [{lo2:.3f},{hi2:.3f}] (base {cb2:.3f} → full {cf2:.3f})")
out["SCANB_dC"]={"delta":obs2,"CI":[lo2,hi2],"base_C":cb2,"full_C":cf2}
# ---------- Nomogram (SCAN-B multivariable fit, presentation only) ----------
Xn=np.column_stack([cu.risk.values, cu[["age_at_diagnosis","er_status","her2_status"]].fillna(cu[["age_at_diagnosis","er_status","her2_status"]].median()).values])
cols=["risk","age","ER","HER2"]
rn=PHReg(T2,Xn,E2).fit(disp=0)
coefs=np.asarray(rn.params).ravel()
# points: scale each variable's (value-min)*|coef| to 0-100 over observed range
pts={}
for j,c in enumerate(cols):
    v=Xn[:,j]; rng_j=v.max()-v.min()
    pts[c]={"coef":float(coefs[j]),"min":float(v.min()),"max":float(v.max()),
            "points_per_unit":float(abs(coefs[j])/rng_j*100) if rng_j else 0.0}
pd.DataFrame(pts).T.to_csv(RES/"nomogram_points.tsv",sep="\t")
# baseline survival at 5y (Breslow approx via KM of prognostic index? use overall KM@5y + PI centering)
order=np.argsort(T2); Tt=T2[order]; Ee=E2[order]
S5=1.0
for t in np.sort(np.unique(Tt[Ee==1])):
    if t<=5:
        at=(Tt>=t).sum(); d=int(((Tt==t)&(Ee==1)).sum()); S5*=(1-d/at) if at else 1
mean_pi=float((Xn@coefs).mean())
log(f"Nomogram base: S0(5y)={S5:.3f} at mean PI {mean_pi:.3f}; coefs {dict(zip(cols,[round(float(x),4) for x in coefs]))}")
out["nomogram"]={"S0_5y":S5,"mean_PI":mean_pi,"coefs":{c:float(coefs[j]) for j,c in enumerate(cols)}}
# figure: nomogram schematic (points axes)
fig,ax=plt.subplots(figsize=(8,4)); ax.axis("off")
ax.text(0.5,0.9,"Nomogram (SCAN-B fit, exploratory presentation — not independently validated)",ha="center",fontsize=10,weight="bold")
y=0.7
for c in cols:
    ax.text(0.05,y,f"{c}: {pts[c]['min']:.2f}–{pts[c]['max']:.2f}  (coef {pts[c]['coef']:.3f})",fontsize=9)
    y-=0.12
ax.text(0.05,0.15,f"Total points → 5y OS via S0={S5:.3f} at mean PI (see calibration plot)",fontsize=8)
plt.tight_layout(); plt.savefig(FIG/"Nomogram.png",dpi=150); plt.savefig(FIG/"Nomogram.pdf")
# ---------- Calibration plot data (SCAN-B deciles, 5y KM observed vs Cox predicted) ----------
pi=Xn@coefs
# predicted 5y per person: S0^exp(pi-mean_pi)
pred=np.power(S5, np.exp(pi-mean_pi))
cu["pred5"]=pred
cu["dec"]=pd.qcut(cu.risk,10,labels=False,duplicates="drop")
rows=[]
for d in sorted(cu.dec.dropna().unique()):
    s=cu[cu.dec==d]; tt=s.overall_survival_days.values/365.25; ee=s.overall_survival_event.values.astype(int)
    o=np.argsort(tt); tt=tt[o]; ee=ee[o]; surv=1.0
    for t in np.sort(np.unique(tt[ee==1])):
        if t<=5:
            at=(tt>=t).sum(); dd=int(((tt==t)&(ee==1)).sum()); surv*=(1-dd/at) if at else 1
    rows.append({"decile":int(d),"n":len(s),"deaths":int(s.overall_survival_event.sum()),"mean_pred_5y":float(s.pred5.mean()),"KM_obs_5y":float(surv)})
cal=pd.DataFrame(rows); cal.to_csv(RES/"calibration_SCANB_deciles.tsv",sep="\t",index=False)
log("Calibration deciles:\n"+cal.to_string(index=False))
plt.figure(figsize=(5,5))
plt.scatter(cal.mean_pred_5y, cal.KM_obs_5y)
plt.plot([0.6,1.0],[0.6,1.0],"--",lw=1); plt.xlabel("Predicted 5y OS (Cox)"); plt.ylabel("KM observed 5y OS")
plt.tight_layout(); plt.savefig(FIG/"Calibration_SCANB.png",dpi=150); plt.savefig(FIG/"Calibration_SCANB.pdf")
# ---------- DCA 5y (known-status simplification disclosed) ----------
known=cu[((cu.overall_survival_days<=5*365.25)&(cu.overall_survival_event==1))|(cu.overall_survival_days>5*365.25)].copy()
y=((known.overall_survival_days<=5*365.25)&(known.overall_survival_event==1)).astype(int).values
p=known.pred5.values  # prob of surviving 5y → risk of death = 1-p
rd=1-p
ths=np.linspace(0.05,0.5,19); nb_model=[]; nb_all=[]
prev=y.mean()
for th in ths:
    tp=((rd>=th)&(y==1)).sum(); fp=((rd>=th)&(y==0)).sum(); n=len(y)
    nb_model.append(tp/n - fp/n*th/(1-th))
    nb_all.append(prev - (1-prev)*th/(1-th))
dca=pd.DataFrame({"threshold":ths,"NB_model":nb_model,"NB_all":nb_all,"NB_none":0.0,"n":len(y),"prev":prev})
dca.to_csv(RES/"dca_SCANB_5y.tsv",sep="\t",index=False)
log(f"DCA 5y (known-status n={len(y)}, prev death={prev:.3f}; censored<5y excluded — simplification disclosed)")
plt.figure(figsize=(6,4))
plt.plot(dca.threshold,dca.NB_model,label="Risk model"); plt.plot(dca.threshold,dca.NB_all,label="Treat all"); plt.axhline(0,ls="--",lw=1,label="Treat none")
plt.xlabel("Threshold (5y death prob)"); plt.ylabel("Net benefit"); plt.legend(fontsize=8)
plt.tight_layout(); plt.savefig(FIG/"DCA_SCANB.png",dpi=150); plt.savefig(FIG/"DCA_SCANB.pdf")
# ---------- Treatment interaction ----------
ti_rows=[]
# SCAN-B chemo + endocrine (binary 0/1, NA→exclude)
for var in ["chemo_treated","endocrine_treated"]:
    s=cu[cu[var].isin([0,1])].copy()
    s["tr"]=pd.to_numeric(s[var]); z=(s.risk.values-s.risk.mean())/s.risk.std()
    X=np.column_stack([z, s.tr.values, z*s.tr.values])
    try:
        r=PHReg(s.overall_survival_days.values/365.25,X,s.overall_survival_event.values).fit(disp=0)
        ti_rows.append({"cohort":"SCANB","var":var,"n":len(s),"events":int(s.overall_survival_event.sum()),"HR_interaction":float(np.exp(r.params[2])),"p_interaction":float(r.pvalues[2]),"HR_risk_main":float(np.exp(r.params[0]))})
    except Exception as e: ti_rows.append({"cohort":"SCANB","var":var,"note":f"fit failed: {e}"})
# GSE20685 adjuvant_chemotherapy yes/no (5 unknowns excluded with disclosure, matching SCAN-B isin rule)
c206=c206[c206.adjuvant_chemotherapy.isin(["yes","no"])].copy()
c206["tr"]=(c206.adjuvant_chemotherapy=="yes").astype(int)
s=c206.copy(); z=(s.risk.values-s.risk.mean())/s.risk.std()
X=np.column_stack([z, s.tr.values, z*s.tr.values])
try:
    r=PHReg(s.follow_up_duration_years.values,X,s.event_death.values).fit(disp=0)
    ti_rows.append({"cohort":"GSE20685","var":"adjuvant_chemotherapy","n":len(s),"events":int(s.event_death.sum()),"HR_interaction":float(np.exp(r.params[2])),"p_interaction":float(r.pvalues[2]),"HR_risk_main":float(np.exp(r.params[0]))})
except Exception as e: ti_rows.append({"cohort":"GSE20685","var":"adjuvant_chemotherapy","note":f"fit failed: {e}"})
ti=pd.DataFrame(ti_rows); ti.to_csv(RES/"treatment_interaction.tsv",sep="\t",index=False)
log("Treatment interaction:\n"+ti.to_string(index=False))
json.dump(out,open(RES/"clinical_value.json","w"),indent=2)
log("WROTE clinical_value outputs (CLINIC-001)")
logf.close()
