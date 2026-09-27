#!/usr/bin/env python3
"""08_figures_tables.py — RUN-ID: FIG-001
Generate all manuscript figures/tables from raw results (no new modeling).
"""
import pandas as pd, numpy as np, json
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import seaborn as sns

ROOT=Path(__file__).resolve().parents[1]
RES=ROOT/"results/raw"; DER=ROOT/"results/derived"; FIG=ROOT/"figures"; TAB=ROOT/"tables"
for d in [FIG,TAB]: d.mkdir(parents=True,exist_ok=True)
logf=open(ROOT/"logs/figures.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()

# --- Fig2 volcano 42568 ---
de=pd.read_csv(RES/"DE_GSE42568_tumor_vs_normal.tsv",sep="\t")
de["neglog10FDR"]=-np.log10(de.p_adj_BH.clip(lower=1e-300))
plt.figure(figsize=(6,5))
sig=(de.p_adj_BH<0.05)&(de.log2FC_median_diff.abs()>1)
plt.scatter(de.log2FC_median_diff, de.neglog10FDR, s=3, alpha=0.4, label=f"ns (n={(~sig).sum()})")
plt.scatter(de[sig].log2FC_median_diff, de[sig].neglog10FDR, s=4, alpha=0.7, label=f"FDR<0.05 & |logFC|>1 (n={sig.sum()})")
plt.axhline(-np.log10(0.05), ls="--", lw=1); plt.axvline(1, ls="--", lw=1); plt.axvline(-1, ls="--", lw=1)
plt.xlabel("log2FC (median tumor - normal)"); plt.ylabel("-log10 FDR (BH, Mann-Whitney)")
plt.legend(fontsize=8); plt.tight_layout()
plt.savefig(FIG/"Fig2_volcano_42568.png",dpi=200); plt.savefig(FIG/"Fig2_volcano_42568.pdf"); log(f"Fig2 volcano: {sig.sum()}/{len(de)} DEGs")

# --- Fig3 training: KM already + ROC + coef forest ---
mv=pd.read_csv(RES/"train_multivariable_cox.tsv",sep="\t")
fig,ax=plt.subplots(figsize=(7,4))
yerr=[mv.HR-mv.HR_lo95, mv.HR_hi95-mv.HR]
ax.errorbar(mv.HR, range(len(mv)), xerr=yerr, fmt="o", capsize=3)
ax.set_yticks(range(len(mv))); ax.set_yticklabels(mv.symbol, fontsize=7)
ax.axvline(1, ls="--", lw=1); ax.set_xlabel("HR (95% CI, multivariable training)"); ax.set_xscale("log")
plt.tight_layout()
plt.savefig(FIG/"Fig3_coef_forest_train.png",dpi=200); plt.savefig(FIG/"Fig3_coef_forest_train.pdf")
# ROC panel (training AUCs already computed; plot illustrative using validation-style? Use train risk file? reconstruct)
# For honest ROC curves, recompute from train risk: need T/E/risk
c425=pd.read_csv(ROOT/"metadata/GSE42568_clinical_curated.tsv",sep="\t")
tr=c425[(c425.tissue=="breast cancer")].copy(); T=tr.OS_time_days.values/365.25; E=tr.OS_event.values.astype(int)
locked=json.load(open(DER/"locked_model.json"))
e425=pd.read_csv(ROOT/"data/processed/GSE42568_expr_gene.tsv",sep="\t",index_col=0)
X=e425.loc[locked["genes"],tr.GSM].T; Z=(X-pd.Series(locked["scaling"]["means"]))/pd.Series(locked["scaling"]["sds"])
risk=Z.values@np.array(locked["coefs"])
from sklearn.metrics import roc_curve, auc
plt.figure(figsize=(6,5))
for t0 in [1,3,5]:
    y=((T<=t0)&(E==1)).astype(int); mask=((T<=t0)&(E==1))|(T>t0)
    fpr,tpr,_=roc_curve(y[mask],risk[mask]); plt.plot(fpr,tpr,label=f"{t0}y AUC={auc(fpr,tpr):.2f} (n={mask.sum()})")
plt.plot([0,1],[0,1],"--",lw=1); plt.xlabel("FPR"); plt.ylabel("TPR"); plt.legend()
plt.tight_layout()
plt.savefig(FIG/"Fig3_ROC_train.png",dpi=200); plt.savefig(FIG/"Fig3_ROC_train.pdf")

# --- Fig4 validation ROC ---
vr=pd.read_csv(RES/"validation_risk_GSE20685.tsv",sep="\t")
plt.figure(figsize=(6,5))
for t0 in [1,3,5]:
    y=((vr.OS_years<=t0)&(vr.OS_event==1)).astype(int); mask=((vr.OS_years<=t0)&(vr.OS_event==1))|(vr.OS_years>t0)
    fpr,tpr,_=roc_curve(y[mask],vr.risk[mask]); plt.plot(fpr,tpr,label=f"{t0}y AUC={auc(fpr,tpr):.2f} (n={mask.sum()})")
plt.plot([0,1],[0,1],"--",lw=1); plt.xlabel("FPR"); plt.ylabel("TPR"); plt.legend()
plt.tight_layout()
plt.savefig(FIG/"Fig4_ROC_valid.png",dpi=200); plt.savefig(FIG/"Fig4_ROC_valid.pdf")

# --- Fig5 immune pathway bars (validation) ---
comp=pd.read_csv(RES/"immune_compare_GSE20685.tsv",sep="\t")
plt.figure(figsize=(8,4))
cols=["red" if (p<0.05 and d>0) else "blue" if (p<0.05 and d<0) else "gray" for p,d in zip(comp.p_adj_BH, comp.median_diff_High_minus_Low)]
plt.barh(comp.pathway, comp.median_diff_High_minus_Low, color=cols)
plt.axvline(0,lw=1); plt.xlabel("Median z-mean score High - Low (validation)")
plt.legend(handles=[plt.Rectangle((0,0),1,1,fc="red"),plt.Rectangle((0,0),1,1,fc="blue"),plt.Rectangle((0,0),1,1,fc="gray")],labels=["FDR<0.05, High>Low","FDR<0.05, High<Low","ns"],fontsize=8); plt.tight_layout()
plt.savefig(FIG/"Fig5_pathway_valid.png",dpi=200); plt.savefig(FIG/"Fig5_pathway_valid.pdf")
# checkpoint box (validation top hits)
chk=pd.read_csv(RES/"checkpoint_expr_GSE20685.tsv",sep="\t")
top=["CD274","CTLA4","LAG3","IDO1","GZMB","PRF1"]
fig,axes=plt.subplots(2,3,figsize=(10,6),sharey=False)
for ax,g in zip(axes.ravel(),top):
    h=chk[chk.group=="High"][g].values; l=chk[chk.group=="Low"][g].values
    ax.boxplot([l,h],labels=["Low","High"],showfliers=False); ax.set_title(g,fontsize=9)
plt.tight_layout()
plt.savefig(FIG/"Fig5_checkpoints_valid.png",dpi=200); plt.savefig(FIG/"Fig5_checkpoints_valid.pdf")

# --- Fig1 flow (schematic; no in-figure titles/captions per BCRT artwork rule; see legends.md) ---
plt.figure(figsize=(9,5)); plt.axis("off")
steps=["GSE42568\n104 tumors + 17 normals\nOS 35 events\n(TRAIN)","GSE20685\n327 tumors\nOS 83 deaths\n(LOCKED VALID)","GSE45827\n130 tumors + 11 normals\n(no survival)\n(BIOLOGY)","SCAN-B GSE96058\nRNA-seq, 3273 patients\nOS 336 deaths\n(LOCKED X-PLATFORM)"]
for i,s in enumerate(steps):
    plt.text(0.09+i*0.273,0.6,s,ha="center",va="center",bbox=dict(boxstyle="round",fc="white",ec="black"),fontsize=8)
plt.tight_layout(); plt.savefig(FIG/"Fig1_flow.png",dpi=200); plt.savefig(FIG/"Fig1_flow.pdf")

# --- Tables ---
# NOTE (2026-09-24): Tab1/Tab3 old versions superseded by Tab1_cohorts_v3/Tab3_performance_v3
# (script 16, EXPORT-001); this script now writes only Tab2. Old files live in tables/archive/.
# Tab2 coefficients (already) -> copy with CIs; display rounding for publication
# (raw full precision retained in results/raw/train_multivariable_cox.tsv)
for c in ["coef", "SE", "HR", "HR_lo95", "HR_hi95"]:
    mv[c] = pd.to_numeric(mv[c], errors="coerce").round(3)
mv["p"] = pd.to_numeric(mv["p"], errors="coerce").round(4)
mv.to_csv(TAB/"Tab2_coefficients.tsv", sep="\t", index=False)
# Tab3 performance
# (Tab3 old version removed here; see note above.)
# Tab3 old version superseded by Tab3_performance_v3 (script 16, EXPORT-001) — no write here.
log("WROTE Figs 1-5 + Tab2")
logf.close()
