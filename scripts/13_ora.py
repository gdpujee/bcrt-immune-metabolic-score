#!/usr/bin/env python3
"""13_ora.py — RUN-ID: ORA-001
Over-representation of 14 locked genes across 11 KEGG sets (hypergeometric, BH over 11).
Background = 21,755 training genes (documented choice; sensitivity with 14,522 intersection in log).
Inputs: results/derived/locked_model.json, results/raw/KEGG_genesets.json, data/processed/GSE42568_expr_gene.tsv (background)
Outputs: results/raw/ORA_14g.tsv, logs/ora.log
"""
import pandas as pd, json
from pathlib import Path
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT/"results/raw"; DER = ROOT/"results/derived"; LOGS = ROOT/"logs"
logf = open(LOGS/"ora.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()

locked = json.load(open(DER/"locked_model.json"))["genes"]
kegg = json.load(open(RES/"KEGG_genesets.json"))["genesets"]
bg = set(pd.read_csv(ROOT/"data/processed/GSE42568_expr_gene.tsv", sep="\t", usecols=["symbol"])["symbol"])
inter = set(pd.read_csv(ROOT/"data/processed/GSE45827_expr_gene.tsv", sep="\t", usecols=["symbol"])["symbol"])
S = set(locked)
log(f"Foreground {len(S)}, background {len(bg)} (intersection-sensitivity {len(inter)})")
rows=[]
for name, members in kegg.items():
    M = set(m for m in members if m in bg)
    k = len(S & M)
    # hypergeometric: N=|bg|, K=|M|, n=|S|, x=k (upper tail)
    p = float(hypergeom.sf(k-1, len(bg), len(M), len(S))) if M else 1.0
    rows.append({"pathway":name,"set_size_in_bg":len(M),"overlap":k,"overlap_genes":";".join(sorted(S&M)),"p_ORA":p})
ora = pd.DataFrame(rows)
_, ora["p_adj_BH"], _, _ = multipletests(ora.p_ORA, method="fdr_bh")
ora = ora.sort_values("p_adj_BH")
ora.to_csv(RES/"ORA_14g.tsv", sep="\t", index=False)
log(ora.to_string(index=False))
# sensitivity: background = intersection
rows2=[]
for name, members in kegg.items():
    M = set(m for m in members if m in inter)
    k = len(S & M)
    p = float(hypergeom.sf(k-1, len(inter), len(M), len(S))) if M else 1.0
    rows2.append((name, len(M), k, p))
log("Intersection-background sensitivity (set_size, overlap, p):")
for r in rows2: log(str(r))
log("WROTE ORA_14g.tsv (ORA-001)")
logf.close()
