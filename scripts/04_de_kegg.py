#!/usr/bin/env python3
"""04_de_kegg.py — RUN-ID: DE-001
DE tumor-vs-normal (Mann-Whitney + BH-FDR) + KEGG immune-metabolic candidate pool.
Inputs: data/processed/*_expr_gene.tsv, metadata/*_curated.tsv
Outputs: results/raw/DE_*.tsv, results/raw/KEGG_*.{json,gmt}, results/raw/candidate_pool.tsv, logs/kegg_fetch.log
KEGG: FROZEN committed snapshot by default (reproducible reruns); --refresh refetches
from rest.kegg.jp with real date + content hash, failing hard on any fetch error.
Pathway sets used in the implemented candidate-pool procedure; the frozen plan does
not establish the pathway identities or relaxed fold-change threshold as prespecified.
"""
import pandas as pd, numpy as np, json, requests, time, argparse, hashlib, datetime
from pathlib import Path
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT/"data/processed"; META = ROOT/"metadata"; RES = ROOT/"results/raw"; LOGS = ROOT/"logs"
# KEGG sets used in the implemented procedure; the frozen plan does not establish
# that these pathway identities were prespecified.
METABOLIC = {"hsa00010":"Glycolysis_Gluconeogenesis","hsa00071":"Fatty_acid_degradation","hsa00061":"Fatty_acid_biosynthesis","hsa00190":"Oxidative_phosphorylation","hsa00030":"Pentose_phosphate","hsa00480":"Glutathione_metabolism"}
IMMUNE = {"hsa04060":"Cytokine_cytokine_receptor","hsa04062":"Chemokine_signaling","hsa04612":"Antigen_processing_presentation","hsa04650":"NK_cytotoxicity","hsa04660":"T_cell_receptor"}
ALL_KEGG = {**METABOLIC, **IMMUNE}

def kegg_genes(hsa):
    # returns set of human gene symbols for pathway via KEGG link + list
    # Step 1: get gene list: http://rest.kegg.jp/link/hsa/<hsa> gives hsa:ENTREZ -> pathway; then convert ENTREZ to symbol via ... use kegg get for each? Faster: use link + precomputed? Instead use: http://rest.kegg.jp/get/<hsa> lists genes with symbols in ORTHOLOGY/GENES section.
    # Robust: fetch get/hsaXXXXX, parse lines starting with gene ENTREZ + symbol.
    r = requests.get(f"https://rest.kegg.jp/get/{hsa}", timeout=60)
    r.raise_for_status()
    genes = set()
    for line in r.text.splitlines():
        # GENE section lines look like: "  10327  AKR1B1; aldo-keto reductase family 1 member B1 [KO:...]"
        # We'll parse tokens after stripping; symbol is second token split by ';'
        s = line.strip()
        if not s: continue
        # heuristic: lines with ';' and starting with digits
        if ";" in s and s[0].isdigit():
            parts = s.split()
            if len(parts) >= 2:
                sym = parts[1].strip(";").strip()
                # filter KO/EC tokens
                if sym and sym[0].isupper() and len(sym) <= 15 and "[" not in sym:
                    genes.add(sym.split(";")[0])
    return genes

# fetch KEGG with provenance (frozen by default; --refresh refetches with fail-hard)
# NOTE: argparse runs BEFORE the log is opened, so --help never truncates the log.
ap = argparse.ArgumentParser()
ap.add_argument("--refresh", action="store_true",
                help="refetch KEGG from rest.kegg.jp (records real date + hash; fails hard on error)")
args = ap.parse_args()
logf = open(LOGS/"kegg_fetch.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()
kegg_release = ""
if not args.refresh:
    frozen = json.load(open(RES / "KEGG_genesets.json"))
    gmt = {k: set(v) for k, v in frozen["genesets"].items()}
    sizes = {k: len(v) for k, v in gmt.items()}
    log(f"KEGG frozen snapshot (no network): {sum(sizes.values())} memberships across {len(gmt)} sets")
else:
    try:
        r = requests.get("https://rest.kegg.jp/info/kegg", timeout=30)
        r.raise_for_status()
        kegg_release = r.text.splitlines()[0] if r.ok else "unknown"
    except Exception as e:
        raise SystemExit(f"REFRESH ABORTED: KEGG info fetch failed: {e}")
    log(f"KEGG info: {kegg_release} (access {datetime.date.today()})")
    gmt = {}; sizes = {}
    for hsa, name in ALL_KEGG.items():
        try:
            g = kegg_genes(hsa)
            time.sleep(0.5)
        except Exception as e:
            raise SystemExit(f"REFRESH ABORTED: {hsa} {name} fetch failed: {e}")
        gmt[name] = sorted(g)
        sizes[name] = len(g)
        log(f"{hsa} {name}: {len(g)} symbols")
    blob = json.dumps({"release": kegg_release, "access": str(datetime.date.today()),
                       "pathways": ALL_KEGG, "sizes": sizes,
                       "genesets": {k: sorted(v) for k, v in gmt.items()}}, indent=2)
    log(f"KEGG refresh sha256: {hashlib.sha256(blob.encode()).hexdigest()}")
    with open(RES / "KEGG_genesets.json", "w") as f:
        f.write(blob)
    gmt = {k: set(v) for k, v in gmt.items()}
# union pools
met_genes = set().union(*[set(gmt[k]) for k in METABOLIC.values()])
imm_genes = set().union(*[set(gmt[k]) for k in IMMUNE.values()])
log(f"Union metabolic={len(met_genes)}, immune={len(imm_genes)}, combined={len(met_genes|imm_genes)}")

def de_tumor_normal(gse, tumor_label, normal_label, tumor_col="tissue"):
    clin = pd.read_csv(META/f"{gse}_clinical_curated.tsv", sep="\t")
    expr = pd.read_csv(PROC/f"{gse}_expr_gene.tsv", sep="\t", index_col=0)
    # map GSM->group
    gsm2group = dict(zip(clin["GSM"], clin[tumor_col]))
    tumors = [c for c in expr.columns if gsm2group.get(c)==tumor_label]
    normals = [c for c in expr.columns if gsm2group.get(c)==normal_label]
    log(f"{gse}: tumors={len(tumors)}, normals={len(normals)}")
    rows = []
    for g in expr.index:
        t = expr.loc[g, tumors].values.astype(float)
        nn = expr.loc[g, normals].values.astype(float)
        t = t[~np.isnan(t)]; nn = nn[~np.isnan(nn)]
        log2fc = float(np.median(t) - np.median(nn))  # log2 scale -> median diff
        try:
            _, p = mannwhitneyu(t, nn, alternative="two-sided")
            if np.isnan(p): p = 1.0
        except Exception:
            p = 1.0
        rows.append((g, log2fc, float(p), float(np.median(t)), float(np.median(nn))))
    de = pd.DataFrame(rows, columns=["symbol","log2FC_median_diff","p_raw","median_tumor","median_normal"])
    _, de["p_adj_BH"], _, _ = multipletests(de["p_raw"], method="fdr_bh")
    de = de.sort_values("p_adj_BH")
    de.to_csv(RES/f"DE_{gse}_tumor_vs_normal.tsv", sep="\t", index=False)
    sig = de[(de.p_adj_BH<0.05)&(de.log2FC_median_diff.abs()>1)]
    log(f"{gse}: DEGs FDR<0.05 & |log2FC|>1: {len(sig)} / {len(de)}")
    return de

de425 = de_tumor_normal("GSE42568", "breast cancer", "normal breast", "tissue")
# GSE45827: diagnosis Breast cancer vs None (normal); curated col 'diagnosis' has 'None (normal)'? check
c458 = pd.read_csv(META/"GSE45827_clinical_curated.tsv", sep="\t")
print("GSE45827 diagnosis values:", c458["diagnosis"].fillna("NA").unique()[:10])
# normalize: map 'None (normal)' variants
c458["tissue2"] = c458["diagnosis"].map({"Breast cancer":"tumor", "None (normal)":"normal"})
# cell lines have other diagnosis? inspect
print(c458[["diagnosis","tumor_subtype","cell_line"]].head(10).to_string())
# For DE use only bio samples (exclude cell lines flagged)
c458_bio = c458[~c458["is_cell_line"]].copy()
c458_bio.to_csv(META/"GSE45827_bio_tmp.tsv", sep="\t", index=False)
# manual DE for 45827 bio
expr458 = pd.read_csv(PROC/"GSE45827_expr_gene.tsv", sep="\t", index_col=0)
gsm2g = dict(zip(c458["GSM"], c458["tissue2"]))
tumors = [c for c in expr458.columns if gsm2g.get(c)=="tumor"]
normals = [c for c in expr458.columns if gsm2g.get(c)=="normal"]
log(f"GSE45827 bio: tumors={len(tumors)}, normals={len(normals)}")
rows=[]
for g in expr458.index:
    t = expr458.loc[g, tumors].values.astype(float); nn = expr458.loc[g, normals].values.astype(float)
    t=t[~np.isnan(t)]; nn=nn[~np.isnan(nn)]
    fc=float(np.median(t)-np.median(nn))
    try:
        _, p = mannwhitneyu(t, nn, alternative="two-sided")
        if np.isnan(p): p = 1.0
    except: p = 1.0
    rows.append((g,fc,float(p)))
de458=pd.DataFrame(rows,columns=["symbol","log2FC_median_diff","p_raw"])
_,de458["p_adj_BH"],_,_=multipletests(de458["p_raw"],method="fdr_bh")
de458=de458.sort_values("p_adj_BH")
de458.to_csv(RES/"DE_GSE45827_tumor_vs_normal.tsv",sep="\t",index=False)
log(f"GSE45827: DEGs FDR<0.05 & |log2FC|>1: {((de458.p_adj_BH<0.05)&(de458.log2FC_median_diff.abs()>1)).sum()} / {len(de458)}")

# candidate pool = DEGs(42568) ∩ KEGG_union ∩ training genes
sig425 = set(de425[(de425.p_adj_BH<0.05)&(de425.log2FC_median_diff.abs()>0.5)]["symbol"])  # implemented relaxed |log2FC|>0.5 pool threshold; prespecification is not established; strict >1 reported separately
kegg_union = met_genes|imm_genes
pool = sig425 & kegg_union
# annotate metabolic vs immune
pool_rows=[]
for g in sorted(pool):
    pool_rows.append({"symbol":g, "in_metabolic":int(g in met_genes), "in_immune":int(g in imm_genes),
                      "log2FC_42568":float(de425.set_index("symbol").loc[g,"log2FC_median_diff"]),
                      "FDR_42568":float(de425.set_index("symbol").loc[g,"p_adj_BH"])})
pool_df=pd.DataFrame(pool_rows)
pool_df.to_csv(RES/"candidate_pool.tsv",sep="\t",index=False)
log(f"Candidate pool (DEG∩KEGG): {len(pool_df)} genes (met-only={(pool_df.in_metabolic==1)&(pool_df.in_immune==0)}.sum(), imm-only={(pool_df.in_metabolic==0)&(pool_df.in_immune==1)}.sum(), both={(pool_df.in_metabolic==1)&(pool_df.in_immune==1)}.sum())")
logf.close()
