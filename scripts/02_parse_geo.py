#!/usr/bin/env python3
"""02_parse_geo.py — RUN-ID: PARSE-001
Parse GEO series matrices + GPL570 annot to gene-level matrices + clinical + SAMPLE_MANIFEST.
Inputs: data/raw/GSE*_series_matrix.txt.gz, data/raw/GPL570.annot.gz
Outputs: data/processed/*_expr_gene.tsv, metadata/*_clinical.tsv, SAMPLE_MANIFEST.tsv, results/raw/parse_summary.json, logs/parse.log
Seed: 42 (deterministic max-mean probe choice). No outcome use.
"""
import gzip, json, hashlib, re, sys
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/"data/raw"; PROC = ROOT/"data/processed"; META = ROOT/"metadata"; RES = ROOT/"results/raw"; LOGS = ROOT/"logs"
for d in [PROC, META, RES, LOGS]: d.mkdir(parents=True, exist_ok=True)
logf = open(LOGS/"parse.log","w")
def log(m):
    print(m); logf.write(m+"\n"); logf.flush()

SERIES = {
 "GSE42568": {"file": RAW/"GSE42568_series_matrix.txt.gz"},
 "GSE45827": {"file": RAW/"GSE45827_series_matrix.txt.gz"},
 "GSE20685": {"file": RAW/"GSE20685_series_matrix.txt.gz"},
}

def sha256(p):
    import hashlib
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda: f.read(1<<20), b""): h.update(b)
    return h.hexdigest()

# --- GPL570 annot: probe -> symbol ---
log("Parsing GPL570.annot.gz ...")
probe2sym = {}
annot_date = "unknown"
with gzip.open(RAW/"GPL570.annot.gz","rt",errors="replace") as f:
    header = None
    for line in f:
        if line.startswith("!Annotation_date"):
            annot_date = line.strip().split("=",1)[1].strip()
        if line.startswith("#") or line.startswith("^") or line.startswith("!"):
            continue
        if header is None:
            header = line.rstrip("\n").split("\t")
            # find ID and Gene symbol columns
            try:
                id_idx = header.index("ID")
                sym_idx = header.index("Gene symbol")
            except ValueError:
                log(f"ANNOT HEADER: {header[:15]}")
                raise
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= max(id_idx, sym_idx): continue
        pid = parts[id_idx].strip().strip('"')
        sym = parts[sym_idx].strip().strip('"')
        if pid and sym and sym not in ("", "---", "NA"):
            # strip whitespace; keep first symbol if multiple (///)
            sym = sym.split("///")[0].strip()
            if sym:
                probe2sym[pid] = sym
log(f"GPL570 annot date={annot_date}, mapped probes={len(probe2sym)}")

def parse_series(gse, fpath):
    log(f"\n=== {gse} ===")
    with gzip.open(fpath,"rt",errors="replace") as f:
        meta_lines = []
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            meta_lines.append(line.rstrip("\n"))
    # sample ids in column order
    gsm_line = [l for l in meta_lines if l.startswith("!Sample_geo_accession")]
    assert len(gsm_line)==1, f"{gse}: expected 1 GSM line, got {len(gsm_line)}"
    gsms = [x.strip().strip('"') for x in gsm_line[0].split("\t")[1:]]
    n = len(gsms)
    log(f"{gse}: n_samples={n}")
    # collect all !Sample_* fields
    fields = {}
    for l in meta_lines:
        if l.startswith("!Sample_"):
            key = l.split("\t",1)[0]
            vals = [x.strip().strip('"') for x in l.split("\t")[1:]]
            fields[key] = vals
    # characteristics: RAGGED across samples (GEO artifact) -> per-sample key:value union
    char_rows = [l for l in meta_lines if l.startswith("!Sample_characteristics_ch1")]
    clin = pd.DataFrame({"GSM": gsms})
    clin["GSE"] = gse
    if "!Sample_title" in fields: clin["title"] = fields["!Sample_title"][:n]
    if "!Sample_source_name_ch1" in fields: clin["source"] = fields["!Sample_source_name_ch1"][:n]
    if "!Sample_platform_id" in fields: clin["platform"] = fields["!Sample_platform_id"][:n]
    # build per-sample dicts
    per_sample = [dict() for _ in range(n)]
    for cr in char_rows:
        vals = [x.strip().strip('"') for x in cr.split("\t")[1:]]
        for j in range(min(n, len(vals))):
            v = vals[j]
            if not v or ":" not in v:
                continue
            k, vv = v.split(":", 1)
            k = re.sub(r"[^a-z0-9]+", "_", k.strip().lower()).strip("_")
            vv = vv.strip()
            if k and k not in per_sample[j]:
                per_sample[j][k] = vv
            elif k and per_sample[j].get(k) in ("", None):
                per_sample[j][k] = vv
    # union keys -> columns
    all_keys = sorted({k for d in per_sample for k in d.keys()})
    log(f"{gse}: characteristics keys ({len(all_keys)}): {all_keys[:20]}")
    for k in all_keys:
        clin[k] = [d.get(k, "") for d in per_sample]
    # expression block
    log(f"{gse}: reading expression block ...")
    expr = pd.read_csv(fpath, sep="\t", compression="gzip", comment="!",
                       index_col=0, low_memory=False)
    # expr index = probe ID (may include quotes); clean
    expr.index = [str(x).strip().strip('"') for x in expr.index]
    expr.columns = [str(x).strip().strip('"') for x in expr.columns]
    # columns should be GSMs; align
    # series matrix columns are GSMs in order; rename positionally to be safe
    if list(expr.columns) != gsms:
        log(f"{gse}: WARNING expr columns != GSM order; renaming positionally ({expr.shape})")
        if expr.shape[1] == n:
            expr.columns = gsms
        else:
            raise ValueError(f"{gse}: expr shape {expr.shape} vs n {n}")
    log(f"{gse}: probe-level shape={expr.shape}, value range [{np.nanmin(expr.values):.2f}, {np.nanmax(expr.values):.2f}], median {np.nanmedian(expr.values):.2f}")
    return clin, expr

all_clin = []
manifest_rows = []
summary = {"annot_date": annot_date, "n_mapped_probes": len(probe2sym), "cohorts": {}}

for gse, info in SERIES.items():
    clin, expr = parse_series(gse, info["file"])
    # save raw clinical
    clin.to_csv(META/f"{gse}_clinical_raw.tsv", sep="\t", index=False)
    # probe->gene: keep mapped probes only, max-mean probe per gene
    expr_mapped = expr[expr.index.isin(probe2sym)]
    log(f"{gse}: mapped probes {expr_mapped.shape[0]}/{expr.shape[0]}")
    # mean per probe for max-mean selection
    probe_means = expr_mapped.mean(axis=1, skipna=True)
    tmp = pd.DataFrame({"probe": expr_mapped.index, "symbol": [probe2sym[p] for p in expr_mapped.index], "mean": probe_means.values})
    best = tmp.sort_values("mean", ascending=False).drop_duplicates("symbol", keep="first").set_index("symbol")["probe"]
    log(f"{gse}: genes after max-mean: {len(best)}")
    gene_expr = expr_mapped.loc[best.values].copy()
    gene_expr.index = best.index  # symbol
    gene_expr.index.name = "symbol"
    gene_expr.to_csv(PROC/f"{gse}_expr_gene.tsv", sep="\t")
    # manifest rows
    for _, r in clin.iterrows():
        manifest_rows.append({"repository_sample_ID": r["GSM"], "cohort": gse,
                              "patient_ID": r["GSM"],  # one sample per patient by design; duplicates audited below
                              "platform": r.get("platform",""),
                              "title": r.get("title",""), "source": r.get("source","")})
    summary["cohorts"][gse] = {"n": len(clin), "n_probes": int(expr.shape[0]),
                                "n_genes": int(len(best)),
                                "sha256": sha256(info["file"]),
                                "columns": list(clin.columns)}
    all_clin.append(clin)

# GSM overlap audit
all_gsm = [r["repository_sample_ID"] for r in manifest_rows]
overlap = len(all_gsm) - len(set(all_gsm))
log(f"\nGSM overlap across cohorts: {overlap} (must be 0)")
# gene intersection
import functools
gene_sets = [set(pd.read_csv(PROC/f"{g}_expr_gene.tsv", sep="\t", nrows=0).columns[1:]) if False else None for g in SERIES]
# simpler: read index
gene_lists = {}
for g in SERIES:
    df = pd.read_csv(PROC/f"{g}_expr_gene.tsv", sep="\t", usecols=["symbol"])
    gene_lists[g] = set(df["symbol"])
inter = set.intersection(*gene_lists.values())
log(f"Gene intersection across 3 cohorts: {len(inter)} (per-cohort: {[(g,len(s)) for g,s in gene_lists.items()]})")
summary["gene_intersection"] = len(inter)
summary["gsm_overlap"] = overlap
summary["per_cohort_genes"] = {g: len(s) for g,s in gene_lists.items()}

man = pd.DataFrame(manifest_rows)
man.to_csv(ROOT/"SAMPLE_MANIFEST.tsv", sep="\t", index=False)
with open(RES/"parse_summary.json","w") as f: json.dump(summary, f, indent=2)
log(f"\nWROTE: SAMPLE_MANIFEST.tsv ({len(man)} rows), gene matrices, parse_summary.json")
log("DONE PARSE-001")
logf.close()
