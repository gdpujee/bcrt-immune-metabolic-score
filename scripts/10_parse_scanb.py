#!/usr/bin/env python3
"""10_parse_scanb.py — RUN-ID: SCANB-001
Parse GSE96058 SOFT clinical (parse stage: dedup by scan-b external ID; title-tagged
replicates RETAINED here and excluded downstream in script 11 — see n_repl_tagged_retained).
Inputs: data/raw/GSE96058_family.soft.gz
Outputs: metadata/GSE96058_clinical_raw.tsv, metadata/GSE96058_clinical_curated.tsv, results/raw/scanb_summary.json
Rule: keep first GSM per scan-b ID; OS = overall survival days/event (all-cause assumed).
"""
import gzip, json, re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/"data/raw"; META = ROOT/"metadata"; RES = ROOT/"results/raw"; LOGS = ROOT/"logs"
META.mkdir(parents=True, exist_ok=True); RES.mkdir(parents=True, exist_ok=True)
logf = open(LOGS/"scanb_parse.log","w")
def log(m): print(m); logf.write(m+"\n"); logf.flush()

# parse SOFT sample blocks: ^SAMPLE = GSM... lines with !Sample_* fields
samples = []; cur = None
with gzip.open(RAW/"GSE96058_family.soft.gz","rt",errors="replace") as f:
    for line in f:
        if line.startswith("^SAMPLE"):
            if cur: samples.append(cur)
            cur = {"GSM": line.strip().split("=")[1].strip(), "chars": []}
        elif line.startswith("!Sample_geo_accession"):
            pass
        elif line.startswith("!Sample_title"):
            if cur is not None: cur["title"] = line.split("=",1)[1].strip()
        elif line.startswith("!Sample_source_name_ch1"):
            if cur is not None: cur["source"] = line.split("=",1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1"):
            if cur is not None: cur["chars"].append(line.split("=",1)[1].strip())
        elif line.startswith("!Sample_platform_id"):
            if cur is not None: cur["platform"] = line.split("=",1)[1].strip()
if cur: samples.append(cur)
log(f"SOFT samples: {len(samples)}")
# pivot chars per sample
rows=[]
for s in samples:
    d = {"GSM": s["GSM"], "title": s.get("title",""), "source": s.get("source",""), "platform": s.get("platform","")}
    for c in s["chars"]:
        if ":" in c:
            k,v = c.split(":",1)
            k = re.sub(r"[^a-z0-9]+","_",k.strip().lower()).strip("_")
            if k and k not in d: d[k]=v.strip()
    rows.append(d)
clin = pd.DataFrame(rows)
clin.to_csv(META/"GSE96058_clinical_raw.tsv", sep="\t", index=False)
log(f"Columns: {clin.columns.tolist()}")
# OS + replicate audit
clin["OS_days"] = pd.to_numeric(clin.get("overall_survival_days"), errors="coerce")
clin["OS_event"] = pd.to_numeric(clin.get("overall_survival_event"), errors="coerce")
clin["OS_years"] = clin["OS_days"]/365.25
# scan-b ID for dedup (col name normalized)
scan_col = [c for c in clin.columns if "scan" in c]
log(f"scan cols: {scan_col}")
key = scan_col[0] if scan_col else "GSM"
dups = clin.duplicated(subset=[key], keep=False).sum()
log(f"Duplicate {key} rows: {dups}; 136 title-tagged replicate rows tracked separately (excluded downstream in script 11)")
cur2 = clin.sort_values("GSM").drop_duplicates(subset=[key], keep="first").copy()
excluded = len(clin)-len(cur2)
n_repl_tagged = int(clin.title.str.contains("repl", na=False).sum())
log(f"Parsed rows (unique by {key}): {len(cur2)} (excluded {excluded} ID-duplicates; {n_repl_tagged} title-tagged replicates RETAINED here, excluded downstream in script 11)")
log(f"OS events in parsed: {(cur2.OS_event==1).sum()} deaths / {(cur2.OS_event==0).sum()} censored; median years {cur2.OS_years.median():.2f}")
# save curated + manifest rows appended later after expression check
cur2.to_csv(META/"GSE96058_clinical_curated.tsv", sep="\t", index=False)
summary = {"n_soft": len(clin), "n_parsed_rows_unique_by_id": len(cur2), "n_id_duplicates_excluded": int(excluded),
 "n_repl_tagged_retained": n_repl_tagged, "n_norepl_titles": int(len(cur2) - n_repl_tagged),
 "n_deaths": int((cur2.OS_event==1).sum()), "n_censored": int((cur2.OS_event==0).sum()),
 "median_years": float(cur2.OS_years.median()), "key": key,
 "pmids": ["32913985","32926574","33937624","35304506"]}
json.dump(summary, open(RES/"scanb_summary.json","w"), indent=2)
log("WROTE GSE96058 curated clinical + scanb_summary.json")
logf.close()
