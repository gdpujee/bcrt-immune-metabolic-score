"""Master sample-manifest builder (MANIFEST-001).

Deterministically rebuilds the full 4012-row SAMPLE_MANIFEST.tsv:
  - GEO cohorts (GSE42568/45827/20685, 603 rows) from clinical_curated TSVs
    (file order = header order; patient_ID = GSM; exclusion empty).
  - SCAN-B (3409 rows) from GSE96058_clinical_curated.tsv (GSM-sorted):
    patient_ID = scan_b_external_id (+ "_repl" iff title-tagged replicate);
    exclusion = "technical replicate" iff repl-tagged else "".
Contract: output byte-identical to committed SAMPLE_MANIFEST.tsv (git HEAD),
else exit nonzero. This closes the S5 rerun gap (script 02 alone rebuilds
only the 603 GEO rows).
"""
import subprocess
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata"
LOGS = ROOT / "logs"
SEED = 42


def log(m):
    print(m, flush=True)
    with open(LOGS / "manifest.log", "a") as f:
        f.write(m + "\n")


rows = []
for gse in ["GSE42568", "GSE45827", "GSE20685"]:
    c = pd.read_csv(META / f"{gse}_clinical_curated.tsv", sep="\t")
    for _, r in c.iterrows():
        rows.append({"repository_sample_ID": r["GSM"], "cohort": gse,
                     "patient_ID": r["GSM"], "platform": r.get("platform", ""),
                     "title": r.get("title", ""), "source": r.get("source", ""),
                     "exclusion": ""})
    log(f"{gse}: {len(c)} rows")
cu = pd.read_csv(META / "GSE96058_clinical_curated.tsv", sep="\t")
for _, r in cu.iterrows():
    repl = bool(pd.notna(r["title"]) and "repl" in str(r["title"]))
    rows.append({"repository_sample_ID": r["GSM"], "cohort": "GSE96058",
                 "patient_ID": str(r["scan_b_external_id"]) + ("_repl" if repl else ""),
                 "platform": r.get("platform", ""), "title": r.get("title", ""),
                 "source": r.get("source", ""),
                 "exclusion": "technical replicate" if repl else ""})
log(f"GSE96058: {len(cu)} rows ({int(cu.title.str.contains('repl', na=False).sum())} repl-tagged)")
man = pd.DataFrame(rows, columns=["repository_sample_ID", "cohort", "patient_ID",
                                  "platform", "title", "source", "exclusion"])
text = man.to_csv(sep="\t", index=False)
ref = subprocess.run(["git", "show", "HEAD:SAMPLE_MANIFEST.tsv"],
                     capture_output=True, text=True, cwd=ROOT).stdout
if text != ref:
    open("/tmp/manifest_new.tsv", "w").write(text)
    raise SystemExit("MISMATCH vs committed SAMPLE_MANIFEST.tsv (see /tmp/manifest_new.tsv)")
open(ROOT / "SAMPLE_MANIFEST.tsv", "w").write(text)
dups = len(man) - man.repository_sample_ID.nunique()
log(f"MANIFEST-001: {len(man)} rows, GSM dups={dups}; exported + verified vs HEAD; PASS")
