"""Recompute the locked-cutoff high/low splits from raw risk files (SPLIT-001).

The transportability table reports how badly the frozen cutoff transports: 48/279
(GSE20685), 156/3117 (SCAN-B) and 1909/71 (METABRIC) patients fall above/below it.
Those three counts were previously typed into the table generator by hand.  Here
they are recomputed from the locked cutoff in results/raw/train_summary.json and
the per-cohort risk files, so the table carries derived numbers instead of magic
constants and any drift in the risk files shows up immediately.

Outputs: results/raw/cutoff_splits.json, logs/cutoff_splits.log
"""
import json
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results/raw"
META = ROOT / "metadata"
logf = open(ROOT / "logs/cutoff_splits.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


cut = json.load(open(RAW / "train_summary.json"))["cutoff"]
log(f"locked cutoff (train_summary.json): {cut!r}")

out = {"locked_cutoff": cut, "rule": "high = risk > cutoff, low = risk <= cutoff"}


def split(label, risk, key):
    hi = int((risk > cut).sum())
    lo = int((risk <= cut).sum())
    out[key] = {"n": int(len(risk)), "high": hi, "low": lo,
                "high_pct": round(100 * hi / len(risk), 2),
                "sd": round(float(risk.std()), 4),
                "text": f"{hi}/{lo}"}
    log(f"{label}: n={len(risk)} high={hi} low={lo} ({out[key]['high_pct']}%) "
        f"SD={out[key]['sd']:.4f} -> '{hi}/{lo}'")


# GSE20685
g = pd.read_csv(META / "GSE20685_clinical_curated.tsv", sep="\t")
gr = pd.read_csv(RAW / "validation_risk_GSE20685.tsv", sep="\t")
g = (g.merge(gr[["GSM", "risk"]], on="GSM")
      .dropna(subset=["follow_up_duration_years", "event_death", "risk"]))
split("GSE20685", g.risk, "GSE20685")

# SCAN-B
s = pd.read_csv(RAW / "rnaseq_risk_GSE96058.tsv", sep="\t")
sc = pd.read_csv(META / "GSE96058_clinical_raw.tsv", sep="\t")
sc = sc[~sc.title.astype(str).str.contains("repl", na=False)]
sc = sc.merge(s[["sample", "risk"]], left_on="title", right_on="sample")
sc["time"] = pd.to_numeric(sc.overall_survival_days, errors="coerce") / 365.25
sc["event"] = pd.to_numeric(sc.overall_survival_event, errors="coerce")
sc = sc.dropna(subset=["time", "event", "risk"])
split("SCAN-B", sc.risk, "SCANB")

# METABRIC
r = pd.read_csv(RAW / "metabric_risk.tsv", sep="\t")
clin = pd.read_csv(META / "METABRIC_clinical_dl.tsv", sep="\t").drop_duplicates("patientId")
mm = r.merge(clin, on="patientId")
mm = mm[pd.to_numeric(mm.OS_MONTHS, errors="coerce").notna()]
split("METABRIC", mm.risk, "METABRIC")

# The values the manuscript already states; a change here means the risk files moved
# and every document quoting the split has to be revisited.
EXPECTED = {"GSE20685": (48, 279), "SCANB": (156, 3117), "METABRIC": (1909, 71)}
for k, (hi, lo) in EXPECTED.items():
    got = (out[k]["high"], out[k]["low"])
    if got != (hi, lo):
        raise SystemExit(f"{k}: recomputed split {got} != documented {hi,lo} — the "
                         f"manuscript, Table 3 and the supplement all quote this split")
log("all three splits reproduce the documented values")

json.dump(out, open(RAW / "cutoff_splits.json", "w"), indent=2)
log("WROTE results/raw/cutoff_splits.json (SPLIT-001)")
logf.close()
