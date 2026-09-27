"""Verify the one-sample-per-patient assumption behind the METABRIC covariate join.

`scripts/33_metabric_adjusted.py` merges SAMPLE-level clinical attributes onto a
patient-level frame keyed by `patientId`.  That is only valid if brca_metabric
contributes one tumour sample per patient to the analysis.  Rather than assume it,
this probe asks the live cBioPortal API how many distinct sampleIds share a
patientId, and records the answer next to the adjusted-Cox results so the
assumption is auditable (external review P1-12, missing-value transparency).

Outputs: results/raw/metabric_sample_multiplicity.json, logs/sample_multiplicity.log
"""
import json
import time
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
LOGS = ROOT / "logs"
logf = open(LOGS / "sample_multiplicity.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


BASE = "https://www.cbioportal.org/api"
ST = "brca_metabric"


def get(url, tries=5):
    last = None
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (a + 1))
    raise RuntimeError(f"GET failed: {url}: {last}")


# ER_STATUS is present for the largest number of samples, so it gives the widest
# view of the sample->patient mapping.
recs, page = [], 0
while True:
    d = get(f"{BASE}/studies/{ST}/clinical-data?clinicalDataType=SAMPLE"
            f"&attributeId=ER_STATUS&projection=SUMMARY&pageSize=5000&pageNumber={page}")
    if not d:
        break
    recs += d
    page += 1
    if page > 12:
        raise RuntimeError("pagination guard tripped")

pairs = [(x["sampleId"], x["patientId"]) for x in recs]
by_patient = Counter(p for _, p in pairs)
multi = {p: c for p, c in by_patient.items() if c > 1}

out = {
    "study": ST,
    "attribute_probed": "ER_STATUS",
    "sample_records": len(recs),
    "distinct_sampleIds": len({s for s, _ in pairs}),
    "distinct_patientIds": len(by_patient),
    "patients_with_multiple_samples": len(multi),
    "max_samples_per_patient": max(by_patient.values()) if by_patient else 0,
    "examples": dict(list(sorted(multi.items()))[:10]),
    "conclusion": ("one tumour sample per patient; the patientId-keyed SAMPLE join in "
                   "scripts/33_metabric_adjusted.py is unambiguous"
                   if not multi else
                   "MULTIPLE SAMPLES PER PATIENT — the patientId-keyed join silently "
                   "keeps one arbitrary sample and must be revisited"),
}
json.dump(out, open(RES / "metabric_sample_multiplicity.json", "w"), indent=2)
for k, v in out.items():
    if k != "examples":
        log(f"{k}: {v}")
log("WROTE results/raw/metabric_sample_multiplicity.json (SAMPLEMULT-001)")
logf.close()
