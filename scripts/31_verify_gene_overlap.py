"""Verify zero gene overlap: locked 14 vs PAM50-50 (OVERLAP-001).

PAM50 50-gene list (Parker et al., J Clin Oncol 2009; PMID 19204204) is
transcribed below with citation. CorePAM is derived by subsetting PAM50
(CorePAM paper, Table 2 coverage + text), so 14∩PAM50=∅ implies
14∩CorePAM=∅. Asserts 0/14 and Jaccard 0.0; fails loudly otherwise.
Outputs: results/raw/gene_overlap.json, logs/gene_overlap.log
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
DER = ROOT / "results/derived"
LOGS = ROOT / "logs"


def log(m):
    print(m, flush=True)
    with open(LOGS / "gene_overlap.log", "a") as f:
        f.write(m + "\n")


PAM50_SRC = ("metadata/pam50_genefu_Parker2009.rda "
             "(BHKLAB/genefu data/pam50.rda, Parker et al. JCO 2009 centroids; "
             "fetched 2026-09-26 from "
             "https://raw.githubusercontent.com/BHKLAB/genefu/master/data/pam50.rda)")
import rdata
_cent = rdata.read_rda(str(ROOT / "metadata/pam50_genefu_Parker2009.rda"))["pam50"]["centroids"]
PAM50 = sorted(str(x) for x in _cent.coords[_cent.dims[0]].values)
assert len(PAM50) == 50, f"PAM50 list must be 50, got {len(PAM50)}"
locked = json.load(open(DER / "locked_model.json"))["genes"]
assert len(locked) == 14
ov = sorted(set(locked) & set(PAM50))
jacc = len(ov) / len(set(locked) | set(PAM50))
log(f"locked14 ∩ PAM50-50 = {ov} ({len(ov)}/14); Jaccard={jacc:.4f}")
assert len(ov) == 0, f"OVERLAP FOUND: {ov}"
out = {"locked14": locked, "pam50_n": 50, "overlap": ov,
       "overlap_frac": "0/14", "jaccard": round(jacc, 4),
       "corepam_inference": "CorePAM subsets PAM50-50; 14∩PAM50=∅ implies 14∩CorePAM=∅",
       "pam50_source": PAM50_SRC}
json.dump(out, open(RES / "gene_overlap.json", "w"), indent=2)
log("WROTE gene_overlap.json (OVERLAP-001)")
