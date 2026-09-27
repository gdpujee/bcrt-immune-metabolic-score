#!/usr/bin/env python3
"""Live audit of the published release (PUBLICAUDIT-001).

Answers the one question no offline gate can: **does the thing the manuscript
points at actually exist and actually match?**  The manuscript states that the
analysis code and derived tables are available at a URL and archived under a
concept DOI.  That statement can be false in ways that are invisible from inside
the repository — the repository was live and returned 200 while its contents were
a snapshot predating the METABRIC validation, the adjusted Cox models and the
whole packaging round, and its archived metadata listed a GitHub login as the
creator, carried no description, no keywords, and declared CC-BY-4.0 for a
repository that had no licence file.

This script is deliberately NOT part of `40_submission_gate.py`.  The gate must
run offline and deterministically; this one needs the network and its result
depends on a third party being up.  Run it before a release and after publishing.

Checks
  1. the cited repository URL resolves
  2. the published tree contains the scripts the manuscript's chain needs
  3. the published tree matches the export manifest — files missing upstream, and
     files upstream that the export no longer publishes
  4. the concept DOI resolves through doi.org
  5. the Zenodo record names both authors (not a GitHub login) and carries a
     description, keywords and an MIT licence
  6. the repository carries a LICENSE file

Exit status: 0 all clear, 1 a real mismatch, 2 something could not be checked.

Outputs: logs/public_release_audit.log
"""
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "logs/public_release_audit.log"

OWNER_REPO = "gdpu11/breast-cancer-immune-metabolic-score"
REPO_URL = f"https://github.com/{OWNER_REPO}"
CONCEPT_DOI = "10.5281/zenodo.22961119"
CONCEPT_RECORD = CONCEPT_DOI.rsplit(".", 1)[1]
MANIFEST = ROOT / "dist/public_release.manifest.json"

# Scripts the paper's reproduction chain cannot do without, checked by presence
# rather than by count so a truncated upload is caught rather than averaged away.
REQUIRED_UPSTREAM = ("scripts/02_parse_geo.py", "scripts/05_prognosis_train.py",
                     "scripts/11_rnaseq_validate.py", "scripts/27_metabric_validate.py",
                     "scripts/33_metabric_adjusted.py", "scripts/40_submission_gate.py",
                     "requirements.txt", "LICENSE", "CITATION.cff", ".zenodo.json",
                     "README.md", "tables/Tab3_performance_v3.tsv")

logf = open(LOG, "w")
fails, unchecked = [], []


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


def check(name, ok, detail=""):
    log(("PASS  " if ok else "FAIL  ") + name + (f" — {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def note(name, detail):
    """A check that could not run.  Never counted as a pass."""
    log(f"SKIP  {name} — {detail}")
    unchecked.append(name)


def fetch(url, accept="application/json", timeout=30):
    req = urllib.request.Request(url, headers={
        "Accept": accept,
        "User-Agent": "bio-dsh-public-release-audit/1.0",
    })
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.status, r.read()


def get_json(url):
    _, body = fetch(url)
    return json.loads(body)


log(f"PUBLICAUDIT-001: auditing {REPO_URL} and {CONCEPT_DOI}\n")

# ------------------------------------------------------------------ 1. repository
meta = None
try:
    meta = get_json(f"https://api.github.com/repos/{OWNER_REPO}")
    check("the cited repository resolves", True,
          f"{meta.get('full_name')} · default branch {meta.get('default_branch')} · "
          f"pushed {meta.get('pushed_at')}")
    # Print the raw field, not the bare value: "PASS ... is public — False" reads
    # as a contradiction, because the API's field is `private`, not `public`.
    check("the repository is public", not meta.get("private"),
          f"private={meta.get('private')}")
    check("the repository declares a licence",
          bool((meta.get("license") or {}).get("spdx_id")),
          f"API says {((meta.get('license') or {}).get('spdx_id'))!r}")
except urllib.error.HTTPError as e:
    check("the cited repository resolves", False, f"HTTP {e.code}")
except Exception as e:                                   # noqa: BLE001
    note("the cited repository resolves", f"{type(e).__name__}: {e}")

# ------------------------------------------------------------- 2/3. tree vs manifest
upstream = None
if meta:
    branch = meta.get("default_branch") or "main"
    try:
        tree = get_json(f"https://api.github.com/repos/{OWNER_REPO}"
                        f"/git/trees/{branch}?recursive=1")
        upstream = {n["path"] for n in tree.get("tree", []) if n["type"] == "blob"}
        check("the published tree could be listed", bool(upstream), f"{len(upstream)} files")
        missing = [p for p in REQUIRED_UPSTREAM if p not in upstream]
        check("the published tree carries every required script", not missing,
              f"missing={missing}")
        if MANIFEST.exists():
            local = {f["path"] for f in json.loads(MANIFEST.read_text())["files"]}
            absent = sorted(local - upstream)
            extra = sorted(upstream - local)
            check("every exported file is published", not absent,
                  f"{len(absent)} not upstream, e.g. {absent[:5]}")
            check("the published tree holds nothing the export dropped", not extra,
                  f"{len(extra)} stale upstream, e.g. {extra[:5]}")
        else:
            note("the published tree matches the export manifest",
                 "run scripts/47_export_public_release.py first")
    except Exception as e:                               # noqa: BLE001
        note("the published tree could be listed", f"{type(e).__name__}: {e}")

# ------------------------------------------------------------------------ 4. DOI
try:
    req = urllib.request.Request(f"https://doi.org/{CONCEPT_DOI}",
                                 headers={"User-Agent": "bio-dsh-public-release-audit/1.0"})
    with urllib.request.urlopen(req, timeout=30) as r:
        check("the concept DOI resolves", r.status == 200, r.geturl())
except urllib.error.HTTPError as e:
    check("the concept DOI resolves", False, f"HTTP {e.code}")
except Exception as e:                                   # noqa: BLE001
    note("the concept DOI resolves", f"{type(e).__name__}: {e}")

# ------------------------------------------------------------------- 5. Zenodo record
try:
    rec = get_json(f"https://zenodo.org/api/records/{CONCEPT_RECORD}")
    md = rec.get("metadata", {})
    names = [c.get("name", "") for c in md.get("creators", [])]
    check("the archived record names the two authors", len(names) == 2
          and all(re.fullmatch(r"[^,]+, .+", n) for n in names), str(names))
    check("no GitHub login is used as a creator",
          all("gdpu11" not in n for n in names), str(names))
    desc = str(md.get("description") or "")
    check("the record carries a description", len(desc) > 200, f"{len(desc)} chars")
    check("the record carries keywords", bool(md.get("keywords")),
          str(md.get("keywords")))
    lic = (md.get("license") or {}).get("id") if isinstance(md.get("license"), dict) \
        else md.get("license")
    check("the record declares the licence the repository ships",
          str(lic).lower() in ("mit", "mit-license"), repr(lic))
    log(f"      record: {md.get('title')!r} version {md.get('version')!r}")
    log(f"      NOTE  this is the LATEST archived version. A record published before "
        f"the metadata fix keeps its old values until edited in the Zenodo web UI.")
except Exception as e:                                   # noqa: BLE001
    note("the Zenodo record could be read", f"{type(e).__name__}: {e}")

# --------------------------------------------------------------------------- summary
log(f"\nPUBLICAUDIT-001: {len(fails)} failure(s), {len(unchecked)} unchecked"
    + (": " + ", ".join(fails) if fails else ""))
logf.close()
sys.exit(1 if fails else (2 if unchecked else 0))
