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
  1. the cited repository URL resolves and GitHub recognizes its MIT licence
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
import hashlib
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "logs/public_release_audit.log"

OWNER_REPO = "gdpujee/bcrt-immune-metabolic-score"
REPO_URL = f"https://github.com/{OWNER_REPO}"
CONCEPT_DOI = "10.5281/zenodo.22994650"
CONCEPT_RECORD = CONCEPT_DOI.rsplit(".", 1)[1]
VERSION_DOI = re.search(r'(?m)^doi:\s*"?([^"\s]+)', LOCAL_CFF)
VERSION_DOI = VERSION_DOI.group(1) if VERSION_DOI else None
VERSION = LOCAL_ZENODO.get("version", "")
GENERATED_UPSTREAM = {"MANIFEST.sha256"}
MANIFEST = ROOT / "dist/public_release.manifest.json"
LOCAL_ZENODO = json.loads((ROOT / ".zenodo.json").read_text())
LOCAL_CFF = (ROOT / "CITATION.cff").read_text()
LOCAL_CFF_TITLE = re.search(r'(?m)^title:\s*"?(.+?)"?\s*$', LOCAL_CFF).group(1)
PAPER_TITLE = re.search(
    r"(?m)^Analysis code and derived data tables accompanying the manuscript \*(.+?)\*",
    (ROOT / "README.md").read_text()).group(1)

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


def git_blob_sha(data):
    """Return Git's SHA-1 for a blob without writing it into a repository."""
    header = b"blob " + str(len(data)).encode("ascii") + b"\0"
    return hashlib.sha1(header + data).hexdigest()


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
    description = str(meta.get("description") or "")
    check("GitHub repository description names the current paper title",
          PAPER_TITLE in description and "manuscript" not in description.lower(),
          repr(description))
    license_info = meta.get("license") or {}
    spdx_id = license_info.get("spdx_id")
    check("GitHub recognizes the repository licence as MIT",
          str(spdx_id).upper() == "MIT",
          f"API says {license_info.get('name')!r} / {spdx_id!r}")
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
        upstream_blobs = {n["path"]: n.get("sha") for n in tree.get("tree", [])
                          if n["type"] == "blob"}
        upstream = set(upstream_blobs)
        check("the published tree could be listed", bool(upstream), f"{len(upstream)} files")
        missing = [p for p in REQUIRED_UPSTREAM if p not in upstream]
        check("the published tree carries every required script", not missing,
              f"missing={missing}")
        if MANIFEST.exists():
            # The hash manifest is generated after the selected payload is
            # copied, so it is intentionally absent from public_release.manifest.json.
            # It is nevertheless part of the release tree and must be current
            # upstream, not misreported as an orphaned stale path.
            local = ({f["path"] for f in json.loads(MANIFEST.read_text())["files"]}
                     | GENERATED_UPSTREAM)
            absent = sorted(local - upstream)
            extra = sorted(upstream - local)
            check("every exported file is published", not absent,
                  f"{len(absent)} not upstream, e.g. {absent[:5]}")
            release_dir = ROOT / "dist/public_release"
            mismatched = []
            for rel in sorted(local & upstream):
                path = release_dir / rel
                if not path.is_file() or git_blob_sha(path.read_bytes()) != upstream_blobs[rel]:
                    mismatched.append(rel)
            check("published file contents match the local export", not mismatched,
                  f"{len(mismatched)} content mismatch(es), e.g. {mismatched[:5]}")
            check("the published tree holds nothing the export dropped", not extra,
                  f"{len(extra)} stale upstream, e.g. {extra[:5]}")
        else:
            note("the published tree matches the export manifest",
                 "run scripts/47_export_public_release.py first")
        raw = f"https://raw.githubusercontent.com/{OWNER_REPO}/{branch}/"
        try:
            remote_readme = fetch(raw + "README.md", accept="text/plain")[1].decode(
                "utf-8", errors="replace")
            check("public README names the current manuscript title",
                  PAPER_TITLE in remote_readme, repr(PAPER_TITLE))
            remote_dois = set(re.findall(r"10\.5281/zenodo\.\d+", remote_readme))
            check("public README cites the current concept DOI",
                  CONCEPT_DOI in remote_dois, str(sorted(remote_dois)))
        except Exception as e:                         # noqa: BLE001
            note("the public README title could be checked", f"{type(e).__name__}: {e}")
        try:
            remote_cff = fetch(raw + "CITATION.cff", accept="text/plain")[1].decode(
                "utf-8", errors="replace")
            title_match = re.search(r'(?m)^title:\s*"?(.+?)"?\s*$', remote_cff)
            check("public CITATION.cff title matches the current release",
                  bool(title_match) and title_match.group(1) == LOCAL_CFF_TITLE,
                  f"expected={LOCAL_CFF_TITLE!r}, observed="
                  f"{title_match.group(1) if title_match else 'missing'}")
            cff_dois = set(re.findall(r"10\.5281/zenodo\.\d+", remote_cff))
            if VERSION_DOI:
                check("public CITATION.cff cites the current version DOI",
                      VERSION_DOI in cff_dois, str(sorted(cff_dois)))
        except Exception as e:                         # noqa: BLE001
            note("the public citation title could be checked", f"{type(e).__name__}: {e}")
        try:
            remote_zenodo = get_json(raw + ".zenodo.json")
            check("public .zenodo.json title matches the current release",
                  remote_zenodo.get("title") == LOCAL_ZENODO.get("title"),
                  f"expected={LOCAL_ZENODO.get('title')!r}, "
                  f"observed={remote_zenodo.get('title')!r}")
        except Exception as e:                         # noqa: BLE001
            note("the public Zenodo metadata file title could be checked",
                 f"{type(e).__name__}: {e}")
        try:
            latest = get_json(f"https://api.github.com/repos/{OWNER_REPO}/releases/latest")
            expected_release_name = f"{VERSION}: {LOCAL_ZENODO['title']}"
            check("latest GitHub release title matches the current version",
                  latest.get("name") == expected_release_name,
                  f"expected={expected_release_name!r}, observed={latest.get('name')!r}")
            check("latest GitHub release uses the manuscript version tag",
                  latest.get("tag_name") == VERSION,
                  f"expected={VERSION!r}, observed={latest.get('tag_name')!r}")
        except Exception as e:                         # noqa: BLE001
            note("the latest GitHub release could be checked",
                 f"{type(e).__name__}: {e}")
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
    check("Zenodo record title matches the current release metadata",
          md.get("title") == LOCAL_ZENODO.get("title"),
          f"expected={LOCAL_ZENODO.get('title')!r}, observed={md.get('title')!r}")
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
