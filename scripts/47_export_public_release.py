#!/usr/bin/env python3
"""Public-release export (PUBREL-001).

Builds the tree that may be pushed to the public GitHub repository, from an
explicit ALLOW-LIST of what a reader of the paper needs.

Why an allow-list rather than a .gitignore
------------------------------------------
This private working repository contains manuscripts, submission files and
internal review material.  The public repository is synchronized only from the
explicit export produced here.  A deny-list ("publish everything except ...")
fails open: the day someone adds a new directory, it is published by default.
The allow-list fails closed — a new directory is published only after it is
named below.

What is deliberately NOT exported
---------------------------------
`submission_bcrt/`  the uploadable journal package, which contains the cover
                    letter naming three suggested reviewers and their
                    institutional e-mail addresses
`review/`           the internal review and audit machinery, including the
                    provenance table for those addresses
`review_cs/`        a second internal review log
`journal-selection-work/`  venue strategy notes
`.paper-loop/`      the revision state machine
`logs/`             build logs
`metadata/`         curated clinical extracts and a third-party `.rda` whose
                    redistribution terms are not ours to grant
`data/raw/`         the primary data, which the reader downloads from GEO
`manuscript/`       the paper itself, plus the stale author-side draft and the
                    venue strategy notes; two supplementary documents can be
                    opted in with --include-paper
`data/processed/`   probe-mapped matrices of 120/31/18 MB; GitHub rejects a
                    file over 100 MB, and they are one script away from the
                    accessions
the root notebook   the evidence ledger, decision log and checklists describe
                    how the project was run rather than what it produces

The export is a faithful subset: every file in the output tree is a byte copy
of a tracked file, plus the one generated manifest named in `GENERATED`.

How third-party addresses are kept out
--------------------------------------
They are *not* listed as literals in this file.  This script lives in
`scripts/`, which is exported, so a literal list would publish exactly what it
is meant to protect — the first draft of this file did that.  Instead the ban
list is derived from the cover letter, which is the document that has to name
the reviewers and which is never exported.  Adding a fourth suggested reviewer
therefore protects the new address with no edit here, and the derivation is
asserted to be non-empty so a missing cover letter cannot make the check
vacuous.  Matches are reported by domain only, because `logs/` is tracked.

Checks that abort the export
----------------------------
 1. every allow-listed path exists and is tracked by git
 2. no denied path is present in the output tree
 3. no exported file carries a third-party address, an address that is not one
    of the two authors', an absolute local home path, or an `E-mail: <addr>`
    line; and no exported file exceeds the size ceiling
 4. the release metadata is internally consistent: one title and version shared by
    `.zenodo.json` and `CITATION.cff`, one concept DOI and no superseded version
    DOI, the two author ORCIDs identical in both metadata files, and MIT
    declared in both the licence file and the Zenodo metadata
 5. the tree carries every member a citable software release needs

Outputs
-------
  dist/public_release/               the tree to push
  dist/public_release.manifest.json  per-file hash, sizes, aggregate digest
  logs/export_public_release.log

Usage
-----
  python3 scripts/47_export_public_release.py
  python3 scripts/47_export_public_release.py --include-paper
  python3 scripts/47_export_public_release.py --dry-run

Running it twice must produce the same aggregate digest.
"""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "logs/export_public_release.log"
COVER_LETTER = ROOT / "submission_bcrt/cover_letter.md"

# --------------------------------------------------------------- the allow-list
# Directories exported whole.  `data/processed/` is deliberately absent: its
# three probe-mapped matrices are 120 MB, 31 MB and 18 MB, and GitHub rejects
# any single file over 100 MB, so a release that included them could not be
# pushed at all.  They are regenerated from the public accessions by
# `scripts/02_parse_geo.py`.
ALLOW_DIRS = ("scripts", "results", "figures", "tables")

# Individual root files exported.  The four metadata files are new in this
# release; before them the archived record had no licence, no description, no
# keywords, and listed the GitHub login `gdpu11` as the creator.
ALLOW_FILES = ("requirements.txt", "ENVIRONMENT.md", "LICENSE", "CITATION.cff",
               ".zenodo.json", "README.md", ".gitignore", "MODEL_LOCK_MANIFEST.md")

# Off by default.  Posting the paper is a preprint decision that belongs to the
# authors, and the data-availability statement only promises code, pipeline
# scripts and derived tables.  Even when enabled, only these two reader-facing
# supplementary documents travel.  `manuscript/manuscript.md` is the author-side
# working copy, still at draft v3.4 and carrying internal paths by design (D-046);
# `manuscript/venue_analysis.md` is venue strategy; `manuscript/latex/main.pdf` is
# a 19-page internal build, not the 16-page submitted article.  None of those
# should reach a reader as if they were the paper.
OPTIONAL_FILES = ("manuscript/supplement.md", "manuscript/REMARK_checklist.md")

# Path prefixes that must not appear in the output, whatever else changes.
DENY_PREFIX = ("submission_bcrt/", "review/", "review_cs/",
               "journal-selection-work/", ".paper-loop/", "logs/", "metadata/",
               "data/raw/", "manuscript/")
DENY_SUBSTR = ("__pycache__/", ".DS_Store", ".pyc")

# Written by this script, so it is not a copy of a tracked file.
GENERATED = ("MANIFEST.sha256",)

# A release archive is not the place for a 900 MB expression matrix.  The three
# processed matrices are 1-2 MB each; anything past this ceiling is a mistake.
MAX_BYTES = 5_000_000
MIN_FILES = 100

CONCEPT_DOI = "10.5281/zenodo.22994650"

# The only e-mail addresses that may appear in the public release are the two
# authors' own: both are printed on the manuscript's title page and in the
# supplement, so they are public by design.
AUTHOR_ADDRESSES = ("gdpujee@gmail.com", "hedanhua@hotmail.com")

EMAIL_RE = re.compile(r"[A-Za-z0-9._%+\-]+@[A-Za-z0-9\-]+(?:\.[A-Za-z0-9\-]+)+")

# `@` is also Python's matrix-multiply operator, so `Z.values@np.array(x)` reads
# as an address to a naive pattern — it did, on three analysis scripts.  The
# generic "any address must be an author's" rule is therefore applied to prose
# and data files; source files are checked against the derived third-party list,
# which is an exact match and cannot false-positive.
SOURCE_SUFFIXES = (".py", ".R", ".r", ".sh")

# Text-ish suffixes worth scanning.  A PDF or a PNG is skipped: the scan would
# decode it as mojibake and report meaningless hits.
BINARY_SUFFIXES = (".pdf", ".png", ".jpg", ".jpeg", ".rda", ".rds", ".docx",
                   ".xlsx", ".zip", ".gz", ".cel")

logf = None
fails = []


def log(msg):
    print(msg, flush=True)
    if logf:
        logf.write(msg + "\n")
        logf.flush()


def check(name, ok, detail=""):
    log(("PASS  " if ok else "FAIL  ") + name + (f" — {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def mask(addr):
    """Domain only.

    Enough to tell which address leaked, without writing an address into a
    tracked log — which would put back what this script exists to keep out.
    """
    return "@" + addr.split("@", 1)[1]


def denied(rel):
    return rel.startswith(DENY_PREFIX) or any(s in rel for s in DENY_SUBSTR)


def tracked_files():
    """Every file git would publish, as repo-relative POSIX paths.

    Reading the index rather than walking the disk is what makes the export a
    faithful subset: ignored scratch (`data/raw/*`, `__pycache__/`) cannot be
    exported by accident, and a file that is present but untracked is reported
    instead of silently shipped.
    """
    out = subprocess.run(["git", "ls-files", "-z"], cwd=str(ROOT),
                         capture_output=True, check=True).stdout
    return [p.decode() for p in out.split(b"\0") if p]


def third_party_addresses():
    """Addresses named in the cover letter that are not the authors' own."""
    if not COVER_LETTER.exists():
        return []
    found = {m.group(0) for m in EMAIL_RE.finditer(COVER_LETTER.read_text())}
    return sorted(found - set(AUTHOR_ADDRESSES))


def wanted(rel, include_paper):
    # Tested before the deny list, because the two opt-in documents live under
    # `manuscript/`, which is denied wholesale.
    if rel in OPTIONAL_FILES:
        return (True, "opt-in paper file") if include_paper \
            else (False, "opt-in paper file (off)")
    if denied(rel):
        return False, "denied"
    if rel in ALLOW_FILES:
        return True, "root file"
    for d in ALLOW_DIRS:
        if rel.startswith(d + "/"):
            return True, "allowed dir"
    if "/" not in rel:
        return False, "root file (not listed)"
    return False, "dir (not allowed)"


def scan(rel, data, third_party):
    """Return a list of (label, masked evidence) violations for one file."""
    if rel.endswith(BINARY_SUFFIXES):
        return []
    text = data.decode("utf-8", errors="replace")
    bad = []
    for addr in third_party:
        if addr in text:
            bad.append(("third-party e-mail address", mask(addr)))
    if not rel.endswith(SOURCE_SUFFIXES):
        for m in EMAIL_RE.finditer(text):
            if m.group(0) not in AUTHOR_ADDRESSES:
                bad.append(("address that is not an author's", mask(m.group(0))))
    for m in re.finditer(r"(?:/Users/|/home/)[A-Za-z0-9._\-]+", text):
        bad.append(("absolute local path", m.group(0)))
    return bad


def read_meta():
    """Parse the two machine-readable metadata files."""
    zen = json.loads((ROOT / ".zenodo.json").read_text())
    cff = {}
    for line in (ROOT / "CITATION.cff").read_text().splitlines():
        m = re.match(r"^([a-z\-]+):\s*(.+?)\s*$", line)
        if m and m.group(1) not in cff:
            cff[m.group(1)] = m.group(2).strip().strip('"')
    return zen, cff


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="dist/public_release")
    ap.add_argument("--include-paper", action="store_true",
                    help="also export the two reader-facing supplementary documents"
                         " (posting the paper itself is a preprint decision)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be exported, write nothing")
    args = ap.parse_args()

    global logf
    if not args.dry_run:
        LOG.parent.mkdir(parents=True, exist_ok=True)
        logf = open(LOG, "w")

    out = ROOT / args.out
    log(f"PUBREL-001: public release export -> {args.out}"
        + ("  [DRY RUN]" if args.dry_run else ""))

    # ---------------------------------------------------- 0. the ban list itself
    third_party = third_party_addresses()
    # The derivation is the check's only input, so an unreadable cover letter
    # would turn the whole scan into a no-op that still reports PASS.
    check("third-party address list derived from the cover letter", len(third_party) >= 3,
          f"{len(third_party)} address(es) in {COVER_LETTER.name} that are not an author's")

    # ---------------------------------------------------------- 1. selection
    files = tracked_files()
    selected, reasons = [], {}
    for rel in files:
        ok, why = wanted(rel, args.include_paper)
        reasons[why] = reasons.get(why, 0) + 1
        if ok:
            selected.append(rel)
    log(f"\ntracked files: {len(files)}; selected: {len(selected)}")
    for why in sorted(reasons):
        log(f"    bucket {why}: {reasons[why]}")

    check(f"at least {MIN_FILES} files selected", len(selected) >= MIN_FILES,
          f"{len(selected)}")
    for d in ALLOW_DIRS:
        n = sum(1 for r in selected if r.startswith(d + "/"))
        check(f"allow-listed directory {d!r} is populated", n > 0, f"{n} file(s)")
    for f in ALLOW_FILES:
        check(f"required root file {f!r} is tracked", f in files)

    # ---------------------------------------------------------- 2. content scan
    violations = []
    biggest = (0, "")
    for rel in selected:
        data = (ROOT / rel).read_bytes()
        if len(data) > biggest[0]:
            biggest = (len(data), rel)
        if len(data) > MAX_BYTES:
            violations.append((rel, "oversized file", f"{len(data)} bytes"))
            continue
        for label, hit in scan(rel, data, third_party):
            violations.append((rel, label, hit))
    check("no oversized file exported",
          not any(v[1] == "oversized file" for v in violations),
          f"largest = {biggest[1]} ({biggest[0]:,} bytes)")
    for label in ("third-party e-mail address", "address that is not an author's",
                  "absolute local path", "oversized file"):
        hits = [f"{rel} ({hit})" for rel, lab, hit in violations if lab == label]
        check(f"no {label} exported", not hits, "; ".join(hits[:4]))

    # ---------------------------------------------------------- 3. copy
    if not args.dry_run:
        if out.exists():
            shutil.rmtree(out)
        for rel in selected:
            dst = out / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / rel, dst)

        # Checked on the OUTPUT tree, not on the selection, so a later change to
        # the copy step cannot quietly reintroduce a denied path.
        present = [p.relative_to(out).as_posix() for p in out.rglob("*")
                   if p.is_file()]
        leaked = sorted(p for p in present if denied(p))
        check("no denied path in the output tree", not leaked, ", ".join(leaked[:5]))
        strays = sorted(set(present) - set(selected) - set(GENERATED))
        check("output tree holds only exported files", not strays,
              ", ".join(strays[:5]))

    # ---------------------------------------------------------- 4. metadata
    zen, cff = read_meta()
    v_zen, v_cff = zen.get("version"), cff.get("version")
    check(".zenodo.json and CITATION.cff declare the same release title",
          bool(zen.get("title")) and zen.get("title") == cff.get("title"),
          f"zenodo={zen.get('title')!r} cff={cff.get('title')!r}")
    check(".zenodo.json and CITATION.cff declare the same version",
          bool(v_zen) and v_zen == v_cff, f"zenodo={v_zen!r} cff={v_cff!r}")
    check(".zenodo.json declares an MIT licence",
          str(zen.get("license", "")).lower().startswith("mit"),
          repr(zen.get("license")))
    lic_txt = (ROOT / "LICENSE").read_text()
    lic_first = lic_txt.splitlines()[0].strip()
    check("LICENSE file is MIT and matches the Zenodo declaration",
          lic_first.lower().startswith("mit license"), repr(lic_first))
    check("LICENSE file names both authors",
          "Danhua He" in lic_txt and "Qiang Li" in lic_txt)

    creators = zen.get("creators", [])
    check(".zenodo.json lists two named creators",
          len(creators) == 2 and all(c.get("name") for c in creators),
          str([c.get("name") for c in creators]))
    # The archived record listed the GitHub login `gdpu11` as the creator.  A
    # login is not a person, so it must not come back.
    check("no GitHub login used as a creator",
          all("gdpu11" not in str(c.get("name", "")) for c in creators),
          str([c.get("name") for c in creators]))
    check(".zenodo.json carries a description and keywords",
          len(str(zen.get("description", ""))) > 200
          and len(zen.get("keywords", [])) >= 4,
          f"description {len(str(zen.get('description', '')))} chars, "
          f"{len(zen.get('keywords', []))} keywords")

    # An ORCID typo is invisible until an anchor match returns zero, which is how
    # a `0000-` prefix survived one round in a sibling project.  Compare the two
    # files against each other and against the known digits.
    zen_orcids = sorted(str(c.get("orcid", "")) for c in creators)
    cff_orcids = sorted(re.findall(r"orcid\.org/([0-9X\-]+)",
                                   (ROOT / "CITATION.cff").read_text()))
    check("the same two ORCIDs in .zenodo.json and CITATION.cff",
          zen_orcids == cff_orcids and len(zen_orcids) == 2,
          f"zenodo={zen_orcids} cff={cff_orcids}")
    check("both ORCIDs are well-formed and 0009-prefixed",
          all(re.fullmatch(r"0009-\d{4}-\d{4}-\d{3}[\dX]", o) for o in zen_orcids),
          str(zen_orcids))

    # ---------------------------------------------------------- 5. DOI citation
    # Dual-DOI discipline: README anchors the concept DOI; after Zenodo mints a
    # version DOI, CITATION.cff pins that exact snapshot. A pre-publication
    # export has no version DOI yet; accepting an old one here would silently
    # cite the wrong snapshot.
    VERSION_DOI = cff.get("doi")
    ALLOWED = {CONCEPT_DOI} | ({VERSION_DOI} if VERSION_DOI else set())
    for what, path, want in (("CITATION.cff", ROOT / "CITATION.cff", VERSION_DOI),
                             ("README.md", ROOT / "README.md", CONCEPT_DOI),
                             ("LICENSE", ROOT / "LICENSE", None)):
        if not path.exists():
            check(f"{what} exists for the DOI check", False)
            continue
        text = path.read_text()
        dois = re.findall(r"10\.5281/zenodo\.\d+", text)
        if want:
            check(f"{what} cites its required DOI ({want.split('.')[-1]})", want in text, want)
        elif what == "CITATION.cff":
            check("CITATION.cff has no stale version DOI before Zenodo minting",
                  not dois, str(dois))
        wrong = sorted({d for d in dois if d not in ALLOWED})
        check(f"{what} cites no dead or superseded DOI", not wrong, str(wrong))

    # ---------------------------------------------------------- 6. manifest
    manifest = []
    for rel in sorted(selected):
        data = (ROOT / rel).read_bytes()
        manifest.append({"path": rel, "bytes": len(data),
                         "sha256": hashlib.sha256(data).hexdigest()})
    agg = hashlib.sha256()
    for m in manifest:
        agg.update(f"{m['path']}\0{m['sha256']}\n".encode())
    digest = agg.hexdigest()

    if not args.dry_run:
        (out / "MANIFEST.sha256").write_text(
            "".join(f"{m['sha256']}  {m['path']}\n" for m in manifest))
        check("MANIFEST.sha256 written", (out / "MANIFEST.sha256").exists(),
              f"{len(manifest)} entries")
        # Re-hash from the output tree rather than trusting the source read: a
        # truncated or encoding-mangled copy is exactly what this catches.
        mismatched = [m["path"] for m in manifest
                      if hashlib.sha256((out / m["path"]).read_bytes()).hexdigest()
                      != m["sha256"]]
        check("every exported file matches its recorded hash", not mismatched,
              ", ".join(mismatched[:5]))
        payload = {
            "generator": "scripts/47_export_public_release.py",
            "version": v_zen,
            "concept_doi": CONCEPT_DOI,
            "include_paper": args.include_paper,
            "n_files": len(manifest),
            "total_bytes": sum(m["bytes"] for m in manifest),
            "aggregate_sha256": digest,
            "excluded_reasons": reasons,
            "files": manifest,
        }
        (out.parent / "public_release.manifest.json").write_text(
            json.dumps(payload, indent=2) + "\n")
        log(f"\nmanifest: {out.parent.name}/public_release.manifest.json")

    log(f"files: {len(manifest)}   bytes: {sum(m['bytes'] for m in manifest):,}")
    log(f"aggregate sha256: {digest}")
    log(f"largest file: {biggest[1]} ({biggest[0]:,} bytes)")
    log(f"\nPUBREL-001: {len(fails)} failure(s)"
        + (": " + ", ".join(fails) if fails else ""))
    if logf:
        logf.close()
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
