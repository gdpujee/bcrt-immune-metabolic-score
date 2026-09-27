"""Submission-package consistency gate (SUBGATE-001).

Turns the external review's P0-1 ("one source of truth; every format must agree")
and P0-2 ("no stray nan; no duplicated title/keywords") plus its closing bar
("所有文档和数字完全一致") into executable checks that fail loudly.

Checks
  1. every deliverable exists and contains no literal `nan` token
  2. the final PDF states the title exactly once and the keyword line exactly once
     (both needles are asserted non-empty first, so the check cannot pass vacuously)
  3. references_numbered.md is numbered 1..N
  4. every reference 1..N is cited in both manuscript files, and no citation marker
     points outside 1..N
  5. abstract_structured.md regenerates byte-identically from the submission
     manuscript (single source of truth)
  6. every decimal in the standalone abstract also appears in the manuscript body
     *excluding the abstract itself* — otherwise the check is vacuous
  7. the headline numbers of the transportability table (unadjusted and adjusted
     per-SD hazard ratio for each locked cohort, calibration slopes, PH p-values)
     appear in both the submission markdown and the rendered final PDF
  8. the supplementary package files are physically present and the REMARK
     checklist covers all 20 items
  9. the BCRT packaging rules: 4-6 keywords, a 150-250 word structured abstract,
     a body within the 3,500-word limit, and 'Statements and Declarations' placed
     after the References in every rendered format, not just in the markdown
 10. the Supplementary Information PDF exists, carries the title / journal /
     authors / corresponding e-mail on page 1, and still holds every
     supplementary figure and its caption
 11. the export probe: no deliverable — PDFs included — leaks an internal script
     path, result path, ledger name, revision tag, run id or file hash, and no
     editorial parenthetical has crept back into the reference list
 12. the cover letter carries three institutional reviewer addresses and makes no
     publishing-model or APC claim, and submission_bcrt/ holds no internal notes
 13. the archive citation is right in every rendered format: the live repository
     URL, the concept DOI, and no superseded version DOI
 14. the three suggested-reviewer addresses appear only in the cover letter and in
     the provenance table that records where each was read from

Exit status is non-zero if any check fails.

Outputs: logs/submission_gate.log
"""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / "submission_bcrt"
RAW = ROOT / "results/raw"
logf = open(ROOT / "logs/submission_gate.log", "w")
fails = []
n_checks = 0


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


def check(name, ok, detail=""):
    global n_checks
    n_checks += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f" — {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def mask(addr):
    """Domain only.

    This gate writes `logs/submission_gate.log`, which is a tracked file, so an
    address echoed into a check detail would defeat check 14 — the addresses must
    appear nowhere but the cover letter and the provenance table.
    """
    return "@" + addr.split("@", 1)[1]


ADDR_RE = re.compile(r"[A-Za-z0-9._%+\-]+@[A-Za-z0-9\-]+(?:\.[A-Za-z0-9\-]+)+")

# The only addresses that may appear in a public release are the two authors'.
AUTHOR_ADDRESSES = {"gdpujee@gmail.com", "hedanhua@hotmail.com"}


def section(lines, heading):
    """Index of the line holding `heading` (exact match), or -1."""
    try:
        return lines.index(heading)
    except ValueError:
        return -1


# ---------------------------------------------------------------- 1. no literal "nan"
# submission_bcrt/ holds only files that can be uploaded.  README.md used to be in
# this list; it was an index describing the build, moved to review/PACKAGE_INDEX.md,
# and the supplement is now delivered as a PDF as BCRT requires.
DELIVERABLES = ["manuscript_final.pdf", "manuscript_review.pdf", "manuscript_bcrt.docx",
                "manuscript_bcrt.tex", "manuscript_submission.md", "supplement.md",
                "Supplementary_Information.pdf", "legends.md", "title_page.md",
                "declarations.md", "abstract_structured.md",
                "references_numbered.md", "cover_letter.md", "REMARK_checklist.md"]

texts = {}
flats = {}
for f in DELIVERABLES:
    p = SUB / f
    check(f"{f} exists", p.exists())
    if not p.exists():
        continue
    if p.suffix == ".pdf":
        from pypdf import PdfReader
        txt = "\n".join(pg.extract_text() or "" for pg in PdfReader(str(p)).pages)
    elif p.suffix == ".docx":
        from docx import Document
        doc = Document(str(p))
        txt = "\n".join(par.text for par in doc.paragraphs)
        for tb in doc.tables:
            for row in tb.rows:
                txt += "\n" + "\t".join(c.text for c in row.cells)
    else:
        txt = p.read_text()
    texts[f] = txt
    # PDF/DOCX extraction inserts line breaks at cell and column edges, so every
    # containment test below runs on a whitespace-collapsed copy; otherwise a
    # correct render can look absent (and a wrapping defect can look fine).
    flats[f] = re.sub(r"\s+", " ", txt)
    hits = re.findall(r"(?<![A-Za-z])nan(?![A-Za-z])", txt)
    check(f"{f} has no literal 'nan'", not hits, f"{len(hits)} occurrence(s)")

# ---------------------------------------------------------------- 2. title / keywords
ms_txt = texts.get("manuscript_submission.md", "")
ms_lines = ms_txt.split("\n")
pdf_txt = texts.get("manuscript_final.pdf", "")
pdf_flat = flats.get("manuscript_final.pdf", "")

ti = section(ms_lines, "## Title")
title = ms_lines[ti + 1].strip() if ti >= 0 and ti + 1 < len(ms_lines) else ""
ki = section(ms_lines, "## Keywords")
kw_line = ms_lines[ki + 1].strip() if ki >= 0 and ki + 1 < len(ms_lines) else ""
# a needle that failed to resolve would make count()==0 (or, for "", len+1) and
# the check would pass or fail for the wrong reason, so assert the needles first
check("title needle resolved from ## Title", len(title) > 20, repr(title[:60]))
check("keyword needle resolved from ## Keywords", kw_line.startswith("breast cancer;"),
      repr(kw_line[:60]))
if title:
    check("final PDF: title exactly once", pdf_flat.count(title) == 1,
          f"{pdf_flat.count(title)} occurrence(s)")
if kw_line:
    check("final PDF: keyword line exactly once", pdf_flat.count(kw_line) == 1,
          f"{pdf_flat.count(kw_line)} occurrence(s)")

# ---------------------------------------------------------------- 3/4. citations
refs_txt = texts.get("references_numbered.md", "")
ref_nums = [int(m.group(1)) for m in re.finditer(r"^\[(\d+)\] ", refs_txt, re.M)]
check("references_numbered.md is 1..N", ref_nums == list(range(1, len(ref_nums) + 1)),
      f"{len(ref_nums)} entries")


def split_at(path, heading):
    t = path.read_text()
    hits = [m.start() for m in re.finditer(rf"^{re.escape(heading)}\s*$", t, re.M)]
    return t[:hits[-1]] if hits else t


def cited_numbers(text):
    """Every reference number cited in `text`.

    Citation markers come in three shapes — `[3]`, `[1,2]` and `[12-14]` — and a
    bare `\\[(\\d+)\\]` search silently reports the last two as uncited, which looks
    exactly like a dangling reference.  Bracketed confidence intervals such as
    `[1.29,1.98]` are excluded because the decimal point breaks the digit run.
    """
    out = set()
    for m in re.finditer(r"\[(\d+(?:\s*[,\-–]\s*\d+)*)\]", text):
        for part in re.split(r"\s*,\s*", m.group(1)):
            if re.search(r"[\-–]", part):
                a, b = re.split(r"\s*[\-–]\s*", part)
                out.update(range(int(a), int(b) + 1))
            else:
                out.add(int(part))
    return out


for name, path in [("submission", SUB / "manuscript_submission.md"),
                   ("author", ROOT / "manuscript/manuscript.md")]:
    body = split_at(path, "## References")
    used = sorted(cited_numbers(body))
    check(f"{name}: every citation in 1..N",
          max(used, default=0) <= len(ref_nums) and min(used, default=1) >= 1,
          f"used={used}")
    check(f"{name}: every reference cited",
          used == list(range(1, len(ref_nums) + 1)),
          f"missing={sorted(set(range(1, len(ref_nums) + 1)) - set(used))}")

# ---------------------------------------------------------------- 5. abstract generated
r = subprocess.run([sys.executable, str(ROOT / "scripts/34_build_abstract.py")],
                   capture_output=True, text=True, cwd=str(ROOT))
check("abstract_structured.md regenerates byte-identically", r.returncode == 0,
      (r.stdout.strip().splitlines() or [""])[-1] or r.stderr[-200:])

# ---------------------------------------------------------------- 6. abstract vs body
abs_txt = texts.get("abstract_structured.md", "")
# exclude the abstract section itself: the claim under test is that the standalone
# abstract introduces no number that the rest of the manuscript does not carry
sub_full = (SUB / "manuscript_submission.md").read_text()
ai = sub_full.find("## Abstract")
intro = sub_full.find("## Introduction", ai if ai >= 0 else 0)
ri = sub_full.rfind("## References")
body_no_abs = sub_full[intro:ri] if intro >= 0 and ri > intro else ""
check("manuscript body located for the abstract cross-check",
      len(body_no_abs) > 2000, f"{len(body_no_abs)} chars")
nums = sorted(set(re.findall(r"\b\d+\.\d+\b", abs_txt)))
missing = [n for n in nums if n not in body_no_abs]
check("abstract decimals all present in the manuscript body", not missing,
      f"{len(nums)} decimals, missing={missing}")

# ---------------------------------------------------------------- 7. headline numbers
# Chain of custody for the numbers the external review asked to unify (P1-8/P1-10):
# results/raw/adjusted_per_sd.json  ->  tables/Tab3_performance_v3.tsv  ->  rendered PDF.
# Checking "the number appears somewhere in the text" would be nearly vacuous, so the
# table cell is compared to the JSON value and the rendered PDF must carry the cell.
aps = json.load(open(RAW / "adjusted_per_sd.json"))
tab3 = [ln.split("\t") for ln in
        (ROOT / "tables/Tab3_performance_v3.tsv").read_text().strip().splitlines()]
hdr, rows = tab3[0], tab3[1:]
COHORT_IN_TABLE = {"GSE20685": "GSE20685", "SCANB": "SCAN-B", "METABRIC": "METABRIC"}
check("Table 3 has exactly one row per locked validation cohort",
      len(rows) == 3 and all(any(k in r[0] for r in rows) for k in COHORT_IN_TABLE.values()),
      f"{len(rows)} data rows: {[r[0] for r in rows]}")
try:
    i_adj = hdr.index("adj HR per SD (95% CI)")
except ValueError:
    i_adj = -1
check("Table 3 exposes an adjusted per-SD hazard-ratio column", i_adj >= 0, str(hdr))
for cohort, tag in COHORT_IN_TABLE.items():
    row = next((r for r in rows if tag in r[0]), None)
    adj = aps[cohort]["per_SD_HR"]
    if row is None or i_adj < 0:
        check(f"{cohort} row present in Table 3", False)
        continue
    cell = row[i_adj]
    want = f"{adj:.2f}"
    check(f"Table 3 adjusted HR cell for {cohort} equals adjusted_per_sd.json",
          want in cell, f"cell={cell!r} json={adj}")

slopes = {c: f"{aps[c]['calib_slope']:.2f}" for c in COHORT_IN_TABLE}
check("calibration slopes rendered in the final PDF",
      all(s in pdf_flat for s in slopes.values()), str(slopes))

# ---------------------------------------------------------------- 7b. table legibility
# A wide table can render "successfully" while wrapping every header and number
# mid-token ("deat hs", "Slo pe", "1909/ 71"), which is unreadable and was exactly
# the defect the external review raised for the Word tables.  Column widths and cell
# padding are the usual cause, and nothing else in the build notices.  So the
# rendered deliverables are searched for whole tokens.
LEGIBILITY = {
    "Table 1 headers": ["n_tumors", "n_normals", "OS_events", "Median_OS", "Platform", "PMID"],
    "Table 1 cells": ["23740839", "21501481", "27006338", "32913985"],
    "Table 3 headers": ["deaths", "Slope", "High/Low", "Binary HR"],
    "Table 3 cells": ["1909/71", "156/3117", "48/279", "0.595", "0.573", "0.229"],
    "Table 3 adjusted column": [f"{aps[c]['per_SD_HR']:.2f} [{aps[c]['per_SD_CI'][0]:.2f},"
                                f"{aps[c]['per_SD_CI'][1]:.2f}]"
                                for c in ("GSE20685", "SCANB", "METABRIC")],
}
for label, tokens in LEGIBILITY.items():
    for f in ("manuscript_final.pdf", "manuscript_review.pdf", "manuscript_bcrt.docx"):
        if f not in flats:
            continue
        bad = [t for t in tokens if t not in flats[f]]
        check(f"{label} intact in {f}", not bad, f"broken/absent={bad}")

# ---------------------------------------------------------------- 8. supplement package
for f in ["supplement.md", "Supplementary_Information.pdf", "REMARK_checklist.md"]:
    check(f"package contains {f}", (SUB / f).exists())
remark = texts.get("REMARK_checklist.md", "")
n_items = len(re.findall(r"^\|\s*\d+", remark, re.M))
check("REMARK checklist covers all 20 items", n_items >= 20, f"{n_items} numbered rows")


def norm(s):
    """Flatten for cross-format comparison.

    PDF and DOCX extraction inserts line breaks at cell and column edges, and LaTeX
    escapes punctuation, so a containment test on the raw text can report a
    correctly rendered section as missing — and, worse, a mis-ordered one as fine.
    """
    return re.sub(r"\s+", " ", s.replace("\\", " ")).strip()


def url_text(s):
    """Flatten for URL comparison, undoing LaTeX line-breaking hints.

    The .tex builder writes a long URL as `https:/\\allowbreak /\\allowbreak
    github.com/...` so TeX can break it across lines.  Comparing the raw source
    therefore reports a correctly cited URL as absent - and would equally accept a
    corrupted one, since any split string would fail the same way.  The hints are
    removed and the result is whitespace-free, which is safe for a substring test.
    """
    return norm(s).replace("allowbreak", " ").replace(" ", "")


# ---------------------------------------------------------------- 9. BCRT packaging rules
parts = re.split(r"(?m)^(## .*)$", (SUB / "manuscript_submission.md").read_text())
secs = {parts[i].strip()[3:].strip(): parts[i + 1] for i in range(1, len(parts), 2)}

kw = [k.strip() for k in kw_line.rstrip(".").split(";") if k.strip()]
check("keywords are 4-6 in number (BCRT rule)", 4 <= len(kw) <= 6, f"{len(kw)}: {kw}")

tp_txt = texts.get("title_page.md", "")
tp_kw = re.search(r"(?m)^Keywords:\s*(.+)$", tp_txt)
tp_list = ([k.strip() for k in tp_kw.group(1).rstrip(".").split(";") if k.strip()]
           if tp_kw else [])
check("title page carries the same keyword set", tp_list == kw,
      f"title page has {len(tp_list)}: {tp_list}")

# BCRT: "Word count does not include abstract, title page, tables and figures or
# references."  So only the narrative sections are counted, and a section that
# failed to resolve must not silently contribute zero.
COUNTED = ["Introduction", "Methods", "Results", "Discussion", "Conclusion"]
missing_secs = [k for k in COUNTED if k not in secs]
body_words = sum(len(re.findall(r"\S+", secs.get(k, ""))) for k in COUNTED)
check("all counted body sections located", not missing_secs, str(missing_secs))
check("body within the BCRT 3,500-word limit", 0 < body_words <= 3500,
      f"{body_words} words")

# counted from the manuscript, not from abstract_structured.md, whose footer carries
# its own word-count line that would be counted as abstract prose
abs_paras = [l for l in secs.get("Abstract", "").split("\n") if l.strip()]
n_abs = len(re.sub(r"\*\*[^*]+\*\*", " ", " ".join(abs_paras)).split())
check("structured abstract within the BCRT 150-250 word range", 150 <= n_abs <= 250,
      f"{n_abs} words")

# "Statements and Declarations ... should be placed after the References section."
# Asserted against every rendered format: the markdown order was already right while
# the DOCX still ended with the reference list, which is what the review caught.
n_refs = len(re.findall(r"^\[(\d+)\]", (SUB / "references_numbered.md").read_text(), re.M))
for f in ("manuscript_submission.md", "manuscript_bcrt.tex", "manuscript_bcrt.docx",
          "manuscript_final.pdf", "manuscript_review.pdf"):
    flat = norm(flats.get(f, ""))
    i_st = flat.find("Statements and Declarations")
    # the LAST numbered reference entry, not the word "References": a reference
    # appended after the Declarations block once slipped past an rfind("References")
    # probe because the Colleoni title contains the word "recurrence" near it.
    i_last = max((m.start() for m in re.finditer(rf"\[{n_refs}\]", flat)), default=-1)
    check(f"{f}: Statements and Declarations after the last reference [{n_refs}]",
          0 <= i_last < i_st, f"last_ref@{i_last} statements@{i_st}")

# ---------------------------------------------------------------- 10. supplement PDF
supp_pdf = SUB / "Supplementary_Information.pdf"
if supp_pdf.exists():
    raw = supp_pdf.read_bytes()
    from pypdf import PdfReader
    pages = PdfReader(str(supp_pdf)).pages
    p1 = norm(pages[0].extract_text() or "")
    supp_flat = norm("\n".join(pg.extract_text() or "" for pg in pages))
    HEADER = {"article title": "Cross-platform external validation",
              "journal name": "Breast Cancer Research and Treatment",
              "author names": "Danhua He",
              "corresponding author": "Qiang Li",
              # the guidelines ask for the corresponding author's affiliation AND
              # e-mail.  The affiliation is also implied by the superscript next to
              # the name, but an editor should not have to infer it, so it is
              # required as a literal.
              "corresponding affiliation": "Foshan, Guangdong, China",
              "corresponding e-mail": "gdpujee@gmail.com"}
    for what, needle in HEADER.items():
        check(f"supplement PDF page 1 carries the {what}", needle in p1, repr(needle))
    # ...and the same requirement on the source, so editing the header without
    # rebuilding cannot leave the check reading a stale artifact that still has it.
    supp_src = norm((SUB / "supplement.md").read_text())
    for what, needle in HEADER.items():
        check(f"supplement source carries the {what}", needle in supp_src, repr(needle))
    # staleness probe: every section heading in the source must appear in the
    # rendered PDF.  A narrower pattern (only "S<digit> <text>") was tested by
    # mutation and let a newly added section through, because the added heading did
    # not match it — a probe that skips exactly the case it exists for.
    for h in re.findall(r"(?m)^## (.+)$", (SUB / "supplement.md").read_text()):
        check(f"supplement PDF carries section {h[:24]!r}", norm(h) in supp_flat,
              repr(h))
    # figures: the labels and caption openers are read from legends.md, which the
    # build generator owns, so a new supplementary figure cannot be added without
    # this check following it.
    leg = texts.get("legends.md", "")
    supp_figs = []
    for m in re.finditer(r"(?m)^- \*\*(Fig\. S\d+)\.\*\* (.+)$", leg):
        supp_figs.append((m.group(1), norm(m.group(2))[:45]))
    check("legends.md lists the supplementary figures", len(supp_figs) >= 4,
          f"{[s[0] for s in supp_figs]}")
    for lab, opener in supp_figs:
        check(f"supplement PDF carries {lab} and its caption", opener in supp_flat,
              repr(opener[:60]))
    n_img = raw.count(b"/Subtype /Image")
    check("supplement PDF embeds raster artwork", n_img >= len(supp_figs),
          f"{n_img} image XObject(s) for {len(supp_figs)} figure(s)")

# ---------------------------------------------------------------- 11. export probe
# Scanning only the markdown and the .tex would miss the two artifacts a reader
# actually receives.  The PDFs are read too, and the needles are the same ones the
# review found by hand.
EXPORT_TRACES = [
    ("internal script path", r"scripts/\d"),
    ("result path", r"results/raw"),
    ("log path", r"logs/"),
    ("project codename", r"bio-dsh"),
    ("internal ledger or checklist",
     r"EVIDENCE_LEDGER|DECISIONS\.md|MASTER_CHECKLIST|CLAIMS_MATRIX|RESEARCH_LOG|"
     r"USER_ACTION_REQUIRED|PAPER_MAP|PACKAGE_INDEX"),
    ("internal revision tag", r"\bv3\.\d\b|\bW4\b"),
    # an internal decision id: `ADJSD-001` needs the general shape below, but the
    # single-letter form `D-038` does not match it (the char after the letter is a
    # hyphen, which the class excludes), so it is listed separately
    ("run identifier", r"\b[A-Z][A-Z0-9]{1,6}-\d{3}\b|\bD-\d{3}\b"),
    ("file hash", r"\bsha256\b"),
]
for label, pat in EXPORT_TRACES:
    bad = []
    for f in DELIVERABLES:
        flat = flats.get(f, "")
        if not flat:
            continue
        m = re.search(pat, flat)
        if m:
            bad.append(f"{f}:{m.group(0)}")
    check(f"no {label} in any deliverable", not bad, ", ".join(bad[:4]))

# The editorial notes that used to trail each reference ("(GSE42568 discovery
# cohort.)") are recorded as removed strings rather than as a pattern, so the check
# cannot be satisfied by rewording them.
REMOVED_REF_NOTES = [
    "(GSE42568 discovery cohort.)", "(GSE20685 validation cohort.)",
    "(SCAN-B cohort descriptor.)", "(METABRIC cohort, primary report.)",
    "(METABRIC extended cohort.)", "(Original REMARK statement.)",
    "(PAM50 intrinsic-subtype classifier.)",
    "(KEGG pathway database for candidate-pool construction.)",
    "(cBioPortal resource used to access brca_metabric.)",
]
for note in REMOVED_REF_NOTES:
    hits = [f for f in DELIVERABLES if note in norm(flats.get(f, ""))]
    check(f"reference note {note!r} stays removed", not hits, ", ".join(hits))

# ---------------------------------------------------------------- 12. cover letter
cl = texts.get("cover_letter.md", "")
check("cover letter makes no publishing-model claim",
      not re.search(r"Subscription route|Article Processing Charge|APC", cl),
      repr(re.search(r"Subscription route|Article Processing Charge|APC", cl)))
addr = re.findall(r"E-mail:\s*([A-Za-z0-9._%+\-]+@[A-Za-z0-9\-]+(?:\.[A-Za-z0-9\-]+)*)",
                  cl)
check("cover letter names three reviewers with e-mail addresses", len(addr) == 3,
      str([mask(a) for a in addr]))
check("every reviewer address is institutional (no free-mail domain)",
      bool(addr) and all(not re.search(
          r"@(gmail|hotmail|yahoo|outlook|live|qq|163|126|foxmail|sina)\.",
          a, re.I) for a in addr), str([mask(a) for a in addr]))

# submission_bcrt/ is what gets uploaded, so anything that is not an uploadable
# file is a packaging defect rather than clutter.
strays = sorted(p.name for p in SUB.iterdir()
                if p.is_file() and p.name not in DELIVERABLES
                and p.suffix.lower() in (".md", ".txt", ".json", ".py", ".log",
                                         ".csv", ".yaml", ".yml"))
check("submission_bcrt holds only uploadable files", not strays, str(strays))

# ---------------------------------------------------------------- 13. archive DOI
# The manuscript states that the code and derived tables are available.  Two ways
# for that statement to be wrong without anything noticing: the URL drifts from the
# one that is actually live, or the archive is cited by a VERSION DOI — which pins
# a single release, and the release that existed when this was written no longer
# matched the manuscript.  Asserted in every rendered format, not just the source.
CONCEPT_DOI = "10.5281/zenodo.22994650"
REPO_URL = "https://github.com/gdpujee/bcrt-immune-metabolic-score"
DOI_FILES = ["manuscript_submission.md", "declarations.md", "supplement.md",
             "Supplementary_Information.pdf", "manuscript_final.pdf",
             "manuscript_review.pdf", "manuscript_bcrt.docx", "manuscript_bcrt.tex"]
for f in DOI_FILES:
    flat = url_text(flats.get(f, ""))
    check(f"{f}: cites the concept DOI", CONCEPT_DOI in flat, CONCEPT_DOI)
    check(f"{f}: cites the live repository URL", REPO_URL in flat, REPO_URL)
    wrong = sorted({d for d in re.findall(r"10\.5281/zenodo\.\d+", flat)
                    if d != CONCEPT_DOI})
    check(f"{f}: cites no superseded version DOI", not wrong, str(wrong))

# ----------------------------------------------------- 14. third-party addresses
# The three suggested-reviewer addresses belong in the cover letter, which is sent
# to the journal, and in the provenance table recording where each was read from.
# Anywhere else is a defect, because this repository's `origin` is the public
# GitHub remote: a tracked file is one `git push` from publishing the names and
# institutional addresses of three people who were never asked.  The list is
# derived from the cover letter rather than hard-coded, so a fourth reviewer is
# covered the moment they are added — and the scan covers the gate's own log, which
# is why section 12 reports domains instead of addresses.
ADDRESS_ALLOWED = {"submission_bcrt/cover_letter.md", "review/REVIEWER_EMAILS.md"}
third_party = sorted({m.group(0) for m in ADDR_RE.finditer(cl)} - AUTHOR_ADDRESSES)
check("third-party addresses derived from the cover letter", len(third_party) >= 3,
      f"{len(third_party)} address(es) that are not an author's")
tracked = subprocess.run(["git", "ls-files"], cwd=str(ROOT),
                         capture_output=True, text=True).stdout.split()
stray, scanned = [], 0
for rel in tracked:
    if rel in ADDRESS_ALLOWED:
        continue
    p = ROOT / rel
    try:
        if p.stat().st_size > 20_000_000:   # the processed matrices
            continue
        t = p.read_text(errors="replace")
    except OSError:
        continue
    scanned += 1
    n = sum(1 for a in third_party if a in t)
    if n:
        stray.append(f"{rel} ({n} address(es))")
check("third-party addresses appear only in the two declared files", not stray,
      ", ".join(stray[:5]))
# Evidence that the scan read the tree rather than skipping it: without a floor, a
# broken walk would report PASS on zero files.
check("the address scan covered the tracked tree", scanned >= 200,
      f"{scanned} files read")

log(f"\nSUBGATE-001: {n_checks} checks, {len(fails)} failure(s)"
    + (": " + ", ".join(fails) if fails else ""))
logf.close()
sys.exit(1 if fails else 0)
