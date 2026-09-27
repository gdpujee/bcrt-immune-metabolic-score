"""Renumber citations in order of first appearance and add 9 verified references.

External review P1-13: the manuscript cited only 10 references, which is thin for
a mature biomarker field.  The added references are the primary PAM50 source, the
METABRIC cohort papers and cBioPortal, the Prosigna analytical validation, the
original REMARK statement, and two external-validation/transportability
methodology papers.  Every entry was verified first-hand against Europe PMC
(title / journal / year / volume / pages / DOI) on 2026-09-26.

The script also enforces the two house rules that the review's P0-1 raised:
  * both manuscript files carry the SAME citation-key sequence, and
  * the printed numbering is order-of-first-appearance.

Outputs: rewrites submission_bcrt/manuscript_submission.md,
         manuscript/manuscript.md, submission_bcrt/references_numbered.md,
         logs/citations.log
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results/raw"
LOGS = ROOT / "logs"
logf = open(LOGS / "citations.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


REFS = {
    "COLLINS2014": "Collins GS, de Groot JA, Dutton S, et al. External validation of multivariable prediction "
                   "models: a systematic review of methodological conduct and reporting. BMC Med Res Methodol. "
                   "2014;14:40. PMID 24645774. DOI https://doi.org/10.1186/1471-2288-14-40.",
    "ALTMAN2000": "Altman DG, Royston P. What do we mean by validating a prognostic model? Stat Med. "
                  "2000;19:453-473. PMID 10694730. DOI "
                  "https://doi.org/10.1002/(sici)1097-0258(20000229)19:4<453::aid-sim350>3.0.co;2-5.",
    "PARKER2009": "Parker JS, Mullins M, Cheang MCU, et al. Supervised risk predictor of breast cancer based on "
                  "intrinsic subtypes. J Clin Oncol. 2009;27:1160-1167. PMID 19204204. DOI "
                  "https://doi.org/10.1200/jco.2008.18.1370.",
    "WIRAPATI": "Wirapati P, Sotiriou C, Kunkel S, et al. Meta-analysis of gene expression profiles in breast "
                "cancer: toward a unified understanding of breast cancer subtyping and prognosis signatures. "
                "Breast Cancer Res. 2008;10:R65. PMID 18662380. DOI https://doi.org/10.1186/bcr2124.",
    "COREPAM": "de Negreiros Botan R, et al. CorePAM: a 24-gene PAM50-derived expression score with "
               "cross-platform external validation for breast cancer prognosis. Breast Cancer Res. "
               "2026;28:142. PMID 42169063. PMCID PMC13474601. DOI "
               "https://doi.org/10.1186/s13058-026-02298-5.",
    "BCRT_TNBC": "Amniouel S, Jafri MS. Development and validation of a highly accurate multigene gene "
                 "expression biomarker to predict chemotherapy response in primary triple-negative breast "
                 "cancer. Breast Cancer Res Treat. 2026;217. PMID 41885969. DOI "
                 "https://doi.org/10.1007/s10549-026-07950-4.",
    "METABSIG": "Wang X, et al. A butyrate metabolism-related gene signature predicts prognosis, immune "
                "landscape, and immunotherapy efficacy in breast cancer. Cancer Med. 2026. PMID 41888914. DOI "
                "https://doi.org/10.1002/cam4.71763.",
    "GSE42568": "Clarke C, Madden SF, Doolan P, et al. Correlating transcriptional networks to breast cancer "
                "survival: a large-scale coexpression analysis. Carcinogenesis. 2013;34:2300-2308. PMID "
                "23740839. DOI https://doi.org/10.1093/carcin/bgt208.",
    "GSE20685": "Kao KJ, Chang KM, Hsu HC, Huang AT. Correlation of microarray-based breast cancer molecular "
                "subtypes and clinical outcomes: implications for treatment optimization. BMC Cancer. "
                "2011;11:143. PMID 21501481. DOI https://doi.org/10.1186/1471-2407-11-143.",
    "GSE45827": "Gruosso T, Mieulet V, Cardon M, et al. Chronic oxidative stress promotes H2AX protein "
                "degradation and enhances chemosensitivity in breast cancer patients. EMBO Mol Med. "
                "2016;8:527-549. PMID 27006338. DOI https://doi.org/10.15252/emmm.201505891.",
    "SCANB": "Brueffer C, Vallon-Christersson J, Grabau D, et al. Clinical value of RNA sequencing-based "
             "classifiers for prediction of the five conventional breast cancer biomarkers: a report from "
             "the population-based multicenter Sweden Cancerome Analysis Network-Breast Initiative. JCO "
             "Precis Oncol. 2018;2:PO.17.00135. PMID 32913985. DOI https://doi.org/10.1200/po.17.00135.",
    "CURTIS2012": "Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 "
                  "breast tumours reveals novel subgroups. Nature. 2012;486:346-352. PMID 22522925. DOI "
                  "https://doi.org/10.1038/nature10983.",
    "PEREIRA2016": "Pereira B, Chin SF, Rueda OM, et al. The somatic mutation profiles of 2,433 breast cancers "
                   "refine their genomic and transcriptomic landscapes. Nat Commun. 2016;7:11479. PMID 27161491. "
                   "DOI https://doi.org/10.1038/ncomms11479.",
    "GAO2013": "Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex cancer genomics and "
               "clinical profiles using the cBioPortal. Sci Signal. 2013;6:pl1. PMID 23550210. DOI "
               "https://doi.org/10.1126/scisignal.2004088.",
    "KEGG": "Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res. "
            "2000;28:27-30. PMID 10592173. DOI https://doi.org/10.1093/nar/28.1.27.",
    "NIELSEN2014": "Nielsen T, Wallden B, Schaper C, et al. Analytical validation of the PAM50-based Prosigna "
                   "Breast Cancer Prognostic Gene Signature Assay and nCounter Analysis System using "
                   "formalin-fixed paraffin-embedded breast tumor specimens. BMC Cancer. 2014;14:177. PMID "
                   "24625003. DOI https://doi.org/10.1186/1471-2407-14-177.",
    "MCSHANE2005": "McShane LM, Altman DG, Sauerbrei W, et al. REporting recommendations for tumour MARKer "
                   "prognostic studies (REMARK). Eur J Cancer. 2005;41:1690-1696. PMID 16043346. DOI "
                   "https://doi.org/10.1016/j.ejca.2005.03.032.",
    "REMARK2023": "Hayes DF, Sauerbrei W, McShane LM. REMARK guidelines for tumour biomarker study reporting: a "
                  "remarkable history. Br J Cancer. 2023;128:443-445. PMID 36476656. DOI "
                  "https://doi.org/10.1038/s41416-022-02046-4.",
    "JUSTICE1999": "Justice AC, Covinsky KE, Berlin JA. Assessing the generalizability of prognostic "
                   "information. Ann Intern Med. 1999;130:515-524. PMID 10075620. DOI "
                   "https://doi.org/10.7326/0003-4819-130-6-199903160-00016.",
}

OLD = {"[1]": ["WIRAPATI"], "[2]": ["COREPAM"], "[3]": ["SCANB"],
       "[4]": ["GSE42568"], "[5]": ["GSE20685"], "[6]": ["GSE45827"],
       "[7]": ["KEGG"], "[8]": ["MCSHANE2005", "REMARK2023"],
       "[9]": ["BCRT_TNBC"], "[10]": ["METABSIG"]}

# (regex anchor, keys appended immediately after the anchor)
INSERTS = [
    (r"Combined models often re-fit in validation or mix endpoints\.",
     ["COLLINS2014", "ALTMAN2000"]),
    (r"varies across PAM50 subtypes", ["PARKER2009"]),
    (r"cBioPortal brca_metabric[;,]\s*1980 analyzable",
     ["CURTIS2012", "PEREIRA2016", "GAO2013"]),
    (r"PAM50 Cox subgroup analysis", ["NIELSEN2014"]),
    (r"single-sample score/threshold portability awaits RNA-seq-platform calibration "
     r"and prospective validation", ["JUSTICE1999"]),
]

FILES = {
    "submission": ROOT / "submission_bcrt" / "manuscript_submission.md",
    "author": ROOT / "manuscript" / "manuscript.md",
}


def split_refs(text):
    """Return (body, refs_heading_start) splitting at the LAST '## References' line."""
    hits = [m.start() for m in re.finditer(r"^## References\s*$", text, re.M)]
    if not hits:
        raise SystemExit("no '## References' heading found")
    cut = hits[-1]
    return text[:cut], text[cut:]


def transform(name, path, fixed_num=None):
    text = path.read_text()
    head, old_refs = split_refs(text)
    start = head.index("## Abstract")
    pre, body = head[:start], head[start:]

    # idempotence guard: the migration is one-shot, so a re-run on an already
    # migrated file is a no-op rather than an error
    if "{{" not in body and len(re.findall(r"^\[\d+\] ", old_refs, re.M)) == len(REFS):
        log(f"[{name}] already migrated ({len(REFS)} references); no-op")
        return None, None, None

    for pat, keys in INSERTS:
        hits = list(re.finditer(pat, body))
        if len(hits) != 1:
            raise SystemExit(f"[{name}] anchor {pat!r} matched {len(hits)} times (need 1)")
        m = hits[0]
        body = body[:m.end()] + "".join("{{%s}}" % k for k in keys) + body[m.end():]

    for old, keys in OLD.items():
        n = body.count(old)
        if n == 0:
            raise SystemExit(f"[{name}] legacy marker {old} not found in body")
        body = body.replace(old, "".join("{{%s}}" % k for k in keys))

    order = []
    for m in re.finditer(r"\{\{([A-Z0-9_]+)\}\}", body):
        if m.group(1) not in order:
            order.append(m.group(1))
    missing = [k for k in REFS if k not in order]
    if missing:
        raise SystemExit(f"[{name}] reference(s) never cited: {missing}")
    # One canonical numbering: derived from the submission document (order of first
    # appearance there) and applied to every file, so the same reference always
    # carries the same number whichever version is opened.
    num = fixed_num or {k: i + 1 for i, k in enumerate(order)}

    def render(m):
        return "[%d]" % num[m.group(1)]

    body = re.sub(r"\{\{([A-Z0-9_]+)\}\}", render, body)

    # collapse adjacent markers [12][13][14] -> [12-14]
    def collapse(m):
        nums = [int(x) for x in m.group(0).strip("[]").split("][")]
        return "[" + (f"{nums[0]}-{nums[-1]}" if len(nums) >= 3 else
                      ",".join(str(n) for n in nums)) + "]"

    body = re.sub(r"(?:\[\d+\])(?:\[\d+\])+", collapse, body)

    # typographic normalisation: a space before the marker, and the marker before
    # the sentence-final period rather than after it
    body = re.sub(r"([A-Za-z])\[(\d)", r"\1 [\2", body)
    body = re.sub(r"\.(\[\d[\d,\-]*\])", r" \1.", body)
    body = body.replace("  ", " ")

    refs_block = "## References\n\n" + "\n".join(
        "[%d] %s" % (num[k], REFS[k]) for k in order) + "\n"
    path.write_text(pre + body.rstrip("\n") + "\n\n" + refs_block)
    log(f"[{name}] {len(order)} references, order: {order}")
    return order, num, body


NUMFILE = RAW / "citation_map.json"


def derive_numbering(path):
    """Recover key -> number from a manuscript's own reference list.

    The numbering migration is one-shot, so on later runs the canonical map has to
    come from the document it was written into rather than from a stored copy that
    could itself drift.  Each printed entry is matched back to exactly one REFS
    value, which also proves the list on disk is still the one REFS describes.
    """
    found = dict(re.findall(r"(?m)^\[(\d+)\] (.+)$", path.read_text()))
    num = {}
    for n, text in found.items():
        matches = [k for k, v in REFS.items() if text.startswith(v[:60])]
        if len(matches) != 1:
            raise SystemExit(f"{path.name}: [{n}] maps to {len(matches)} references "
                             f"{matches}: {text[:70]!r}")
        num[matches[0]] = int(n)
    if sorted(num.values()) != list(range(1, len(REFS) + 1)):
        raise SystemExit(f"{path.name}: recovered numbering is not 1..{len(REFS)}: "
                         f"{sorted(num.values())}")
    return num


def save_numbering(order, num):
    RAW.mkdir(parents=True, exist_ok=True)
    NUMFILE.write_text(json.dumps({"numbering": {k: num[k] for k in order}},
                                  indent=1) + "\n")
    log(f"WROTE {NUMFILE.relative_to(ROOT)} ({len(order)} entries)")


def rewrite_refs_block(path, order, num):
    """Rewrite `## References` from the canonical strings in REFS.

    Makes the printed list a derived artifact in both carriers, so a removed
    editorial note or a corrected DOI cannot survive in one file only.  Whatever
    section follows the references is preserved — the submission manuscript now
    ends with Statements and Declarations, which must stay after the references.
    """
    text = path.read_text()
    hits = [m.start() for m in re.finditer(r"(?m)^## References\s*$", text)]
    if not hits:
        raise SystemExit(f"{path.name}: no '## References' heading")
    cut = hits[-1]
    nxt = re.search(r"(?m)^## ", text[cut + 1:])
    tail = text[cut + 1 + nxt.start():] if nxt else ""
    block = ("## References\n\n"
             + "\n".join("[%d] %s" % (num[k], REFS[k]) for k in order) + "\n")
    new = text[:cut] + block + ("\n" + tail if tail else "")
    if new == text:
        log(f"[{path.name}] references block already current ({len(order)} entries)")
        return
    path.write_text(new)
    log(f"[{path.name}] references block rewritten ({len(order)} entries, "
        f"{len(text)} -> {len(new)} chars)")


def sentence_for(body, marker):
    """The sentence (or clause) containing `marker`, for the internal map."""
    for m in re.finditer(r"[^.\n]*" + re.escape(marker) + r"[^.\n]*", body):
        s = " ".join(m.group(0).split())
        if s:
            return s
    return ""


results = {}
results["submission"] = transform("submission", FILES["submission"])
migrated = results["submission"][0] is not None

if migrated:
    sub_order, sub_num, _ = results["submission"]
else:
    log("CITATION-SYNC-001: already migrated; numbering recovered from the manuscript")
    sub_num = derive_numbering(FILES["submission"])
    sub_order = [k for k, _ in sorted(sub_num.items(), key=lambda kv: kv[1])]

if NUMFILE.exists():
    persisted = json.load(open(NUMFILE))["numbering"]
    if persisted != sub_num:
        raise SystemExit("citation_map.json disagrees with the manuscript's numbering:\n"
                         f"  persisted:  {persisted}\n  manuscript: {sub_num}")
    log(f"citation_map.json agrees with the manuscript ({len(sub_num)} entries)")
save_numbering(sub_order, sub_num)

if migrated:
    results["author"] = transform("author", FILES["author"], fixed_num=sub_num)
    seqs = {name: tuple(r[0]) for name, r in results.items()}
    if set(seqs["submission"]) != set(seqs["author"]):
        raise SystemExit("the two manuscript files do not cite the same reference set: "
                         f"{seqs}")
    for name, (order, num, _b) in results.items():
        if num != sub_num:
            raise SystemExit(f"[{name}] numbering diverges from the canonical map")
    log(f"CITATION-SYNC-001: both files cite the same {len(sub_order)} references under "
        "one canonical numbering (order of first appearance in the submission manuscript)")

# Both carriers get the reference block re-emitted, so the printed list is derived.
for path in FILES.values():
    rewrite_refs_block(path, sub_order, sub_num)

# Reader-facing numbered list.
out = ROOT / "submission_bcrt" / "references_numbered.md"
out.write_text(
    "# Numbered references (order of first appearance; each verified first-hand "
    "against Europe PMC)\n\n"
    + "\n".join("[%d] %s" % (sub_num[k], REFS[k]) for k in sub_order) + "\n")
log(f"WROTE {out.relative_to(ROOT)} with {len(sub_order)} entries")

# Internal navigation aid.  It used to be appended to the file above, where its
# truncated fragments and wrong anchors shipped inside the submission package.
sub_body = FILES["submission"].read_text()
sub_body = sub_body[:sub_body.rfind("## References")]
review = ROOT / "review" / "CITATION_MAP.md"
review.write_text(
    "# Citation map (internal working note; not part of the submission package)\n\n"
    "Reference number -> the sentence in `submission_bcrt/manuscript_submission.md` "
    "that first cites it.\n\n"
    + "\n".join(f"- [{sub_num[k]}] {sentence_for(sub_body, '[%d]' % sub_num[k])}"
                for k in sub_order) + "\n")
log(f"WROTE {review.relative_to(ROOT)} ({review.stat().st_size} bytes)")
logf.close()
