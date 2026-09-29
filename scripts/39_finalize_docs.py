"""Inject the new evidence into every manuscript document, once, with assertions.

Covers the external-review items that are text-level:
  P1-6/P1-12  METABRIC clinical-adjusted Cox + missingness transparency
  P1-7        formal time-varying (non-PH) model + HR(t) figure
  P1-9        absolute cutoff/calibration framed as a secondary analysis
  P1-11       derivation cohort framed as construction, not inference
  P1-10       transportability table pointer
  P1-14       supplement + REMARK checklist physically present in the package

Every edit is an (anchor, replacement) pair asserted to match exactly once, so a
stale anchor fails loudly instead of silently dropping a change.

Outputs: rewrites submission_bcrt/manuscript_submission.md, manuscript/manuscript.md,
         regenerates submission_bcrt/legends.md from the manuscript,
         submission_bcrt/title_page.md, submission_bcrt/declarations.md,
         manuscript/supplement.md and its package copy submission_bcrt/supplement.md,
         submission_bcrt/REMARK_checklist.md, review/PACKAGE_INDEX.md,
         logs/finalize_docs.log

The supplement is a reader-facing document, so its three generated regions are
built here from the raw result files: S3 (Table S1), S6 (supporting evidence
summary) and S7 (figure plates).  Nothing in them names an internal script, run
id, hash or result path — external review P0-5/P0-6.
"""
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / "submission_bcrt"
MS = ROOT / "manuscript"
RAW = ROOT / "results/raw"
TAB = ROOT / "tables"
LOGS = ROOT / "logs"
logf = open(LOGS / "finalize_docs.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


tvc = json.load(open(RAW / "ph_timevarying.json"))
madj = json.load(open(RAW / "metabric_adjusted.json"))
aps = json.load(open(RAW / "adjusted_per_sd.json"))
phex = json.load(open(RAW / "ph_diagnostics_exact.json"))
splits = json.load(open(RAW / "cutoff_splits.json"))

cc = madj["complete_case"]
ccr = cc["terms"]["risk_sd"]
miss = madj["missingness"]
mi = madj["imputed_sensitivity"]["terms"]["risk_sd"]
ci_inc = madj["complete_case_incremental"]
red = madj["reduced_core"]
red_inc = madj["reduced_core_incremental"]


def hr(term):
    return f"{term['HR']:.2f} [{term['CI'][0]:.2f},{term['CI'][1]:.2f}], p={term['p']:.1e}"


def miss_txt():
    return ", ".join(f"{k} {v['missing']}/{v['n_total']} ({v['pct']}%)"
                     for k, v in miss.items())


def inc_txt(inc, prefix="ΔC"):
    """Incremental-value phrase, with the bootstrap interval read rather than
    paraphrased: an interval that spans zero must not be described as an
    improvement, however small the point estimate looks."""
    lo, hi = inc["delta_C_CI95_boot500"]
    spans0 = lo <= 0 <= hi
    fmt_bound = lambda x: f"{x:+.4f}" if abs(x) < 0.001 else f"{x:+.3f}"
    base = (f"{prefix} {inc['delta_C']:+.3f} "
            f"[{fmt_bound(lo)},{fmt_bound(hi)}] "
            f"(bootstrap 500, seed 42)")
    if spans0:
        base += " — the interval includes zero, so the added discrimination is not " \
                "distinguishable from none"
    return base


JOBS = []
PROBE_LEN = 80


def probe_of(old, new):
    """A distinctive slice of the injected text used to detect a re-run.

    Taken from just after the longest common prefix of anchor and replacement, so
    it works for a replacement that quotes the anchor verbatim *and* for one that
    rewrites the anchor's tail (several pairs change a trailing '.' to ';' and
    append).  When that leaves too few characters to be distinctive — the
    declarations pair only differs in its last digit — the tail of the replacement
    is used instead.
    """
    i = 0
    while i < min(len(old), len(new)) and old[i] == new[i]:
        i += 1
    probe = new[i:i + PROBE_LEN].strip()
    if len(probe) < 24:
        probe = new[-PROBE_LEN:].strip()
    return probe


def classify(t, old, new):
    """Return 'apply', 'already' or 'ambiguous' for one (old, new) pair.

    Re-running an injector is the classic way to silently duplicate prose, so the
    state is decided explicitly rather than by "the anchor matched, replace it".
    Anything that is neither cleanly un-applied nor cleanly applied aborts.
    """
    probe = probe_of(old, new)
    if not probe:
        raise SystemExit(f"empty idempotence probe for {old[:60]!r}")
    n_old, n_probe = t.count(old), t.count(probe)
    if n_probe == 0:
        return "apply" if n_old == 1 else "ambiguous"
    if n_old == 0:
        return "already"
    if old in new:
        # the replacement quotes the anchor verbatim (some prepend to it, some append
        # to it), so the anchor necessarily survives inside the injected text — that
        # is the applied state, not drift
        return "already"
    return "ambiguous"


def stage(path, pairs, label):
    """Queue one edit batch for `path`.

    Nothing is written until every anchor of every queued batch has been resolved,
    so a stale anchor anywhere leaves the working tree untouched instead of
    half-rewritten — a partial injection is much harder to notice than a loud abort.
    """
    t = path.read_text()
    todo = []
    for old, new in pairs:
        verdict = classify(t, old, new)
        if verdict == "ambiguous":
            raise SystemExit(f"[{label}] anchor is neither absent-with-new-text-present "
                             f"nor present-once (old x{t.count(old)}, "
                             f"new-text x{t.count(probe_of(old, new))}): {old[:90]!r}")
        if verdict == "apply":
            todo.append((old, new))
    if todo:
        JOBS.append((path, todo, label))
    else:
        log(f"[{label}] already applied — nothing to do ({len(pairs)} pairs)")


def commit_jobs():
    for path, pairs, label in JOBS:
        t = path.read_text()
        for old, new in pairs:
            n = t.count(old)
            if n != 1:
                raise SystemExit(
                    f"[{label}] anchor matched {n}x at write time (need 1): {old[:90]!r}")
            t = t.replace(old, new)
        path.write_text(t)
        log(f"[{label}] applied {len(pairs)} edits")


def split_sections(text):
    """Split a markdown document into (preamble, [(heading, body), ...]).

    `^## ` (two hashes + space) is the section marker; `### ` subsections stay
    inside their parent's body because the third character is not a space.
    """
    parts = re.split(r"(?m)^(## .*)$", text)
    return parts[0], [(parts[i].strip(), parts[i + 1]) for i in range(1, len(parts), 2)]


def move_section_after(path, section, target):
    """Move `## section` so it directly follows `## target`, asserting the result.

    BCRT: the 'Statements and Declarations' heading "should be placed after the
    References section", and a submission that omits required statements is
    returned as incomplete.  The section is moved as a block rather than retyped,
    so its text keeps exactly one owner, and the new order is re-read from disk
    after the write — asserting on the in-memory copy would not prove anything.
    """
    pre, secs = split_sections(path.read_text())
    names = [h[3:].strip() for h, _ in secs]
    if names.count(section) != 1 or names.count(target) != 1:
        raise SystemExit(f"[order] {path.name}: need exactly one '{section}' and one "
                         f"'{target}' section, found {names.count(section)}/"
                         f"{names.count(target)}")
    if names.index(section) == names.index(target) + 1:
        log(f"[order] '{section}' already directly after '{target}' — nothing to do")
        return
    item = secs.pop(names.index(section))
    names.remove(section)
    secs.insert(names.index(target) + 1, item)
    path.write_text(pre + "".join(h + b for h, b in secs))
    _, chk = split_sections(path.read_text())
    order = [h[3:].strip() for h, _ in chk]
    if order.index(section) != order.index(target) + 1:
        raise SystemExit(f"[order] {path.name}: '{section}' did not land after "
                         f"'{target}' (order now {order})")
    log(f"[order] '{section}' moved to directly after '{target}'")


# Fail before touching anything if the supplementary checklist is absent: it is
# a required part of the submission package (P1-14).
if not (MS / "REMARK_checklist.md").exists():
    raise SystemExit("REMARK_checklist.md missing from manuscript/ — cannot ship the package")


# ------------------------------------------------------------------ submission md
sub = SUB / "manuscript_submission.md"
pairs = []

# --- Methods: METABRIC adjusted model + missingness ---
pairs.append((
    "covariate-adjusted Cox (GSE20685: age+T+N; SCAN-B: age+ER+HER2)",
    "covariate-adjusted Cox (GSE20685: age+T+N; SCAN-B: age+ER+HER2; METABRIC: "
    "age+positive lymph nodes+grade+ER+HER2+tumor-size category derived from recorded "
    "millimetres)"))

pairs.append((
    "for SCAN-B, age had zero missing values, ER status was missing in 200 patients, "
    "and HER2 status was missing in 122 patients (complete-case n=2,963, 286 deaths).",
    "for SCAN-B, age had zero missing values, ER status was missing in 200 patients, "
    "and HER2 status was missing in 122 patients (complete-case n=2,963, 286 deaths); "
    f"for METABRIC, per-variable missingness was {miss_txt()} "
    f"(complete-case n={cc['n']}, {cc['deaths']} deaths)."))

# --- Methods: formal time-varying model ---
pairs.append((
    "whereas SCAN-B risk association was stable over time (Fig. S3).",
    "whereas SCAN-B risk association was stable over time (Fig. S3).  A formal "
    "time-varying-coefficient model, h(t | risk) = h0(t)·exp(β1·risk + β2·risk·log t), "
    "was fitted to each validation cohort (`scripts/35_ph_timevarying.R`; "
    "`results/raw/ph_timevarying.json`), giving HR(t) = exp(β1 + β2·log t) with a full "
    "covariance matrix; the reported likelihood-ratio test compares that model with the "
    "proportional-hazards model, and HR(t) curves with 95% confidence bands are shown "
    "in Fig. S4.  Absolute threshold behaviour is a secondary, transportability-focused "
    "analysis throughout: the primary estimand is the continuous per-SD association."))

# --- Results: METABRIC adjusted Cox ---
_sub_now = sub.read_text()
_metabric_adjustment_markers = (
    "Clinical adjustment of METABRIC (age+positive nodes+grade+ER+HER2+tumor-size category",
    "clinical-only C 0.663 to clinical+score C 0.666",
    "median/mode-imputation sensitivity gave a similar adjusted estimate",
)
if all(marker in _sub_now for marker in _metabric_adjustment_markers):
    log("[results-metabric-adjusted] current detailed analysis already present")
elif not any(marker in _sub_now for marker in _metabric_adjustment_markers):
    pairs.append((
        "so PAM50 heterogeneity remains a SCAN-B-only exploratory finding.",
        "so PAM50 heterogeneity remains a SCAN-B-only exploratory finding. Clinical "
        f"adjustment of METABRIC (age+positive nodes+grade+ER+HER2+tumor size; "
        f"complete-case n={cc['n']}, {cc['deaths']} deaths) attenuated but did not remove the "
        f"association: adjusted HR {hr(ccr)} per SD. Discrimination barely moved "
        f"(clinical-only C {ci_inc['C_clinical_only']:.3f} to clinical+score C "
        f"{ci_inc['C_clinical_plus_risk']:.3f}; {inc_txt(ci_inc)}; likelihood-ratio "
        f"p={ci_inc['LRT_p']:.3g}), and the median/mode-imputation "
        f"sensitivity gave a similar adjusted estimate (HR {mi['HR']:.2f} "
        f"[{mi['CI'][0]:.2f},{mi['CI'][1]:.2f}])."))
else:
    raise SystemExit("[results-metabric-adjusted] current analysis is partial; review it manually")

# --- Results: formal time-varying results ---
pairs.append((
    "while SCAN-B risk association remained approximately constant over follow-up (Fig. S3).",
    "while SCAN-B risk association remained approximately constant over follow-up "
    "(Fig. S3). The formal time-varying-coefficient model agreed (Fig. S4): the "
    f"score×log(time) term was significant in GSE20685 (β2="
    f"{tvc['GSE20685']['beta_tt']:.3f}, LRT χ²={tvc['GSE20685']['LRT_chisq']:.2f}, df=1, "
    f"p={tvc['GSE20685']['LRT_p']:.3g}) and in METABRIC (β2={tvc['METABRIC']['beta_tt']:.3f}, "
    f"LRT χ²={tvc['METABRIC']['LRT_chisq']:.2f}, p={tvc['METABRIC']['LRT_p']:.1e}) but not "
    f"in SCAN-B (β2={tvc['SCANB']['beta_tt']:.3f}, LRT χ²={tvc['SCANB']['LRT_chisq']:.2f}, "
    f"p={tvc['SCANB']['LRT_p']:.2f}). The implied hazard ratio per SD fell from "
    f"{tvc['GSE20685']['hr_at']['t1']:.2f} "
    f"[{tvc['GSE20685']['ci_lo']['t1']:.2f}-{tvc['GSE20685']['ci_hi']['t1']:.2f}] at 1 year to "
    f"{tvc['GSE20685']['hr_at']['t10']:.2f} "
    f"[{tvc['GSE20685']['ci_lo']['t10']:.2f}-{tvc['GSE20685']['ci_hi']['t10']:.2f}] at 10 years "
    f"in GSE20685 and from {tvc['METABRIC']['hr_at']['t1']:.2f} "
    f"[{tvc['METABRIC']['ci_lo']['t1']:.2f}-{tvc['METABRIC']['ci_hi']['t1']:.2f}] to "
    f"{tvc['METABRIC']['hr_at']['t10']:.2f} "
    f"[{tvc['METABRIC']['ci_lo']['t10']:.2f}-{tvc['METABRIC']['ci_hi']['t10']:.2f}] in METABRIC, "
    f"versus only {tvc['SCANB']['hr_at']['t1']:.2f} to {tvc['SCANB']['hr_at']['t10']:.2f} in "
    "SCAN-B; the locked score therefore carries a reproducible but time-dependent "
    "prognostic association that is strongest during the early years of follow-up."))

# --- Results: incremental value on a common scale ---
pairs.append((
    "Incremental discrimination over clinical-only models was modest in exploratory "
    "same-cohort fits (ΔC +0.028 [0.003,0.072] and +0.024 [0.011,0.039]; Table 4);",
    "Incremental discrimination over clinical-only models was modest in exploratory "
    "same-cohort fits (ΔC +0.028 [0.003,0.072] and +0.024 [0.011,0.039]; Table 4); "
    f"in METABRIC the corresponding complete-case increment was {inc_txt(ci_inc)} "
    f"(likelihood-ratio p={ci_inc['LRT_p']:.3g}), so the score adds little measurable "
    "discrimination once routine clinical variables are accounted for and no "
    "clinical-utility claim is made;"))

# --- Abstract: the clinical-adjusted increment, so the abstract cannot imply utility ---
# --- Table legends: Table 3 is now the transportability summary; Table S1 exists ---
pairs.append((
    "- **Table 3.** Overall survival performance of the locked 14-gene score across "
    "cohorts (C-index, continuous per-SD HR, and binary HR with 95% confidence intervals "
    "and p-values).",
    "- **Table 3.** Evaluation summary of the recorded 14-gene score across GSE20685 "
    "and the subsequent SCAN-B and METABRIC cohorts: one row per model-evaluation "
    "cohort with n, deaths, unadjusted and clinically adjusted hazard ratio per cohort "
    "SD, C-index, raw-score Cox coefficient, locked-cutoff high/low split, binary hazard ratio, "
    "and the proportional-hazards test p-value for the score. All hazard ratios are on "
    "the per-SD scale so that adjusted and unadjusted estimates are directly comparable; "
    "per-1-score-unit estimates and the non-inferential derivation-cohort row are in "
    "Table S1 of Online Resource 1. Primary adjusted models use median imputation for "
    "GSE20685 (n=327) and SCAN-B (n=3,273), and complete-case analysis for METABRIC "
    "(n=1,815; unrounded 95% CI [1.004, 1.168], p=0.040; full-cohort imputation "
    "sensitivity n=1,980, HR 1.08 [1.01, 1.16]).\n"
    "- **Table S1.** Derivation-cohort estimates (non-inferential) and per-1-score-unit "
    "adjusted estimates (Online Resource 1)."))

apply_label = "submission"
stage(sub, pairs, apply_label)

# ------------------------------------------------------------------ author md
auth = MS / "manuscript.md"
pairs = [
    ("PAM50 interaction p=0.42 (SCAN-B heterogeneity not replicated).",
     "PAM50 interaction p=0.42 (SCAN-B heterogeneity not replicated). Clinical "
     f"adjustment (six-covariate set; complete-case n={cc['n']}, "
     f"{cc['deaths']} deaths; missingness {miss_txt()}) gave adjusted HR {hr(ccr)} per SD; "
     f"clinical-only C {ci_inc['C_clinical_only']:.3f} to clinical+score C "
     f"{ci_inc['C_clinical_plus_risk']:.3f} ({inc_txt(ci_inc)}, "
     f"likelihood-ratio p={ci_inc['LRT_p']:.3g}); imputation sensitivity HR {mi['HR']:.2f} "
     f"[{mi['CI'][0]:.2f},{mi['CI'][1]:.2f}]."),
    ("while SCAN-B risk association was stable over time (≤5y HR 1.447 vs >5y 1.415; "
     "time-stratified Cox).",
     "while SCAN-B risk association was stable over time (≤5y HR 1.447 vs >5y 1.415; "
     "time-stratified Cox). Formal time-varying-coefficient Cox "
     "(`scripts/35_ph_timevarying.R`; HR(t) curves in Fig. S4) agrees: score×log(time) "
     f"LRT χ²={tvc['GSE20685']['LRT_chisq']:.2f} p={tvc['GSE20685']['LRT_p']:.3g} "
     f"(GSE20685), χ²={tvc['METABRIC']['LRT_chisq']:.2f} p={tvc['METABRIC']['LRT_p']:.1e} "
     f"(METABRIC), χ²={tvc['SCANB']['LRT_chisq']:.2f} p={tvc['SCANB']['LRT_p']:.2f} "
     f"(SCAN-B); HR(t) per SD falls {tvc['GSE20685']['hr_at']['t1']:.2f}->"
     f"{tvc['GSE20685']['hr_at']['t10']:.2f} (GSE20685), {tvc['METABRIC']['hr_at']['t1']:.2f}->"
     f"{tvc['METABRIC']['hr_at']['t10']:.2f} (METABRIC), {tvc['SCANB']['hr_at']['t1']:.2f}->"
     f"{tvc['SCANB']['hr_at']['t10']:.2f} (SCAN-B)."),
    ("single-sample score/threshold portability awaits RNA-seq-platform calibration and "
     "prospective validation [19];",
     "the derivation cohort served for construction of the frozen score rather than as "
     "inferential evidence, so prognostic inference rests on the external validations; "
     "absolute threshold transport is a secondary analysis and single-sample "
     "score/threshold portability awaits RNA-seq-platform calibration and prospective "
     "validation [19];"),
]
stage(auth, pairs, "author")

# ------------------------------------------------------------------ legends.md
# Generated from the manuscript, so the caption text a reader sees in this file is
# the same text the manuscript carries and cannot drift away from it.  The
# generator/source annotations that used to be appended to each entry were internal
# engineering traces in a file that ships with the submission (external review
# P0-6); the provenance map now lives in review/FIGURE_PROVENANCE.md.
leg = SUB / "legends.md"

FIG_FILES = {
    "Fig. 1": "Fig1_flow", "Fig. 2": "Fig2_volcano_42568",
    "Fig. 3": "Fig3_coef_forest_train", "Fig. 4": "train_KM",
    "Fig. 5": "Fig3_ROC_train", "Fig. 6": "valid_KM",
    "Fig. 7": "Fig4_ROC_valid", "Fig. 8": "Fig5_pathway_valid",
    "Fig. 9": "Fig5_checkpoints_valid", "Fig. 10": "RNASEQ_KM",
    "Fig. 11": "Fig_PAM50_forest", "Fig. 12": "Fig12_crosscohort_forest",
    "Fig. S1": "QC_followup_expr", "Fig. S2": "risk_distributions",
    "Fig. S3": "PH_schoenfeld", "Fig. S4": "HR_t_curves",
}
TBL_FILES = {
    "Table 1": "Tab1_cohorts_v3", "Table 2": "Tab2_coefficients",
    "Table 3": "Tab3_performance_v3", "Table 4": "Tab4_incremental_value",
    "Table S1": "TabS1_derivation_and_units", "Table S2": None, "Table S3": None,
}
# The order a reader expects: main figures 1..12, then supplementary S1..S4.
EXPECTED_ORDER = ([f"Fig. {i}" for i in range(1, 13)]
                  + [f"Fig. S{i}" for i in range(1, 5)])


def build_legends():
    lines = sub.read_text().split("\n")

    def bullets(heading):
        out, on = [], False
        for ln in lines:
            if ln.strip() == heading:
                on = True
                continue
            if on and ln.startswith("## "):
                break
            if on and ln.startswith("- **"):
                out.append(ln.strip())
        return out

    figs, tbls = bullets("## Figure legends"), bullets("## Tables")
    got = [re.match(r"- \*\*(Fig\.? ?S?\d+)\*\*", b).group(1).replace("Fig ", "Fig. ")
           for b in figs]
    if got != EXPECTED_ORDER:
        raise SystemExit(f"figure legends are out of order or incomplete: {got}")
    got_t = [re.match(r"- \*\*(Table ?S?\d+)\.\*\*", b).group(1) for b in tbls]
    if got_t != list(TBL_FILES):
        raise SystemExit(f"table titles are out of order or incomplete: {got_t}")

    L = ["# Figure legends and table titles", "",
         "This file repeats the legend text of `manuscript_submission.md` beside the "
         "file name of each figure, so artwork can be uploaded against its caption. "
         "Figures S1-S4 and Tables S1-S3 are included in Online Resource 1.", "",
         "## Figures", ""]
    for b in figs:
        m = re.match(r"- \*\*(Fig\.? ?S?\d+)\*\* (.*)", b)
        label = m.group(1).replace("Fig ", "Fig. ")
        stem = FIG_FILES[label]
        L.append(f"- **{label}** {m.group(2)}  ")
        L.append(f"  Files: `figures/{stem}.pdf` (vector), `figures/{stem}.png` (raster).")
    L += ["", "## Tables", ""]
    for b in tbls:
        m = re.match(r"- \*\*(Table ?S?\d+)\.\*\* (.*)", b)
        label = m.group(1)
        L.append(f"- **{label}.** {m.group(2)}  ")
        if TBL_FILES.get(label):
            L.append(f"  File: `tables/{TBL_FILES[label]}.tsv`.")
    L.append("")
    leg.write_text("\n".join(L))
    log(f"[legends] regenerated from the manuscript ({len(figs)} figures, "
        f"{len(tbls)} tables; {leg.stat().st_size} bytes)")

# ------------------------------------------------------------------ keywords (both sources)
# BCRT: "Please provide 4 to 6 keywords which can be used for indexing purposes."
# The previous revision carried nine, which is out of range.  Six are kept, and
# they are chosen to still cover the design (external validation, transportability)
# and the subtype analysis (PAM50) without listing individual cohort names —
# SCAN-B and METABRIC are not indexing terms a reader would search.
KW_OLD = ("breast cancer; prognostic signature; overall survival; external validation; "
          "cross-platform transportability; calibration; PAM50; SCAN-B; METABRIC.")
KW_NEW = ("breast cancer; prognostic biomarker; overall survival; external validation; "
          "PAM50; transportability.")
stage(sub, [("## Keywords\n" + KW_OLD, "## Keywords\n" + KW_NEW)], "submission-keywords")

tp = SUB / "title_page.md"
_ms_title = re.search(r"(?m)^## Title\s*\n([^\n]+)", sub.read_text())
if not _ms_title:
    raise SystemExit("submission manuscript has no ## Title section")
TITLE_OLD = "Cross-platform external validation of an immune–metabolic transcriptional score for breast cancer overall survival"
TITLE_NEW = _ms_title.group(1).strip()
stage(tp, [("Title: " + TITLE_OLD, "Title: " + TITLE_NEW)], "title-page-title")
stage(tp, [("Keywords: " + KW_OLD, "Keywords: " + KW_NEW)], "title-page-keywords")

# ------------------------------------------------------------------ declarations
dec = SUB / "declarations.md"
if "corrected conditional-alpha bootstrap" not in dec.read_text():
    raise SystemExit("declarations omit the archived bootstrap analysis")
log("[declarations] current archive contents described")

# ------------------------------------------------------------------ archived release
# The paper promises that the code and the derived tables are available.  The
# release is cited by its CONCEPT DOI, which always resolves to the most recent
# version; a version DOI is minted per release, so citing one would send a reader
# to the superseded v3.4.2 archive - the snapshot whose contents no longer matched
# the manuscript.  Anchors are per-sentence because `declarations.md` states
# availability twice (data and code) with the same trailing phrase.
REPO = "https://github.com/gdpujee/bcrt-immune-metabolic-score"
# Version DOI pins the archived snapshot that matches this manuscript. It is
# minted by Zenodo after the GitHub release tag exists, then written to
# CITATION.cff before the submission documents are rebuilt.
_cff_text = (ROOT / "CITATION.cff").read_text()
_version_doi = re.search(r'(?m)^doi:\s*"?(10\.5281/zenodo\.\d+)', _cff_text)
if not _version_doi:
    raise SystemExit("CITATION.cff has no minted version DOI; finalize the release first")
ARCHIVE = f"https://doi.org/{_version_doi.group(1)}"
CONCEPT = "https://doi.org/10.5281/zenodo.22994650"
AVAIL_TAIL = f"publicly available at {REPO}."

VERSION_DOI = _version_doi.group(1)
CONCEPT_DOI = re.search(r"10\.5281/zenodo\.\d+", CONCEPT).group(0)


def sync_version_doi(path):
    """Update only this active submission source's Zenodo version DOI."""
    text = path.read_text()
    found = re.findall(r"10\.5281/zenodo\.\d+", text)
    replaced = 0

    def replace(match):
        nonlocal replaced
        if match.group(0) == CONCEPT_DOI:
            return match.group(0)
        replaced += 1
        return VERSION_DOI

    updated = re.sub(r"10\.5281/zenodo\.\d+", replace, text)
    if replaced:
        path.write_text(updated)
    log(f"[version-doi] {path.relative_to(ROOT)}: {replaced} updated; "
        f"{len(found) - replaced} concept reference(s) kept")


for _path in (sub, dec, SUB / "cover_letter.md", MS / "supplement.md"):
    sync_version_doi(_path)

# ------------------------------------------------------------------ supplement + REMARK into the package
# The supplement is a reader-facing document and carries three generator-owned
# regions: S3 (Table S1), S6 (supporting evidence summary) and S7 (figure plates).
# Every number in them is read from the raw result files, so the supplement cannot
# drift away from the analyses the way a hand-written summary would.  The markers
# are the contract; because each builder is a pure function of those files,
# re-running replaces each region in place and is idempotent.
SUP_SRC = MS / "supplement.md"
SUP = SUB / "supplement.md"

_availability = SUP_SRC.read_text()
if (VERSION_DOI not in _availability or CONCEPT_DOI not in _availability
        or "corrected conditional-alpha bootstrap" not in _availability
        or "not yet in the public" in _availability):
    raise SystemExit("supplement availability statement is not synchronized to the public release")
log("[supplement-availability] current Zenodo version and public analyses cited")

SUPP_FIG_FILES = {
    "Fig. S1": "QC_followup_expr", "Fig. S2": "risk_distributions",
    "Fig. S3": "PH_schoenfeld", "Fig. S4": "HR_t_curves",
}
REGIONS = ("S3", "S6", "S7")


def fill_region(text, name, builder, path):
    """Replace one `<!-- NAME-GENERATED-START/END -->` region with built text.

    A missing, duplicated or inverted marker pair aborts rather than passing
    through: a region that is silently skipped renders as nothing at all in the
    supplement PDF, which is the failure mode hardest to notice.
    """
    start, end = f"<!-- {name}-GENERATED-START -->", f"<!-- {name}-GENERATED-END -->"
    ns, ne = text.count(start), text.count(end)
    if ns != 1 or ne != 1:
        raise SystemExit(f"{path.name}: region {name} needs exactly one marker pair "
                         f"(found {ns} start / {ne} end)")
    a, b = text.index(start) + len(start), text.index(end)
    if a > b:
        raise SystemExit(f"{path.name}: region {name} markers are out of order")
    body = builder().strip()
    if not body:
        raise SystemExit(f"{path.name}: region {name} built empty text")
    return text[:a] + "\n" + body + "\n" + text[b:]


def manuscript_legends():
    """Figure legends of the submission manuscript, keyed by label.

    The supplement restates these captions rather than paraphrasing them, so the
    two documents cannot describe the same figure differently.
    """
    ms = (SUB / "manuscript_submission.md").read_text()
    out = {}
    for m in re.finditer(r"(?m)^- \*\*(Fig\.? ?S?\d+)\*\* (.*)$", ms):
        out[m.group(1).replace("Fig ", "Fig. ")] = m.group(2).strip()
    return out


def build_s3():
    """Table S1, rendered from the TSV that ships beside it."""
    rows = [r.split("\t") for r in
            (TAB / "TabS1_derivation_and_units.tsv").read_text().strip().splitlines()]
    head, data = rows[0], rows[1:]
    if head != ["Analysis", "n", "deaths", "Estimate", "Note"]:
        raise SystemExit(f"Table S1 columns changed: {head}")
    if len(data) != 4:
        raise SystemExit(f"Table S1 has {len(data)} data rows (expected 4)")
    L = ["**Table S1.** Derivation-cohort estimates (non-inferential) and "
         "per-1-score-unit adjusted estimates. Validation rows are on the "
         "per-1-score-unit scale except METABRIC, whose model enters "
         "the score as its per-SD linear predictor.", "",
         "| " + " | ".join(head) + " |",
         "|" + "---|" * len(head)]
    L += ["| " + " | ".join(r) + " |" for r in data]
    L += ["",
          "The per-SD equivalents used in the main Table 3 follow by rescaling with "
          "each cohort's score standard deviation (S6.3)."]
    return "\n".join(L)


def build_s6():
    c = madj["complete_case"]
    cr = c["terms"]["risk_sd"]
    inc = madj["complete_case_incremental"]
    rinc = madj["reduced_core_incremental"]
    red = madj["reduced_core"]["terms"]["risk_sd"]
    corr = madj["coding_correction"]
    sup_old = corr["superseded_complete_case"]
    ex = phex
    smp = json.load(open(RAW / "metabric_sample_multiplicity.json"))
    spl = json.load(open(RAW / "cutoff_splits.json"))
    corr_s = json.load(open(RAW / "corrective_summary.json"))
    rna_s = json.load(open(RAW / "rnaseq_summary.json"))
    met_s = json.load(open(RAW / "metabric_summary.json"))

    def hrv(t):
        return f"{t['HR']:.2f} [{t['CI'][0]:.2f},{t['CI'][1]:.2f}], p={t['p']:.3g}"

    L = []
    L.append("### S6.1 METABRIC clinical-adjusted Cox")
    L.append("The METABRIC adjusted model included age, positive-node count, grade, ER "
             "status, HER2 status and tumor-size category (T category). The locked "
             "score enters as the per-SD linear predictor and is never refit. The "
             "materials in this review package do not establish when the covariate "
             "set was selected.")
    L.append("")
    L.append(f"- Complete-case model (PRIMARY), n={c['n']}, {c['deaths']} deaths: "
             f"adjusted HR per SD **{hrv(cr)}** (unrounded 95% CI "
             f"[{cr['CI'][0]:.3f}, {cr['CI'][1]:.3f}]; rounded to 2 decimal places as "
             f"[{cr['CI'][0]:.2f}, {cr['CI'][1]:.2f}] in Table 3).")
    for term in madj["model_covariates"]:
        t = c["terms"][term]
        L.append(f"  - {term}: HR {t['HR']:.3f} [{t['CI'][0]:.3f},{t['CI'][1]:.3f}], "
                 f"p={t['p']:.3g}")
    L.append("- Per-variable missingness in the 1,980-patient analysis frame: "
             + ", ".join(f"{k} {v['missing']}/{v['n_total']} ({v['pct']}%)"
                         for k, v in madj["missingness"].items()) + ".")
    L.append(f"- Incremental discrimination over the same clinical variables: clinical-only "
             f"C {inc['C_clinical_only']:.3f} -> clinical+score C "
             f"{inc['C_clinical_plus_risk']:.3f}; delta C {inc['delta_C']:+.3f} "
             f"[{inc['delta_C_CI95_boot500'][0]:+.3f},{inc['delta_C_CI95_boot500'][1]:+.3f}] "
             f"(bootstrap 500, seed 42), likelihood-ratio chi-square={inc['LRT_chi2']:.2f} "
             f"(df={inc['LRT_df']}), p={inc['LRT_p']:.3g}.")
    L.append(f"- Median/mode-imputation sensitivity (full covariate set, "
             f"n={madj['imputed_sensitivity']['n']}, "
             f"{madj['imputed_sensitivity']['deaths']} deaths): adjusted HR per SD "
             f"**{hrv(madj['imputed_sensitivity']['terms']['risk_sd'])}**.")
    L.append(f"- Reduced core-set sensitivity (age+grade+ER+HER2, complete case, "
             f"n={madj['reduced_core']['n']}, {madj['reduced_core']['deaths']} deaths): "
             f"adjusted HR per SD **{hrv(red)}**; delta C {rinc['delta_C']:+.3f} "
             f"[{rinc['delta_C_CI95_boot500'][0]:+.3f},{rinc['delta_C_CI95_boot500'][1]:+.3f}], "
             f"LRT p={rinc['LRT_p']:.3g}.")
    L.append(f"- **Covariate coding correction disclosed.** {corr['issue']}. "
             f"The superseded complete-case estimate was n={sup_old['n']}, "
             f"{sup_old['deaths']} deaths, score HR {sup_old['risk_sd_HR']:.3f} "
             f"[{sup_old['risk_sd_CI'][0]:.3f},{sup_old['risk_sd_CI'][1]:.3f}] "
             f"(p={sup_old['risk_sd_p']:.3g}) with tumor_size HR "
             f"{sup_old['tumor_size_HR']:.3f} "
             f"[{sup_old['tumor_size_CI'][0]:.3f},{sup_old['tumor_size_CI'][1]:.3f}] "
             f"(p={sup_old['tumor_size_p']:.1e}). Corrected coding: "
             f"{madj['coding']['tumor_size']}. The covariate set itself is unchanged.")
    L.append(f"- Sample-to-patient mapping verified against the live cBioPortal API: "
             f"{smp['distinct_sampleIds']} distinct sampleIds over "
             f"{smp['distinct_patientIds']} patientIds for the widest-coverage attribute "
             f"({smp['patients_with_multiple_samples']} patients with more than one "
             f"sample; max {smp['max_samples_per_patient']}), so the patientId-keyed "
             f"SAMPLE join is unambiguous.")
    L.append("")
    L.append("### S6.2 Formal time-varying (non-proportional-hazards) models")
    L.append("`survival::coxph(Surv(time,event) ~ risk + tt(risk))` with "
             "`tt = risk*log t`, Breslow ties, R 4.6.1 / `survival` "
             f"{tvc['survival_version']}. The likelihood-ratio test compares this model "
             "with the proportional-hazards model (1 df), not with the null model.")
    L.append("")
    L.append("| Cohort | n | deaths | beta(score x log t) | LRT chi-square (df=1) | LRT p | HR(1 y) [95% CI] | HR(10 y) [95% CI] |")
    L.append("|---|---|---|---|---|---|---|---|")
    for key, label in (("GSE20685", "GSE20685"), ("SCANB", "SCAN-B"), ("METABRIC", "METABRIC")):
        d = tvc[key]
        L.append(f"| {label} | {d['n']} | {d['deaths']} | {d['beta_tt']:.3f} | "
                 f"{d['LRT_chisq']:.2f} | {d['LRT_p']:.3g} | "
                 f"{d['hr_at']['t1']:.2f} [{d['ci_lo']['t1']:.2f},{d['ci_hi']['t1']:.2f}] | "
                 f"{d['hr_at']['t10']:.2f} [{d['ci_lo']['t10']:.2f},{d['ci_hi']['t10']:.2f}] |")
    L.append("")
    L.append("*(One patient with zero follow-up was excluded from the log-time "
             "interaction model, yielding n=1,979).*")
    L.append("")
    L.append("Exact Grambsch-Therneau ranked-time tests for the same score: "
             + ", ".join(f"{k.replace('_cont', '')} p={ex[k]['vars']['risk']['p']:.4g}"
                         for k in ("GSE20685_cont", "SCANB_cont", "METABRIC_cont"))
             + ". " + (
                 "Both analyses support stronger early associations in GSE20685 and "
                 "METABRIC; no statistically detectable time variation was observed in "
                 "SCAN-B. Each overall hazard ratio is therefore an average over "
                 "follow-up."))
    mazp = json.load(open(RAW / "metabric_adjusted_zph.json"))
    L.append(f"For the METABRIC complete-case adjusted model above, the same exact "
             f"test gives score p={mazp['vars']['risk_sd']['p']:.4g} and global "
             f"p={mazp['global']['p']:.3g}, confirming that the departure is a property "
             f"of the follow-up distribution rather than an artefact of the unadjusted "
             f"fit.")
    L.append("")
    L.append("### S6.3 Common hazard-ratio scale and locked-cutoff transport")
    L.append("Every adjusted hazard ratio is re-expressed per one SD of the locked score "
             "by refitting the documented model and asserting that the published per-unit "
             "point estimate reproduces. Per-1-score-unit values are in Table S1. The "
             "n/deaths columns describe the full analysis cohort, which is what the "
             "unadjusted HR, raw-score Cox coefficient and cutoff split are computed on; "
             "the adjusted column alone uses the METABRIC complete-case analysis "
             f"(n={aps['METABRIC']['n']}, {aps['METABRIC']['deaths']} deaths).")
    L.append("")
    L.append("| Cohort | n | deaths | HR per SD (unadjusted) | HR per SD (adjusted) | "
             "Raw-score Cox coefficient | Above/below locked cutoff | Cohort SD |")
    L.append("|---|---|---|---|---|---|---|---|")
    unadj = {
        "GSE20685": corr_s["valid_continuous_perSD_HR"],
        "SCANB": rna_s["hr_cont"],
        "METABRIC": [met_s["cont_HR"], met_s["cont_CI"][0], met_s["cont_CI"][1]],
    }
    for key, label, skey in (("GSE20685", "GSE20685", "GSE20685"),
                             ("SCANB", "SCAN-B", "SCANB"),
                             ("METABRIC", "METABRIC", "METABRIC")):
        a = aps[key]
        u = unadj[key]
        val_s = json.load(open(RAW / "validation_summary.json"))
        nrow = {"GSE20685": (val_s["n"], val_s["events"]),
                "SCANB": (rna_s["n"], rna_s["events"]),
                "METABRIC": (met_s["n"], met_s["deaths"])}[key]
        L.append(f"| {label} | {nrow[0]} | {nrow[1]} | "
                 f"{u[0]:.2f} [{u[1]:.2f},{u[2]:.2f}] | "
                 f"{a['per_SD_HR']:.2f} [{a['per_SD_CI'][0]:.2f},{a['per_SD_CI'][1]:.2f}] | "
                 f"{a['raw_score_coefficient']:.2f} | {spl[skey]['text']} "
                 f"({spl[skey]['high_pct']}% above) | {spl[skey]['sd']:.4f} |")
    L.append("")
    L.append(f"The locked cutoff is {spl['locked_cutoff']:.6f}; it was derived on the "
             "discovery cohort and is applied unchanged. The split is severely "
             "imbalanced in opposite directions across platforms, which is the measured "
             "content of the transportability claim, and is why the absolute threshold "
             "is reported only as a secondary analysis.")
    return "\n".join(L)


def build_s7():
    """Figure plates for Figs. S1-S4, captioned with the manuscript's own text."""
    caps = manuscript_legends()
    L = []
    for lab in ("Fig. S1", "Fig. S2", "Fig. S3", "Fig. S4"):
        if lab not in caps:
            raise SystemExit(f"no manuscript legend for {lab}; the supplement plate "
                             "would render without a caption")
        stem = SUPP_FIG_FILES[lab]
        raster = ROOT / "figures" / f"{stem}.png"
        if not raster.exists():
            raise SystemExit(f"missing raster for {lab}: figures/{stem}.png")
        n = lab.replace("Fig. S", "")
        L += [f"### Figure S{n}", "", caps[lab], "",
              f"![Figure S{n}](../figures/{stem}.png)", ""]
    return "\n".join(L)


# ---- everything verified: write it all now ----
commit_jobs()

# Section order is a packaging requirement, not content, so it is asserted here
# rather than trusted to the file's current shape (BCRT: Statements and
# Declarations after the References).  Runs after commit_jobs() so it sees the
# final text, and is a no-op once the order is right.
move_section_after(sub, "Statements and Declarations", "References")

# The supplement regions are filled only now: S7 restates the manuscript figure
# captions, so it has to read the manuscript *after* every staged edit has landed.
_src = SUP_SRC.read_text()
for _name, _builder in (("S3", build_s3), ("S6", build_s6), ("S7", build_s7)):
    _src = fill_region(_src, _name, _builder, SUP_SRC)
for _name in REGIONS:
    _m = re.search(rf"<!-- {_name}-GENERATED-START -->(.*?)<!-- {_name}-GENERATED-END -->",
                   _src, re.S)
    if not _m or len(_m.group(1).strip()) < 100:
        raise SystemExit(f"supplement region {_name} is empty after filling")
SUP_SRC.write_text(_src)
shutil.copyfile(SUP_SRC, SUP)
log(f"[supplement] {len(REGIONS)} generated regions rebuilt; source "
    f"{len(_src)} chars, package copy {SUP.stat().st_size} bytes")

remark_src = MS / "REMARK_checklist.md"
remark = SUB / "REMARK_checklist.md"
shutil.copyfile(remark_src, remark)
log(f"[REMARK] checklist copied into the submission package ({remark.stat().st_size} bytes)")

build_legends()

# ------------------------------------------------------------------ package index
# The index used to be submission_bcrt/README.md, i.e. a file describing the build
# sitting inside the directory that gets uploaded, and it named internal scripts and
# a user-action ledger.  submission_bcrt/ now holds only uploadable files; the index
# lives with the other working notes.
index = ROOT / "review" / "PACKAGE_INDEX.md"
index.write_text(
    "# BCRT submission package index (internal working note; do not upload)\n\n"
    "`submission_bcrt/` contains only files that can be uploaded. Everything in this "
    "directory is working material.\n\n"
    "## Uploadable files\n\n"
    "- `cover_letter.md` - cover letter, incl. three suggested reviewers.\n"
    "- `title_page.md` - title, authors, affiliations, ORCIDs, keywords.\n"
    "- `manuscript_submission.md` - the single source of truth for every rendered file.\n"
    "- `manuscript_final.pdf` - publication-grade PDF.\n"
    "- `manuscript_review.pdf` - same content, review rendering with figures inline.\n"
    "- `manuscript_bcrt.doc` - Word manuscript, tables on landscape pages.\n"
    "- `manuscript_bcrt.tex` - LaTeX source.\n"
    "- `abstract_structured.md` - structured abstract (generated from the manuscript).\n"
    "- `references_numbered.md` - numbered reference list, order of first appearance.\n"
    "- `legends.md` - figure legends and table titles, beside each file name.\n"
    "- `declarations.md` - Statements and Declarations.\n"
    "- `REMARK_checklist.md` - completed REMARK checklist (all 20 items).\n"
    "- `supplement.md` - Supplementary Material source.\n"
    "- `ESM_1.pdf` - Supplementary Information (Online Resource 1), as uploaded.\n\n"
    "## Regeneration order\n\n"
    "`02`-`42` produce the results and tables; then `34_build_abstract.py`, "
    "`36_citations.py`, `44_build_remark.py`, `39_finalize_docs.py`, "
    "`23_build_latex.py`, `24_build_docx.py`, `25_build_final_pdf.py`, "
    "`22_build_review_pdf.py`, `46_build_supplement_pdf.py`, "
    "`40_submission_gate.py`.\n\n"
    "## Working notes\n\n"
    "- `FIGURE_PROVENANCE.md` - figure/table file to generator and source.\n"
    "- `CITATION_MAP.md` - reference number to the sentence that cites it.\n\n"
    "Status: submission-ready for Breast Cancer Research and Treatment. Remaining "
    "author-only items are tracked in the project's user-action ledger.\n")
log(f"[index] regenerated review/PACKAGE_INDEX.md ({index.stat().st_size} bytes)")

# the superseded in-package README must not be left behind
stale = SUB / "README.md"
if stale.exists():
    stale.unlink()
    log("[index] removed the superseded submission_bcrt/README.md")

logf.close()
