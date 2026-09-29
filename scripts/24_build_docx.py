"""Build the BCRT Word manuscript (WORDSRC-001).

python-docx: headings, paragraphs with **bold** runs, bullets, 4 tables from
 tables/*.tsv, 12 main-figure PNGs at their legends. Review/submission rendering
only; no scientific content is created or altered here. A temporary DOCX is
converted to the journal's requested legacy .doc format with LibreOffice; the
DOC is round-tripped to DOCX to verify that tables and figures survived.
Outputs: submission_bcrt/manuscript_bcrt.doc, logs/word_src.log
"""
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.section import WD_ORIENT, WD_SECTION
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "submission_bcrt" / "manuscript_submission.md"
FIG = ROOT / "figures"
TAB = ROOT / "tables"
OUT = ROOT / "submission_bcrt" / "manuscript_bcrt.doc"
LOGS = ROOT / "logs"
logf = open(LOGS / "word_src.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


def rich_para(p, t):
    for tok in re.split(r"(\*\*.+?\*\*)", t):
        if tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            r = p.add_run(tok[2:-2])
            r.bold = True
        else:
            p.add_run(tok)


FIG_MAP = {
    "Fig. 1": "Fig1_flow.png",
    "Fig. 2": "Fig2_volcano_42568.png",
    "Fig. 3": "Fig3_coef_forest_train.png",
    "Fig. 4": "train_KM.png",
    "Fig. 5": "Fig3_ROC_train.png",
    "Fig. 6": "valid_KM.png",
    "Fig. 7": "Fig4_ROC_valid.png",
    "Fig. 8": "Fig5_pathway_valid.png",
    "Fig. 9": "Fig5_checkpoints_valid.png",
    "Fig. 10": "RNASEQ_KM.png",
    "Fig. 11": "Fig_PAM50_forest.png",
    "Fig. 12": "Fig12_crosscohort_forest.png",
    "Fig. S1": "QC_followup_expr.png",
    "Fig. S2": "risk_distributions.png",
    "Fig. S3": "PH_schoenfeld.png",
    "Fig. S4": "HR_t_curves.png",
}

TBL_MAP = {
    "Table 1": "Tab1_cohorts_v3.tsv",
    "Table 2": "Tab2_coefficients.tsv",
    "Table 3": "Tab3_performance_v3.tsv",
    "Table 4": "Tab4_incremental_value.tsv",
}

doc = Document()
st = doc.styles["Normal"]
st.font.name = "Times New Roman"
st._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
st.font.size = Pt(10)
nfig = ntbl = 0


def set_orientation(section, landscape):
    """Word submission tables are unreadable in portrait (external review P0-3),
    so the Tables section is laid out landscape and the rest stays portrait."""
    w, h = section.page_width, section.page_height
    if landscape:
        section.orientation = WD_ORIENT.LANDSCAPE
        section.page_width, section.page_height = max(w, h), min(w, h)
    else:
        section.orientation = WD_ORIENT.PORTRAIT
        section.page_width, section.page_height = min(w, h), max(w, h)


def style_table(tb, ncols):
    tb.autofit = True
    for row in tb.rows:
        for c in row.cells:
            for p in c.paragraphs:
                p.paragraph_format.space_before = Pt(1)
                p.paragraph_format.space_after = Pt(1)
                for r in p.runs:
                    r.font.size = Pt(8)


lines = [ln for ln in MS.read_text().split("\n")
         if not re.match(r"^- \*\*Fig\. S\d+\*\*", ln)]
for ln in lines:
    ln = ln.rstrip()
    if ln.startswith("# "):
        doc.add_heading(ln[2:], level=0)
    elif ln.startswith("### "):
        doc.add_heading(ln[4:], level=2)
    elif ln.startswith("## "):
        sec_name = ln[3:].strip()
        if sec_name == "Tables":
            set_orientation(doc.add_section(WD_SECTION.NEW_PAGE), True)
        elif sec_name == "References":
            set_orientation(doc.add_section(WD_SECTION.NEW_PAGE), False)
        doc.add_heading(sec_name, level=1)
    elif ln.startswith("- **Fig"):
        m = re.match(r"- \*\*(Fig\.? ?S?\d+)\*\* (.*)", ln)
        if m:
            fig_id = m.group(1).replace("Fig ", "Fig. ")
            p = doc.add_paragraph()
            p.paragraph_format.keep_with_next = True
            r = p.add_run(fig_id + " ")
            r.bold = True
            rich_para(p, m.group(2))
            
            fig_file = FIG_MAP.get(fig_id)
            if not fig_file:
                f = re.search(r"Files: ([A-Za-z0-9_]+\.png)", ln)
                if f:
                    fig_file = f.group(1)
            
            if not fig_file or not (FIG / fig_file).exists():
                raise SystemExit(f"no image for {fig_id} (FIG_MAP -> {fig_file!r}); a "
                                 "legend must not render as a caption with no figure")
            doc.add_picture(str(FIG / fig_file), width=Inches(6.0))
            doc.paragraphs[-1].paragraph_format.keep_together = True
            nfig += 1
    elif ln.startswith("- **Table"):
        m = re.match(r"- \*\*(Table \d+)\.\*\* (.*)", ln)
        if m:
            tbl_id = m.group(1)
            p = doc.add_paragraph()
            p.paragraph_format.keep_with_next = True
            r = p.add_run(tbl_id + ". ")
            r.bold = True
            rich_para(p, m.group(2))
            
            tbl_file = TBL_MAP.get(tbl_id)
            if not tbl_file:
                f = re.search(r"File: tables/([A-Za-z0-9_]+\.tsv)", ln)
                if f:
                    tbl_file = f.group(1)
            
            if tbl_file and (TAB / tbl_file).exists():
                d = pd.read_csv(TAB / tbl_file, sep="\t", keep_default_na=False,
                                dtype=str).replace({"": "NA"})
                notes = []
                if tbl_id == "Table 3" and "Note" in d.columns:
                    for _, row in d.iterrows():
                        if row["Note"]:
                            note_str = row["Note"].replace("LumA/B null", "LumA/B no sig association")
                            notes.append(f"{row['Cohort']}: {note_str}")
                    d_disp = d.drop(columns=["Note"])
                else:
                    d_disp = d
                
                tb = doc.add_table(rows=1 + len(d_disp), cols=len(d_disp.columns))
                tb.style = "Table Grid"
                for j, c in enumerate(d_disp.columns):
                    tb.rows[0].cells[j].text = c
                for i, (_, r_data) in enumerate(d_disp.iterrows()):
                    for j, x in enumerate(r_data.tolist()):
                        tb.rows[i + 1].cells[j].text = x
                style_table(tb, len(d_disp.columns))
                
                if notes:
                    np = doc.add_paragraph()
                    nr = np.add_run("Table 3 Notes: ")
                    nr.bold = True
                    np.add_run("; ".join(notes))
                
                ntbl += 1
    elif ln.startswith("- "):
        doc.add_paragraph(ln[2:], style="List Bullet")
    elif ln.strip() == "":
        continue
    else:
        p = doc.add_paragraph()
        rich_para(p, ln)

def normalize_zip(path, stamp=(1980, 1, 1, 0, 0, 0)):
    """Rewrite the .docx so its bytes do not depend on the build time.

    python-docx stamps every zip entry with the current local time, so two
    identical builds produced different sha256 even though all 33 member payloads
    were byte-identical (verified by hashing each member of two consecutive
    builds). Rewriting the archive with a fixed ZipInfo timestamp makes the .docx
    byte-reproducible, like the PDFs, so its hash can be pinned in the ledger.
    """
    tmp = path.with_suffix(path.suffix + ".tmp")
    with zipfile.ZipFile(path) as src, \
            zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as dst:
        for info in src.infolist():
            new = zipfile.ZipInfo(info.filename, date_time=stamp)
            new.compress_type = info.compress_type
            new.external_attr = info.external_attr
            new.internal_attr = info.internal_attr
            new.create_system = info.create_system
            dst.writestr(new, src.read(info.filename))
    tmp.replace(path)


n_legends = sum(1 for l in lines if l.startswith("- **Fig"))
if nfig != n_legends:
    raise SystemExit(f"{nfig} image(s) embedded for {n_legends} figure legend(s)")

# BCRT requests automatic page numbering. Use a PAGE field in each section so
# numbering remains live when editors reflow the Word document.
for section in doc.sections:
    section.footer.is_linked_to_previous = False
    footer = section.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    page = OxmlElement("w:fldSimple")
    page.set(qn("w:instr"), "PAGE")
    run = OxmlElement("w:r")
    txt = OxmlElement("w:t")
    txt.text = "1"
    run.append(txt)
    page.append(run)
    footer._p.append(page)

title_match = re.search(r"(?m)^## Title\s*\n([^\n]+)", MS.read_text())
if not title_match:
    raise SystemExit("submission manuscript has no ## Title section")
expected_title = title_match.group(1).strip()

with tempfile.TemporaryDirectory(prefix="bio-dsh-word-") as tmp:
    tmp = Path(tmp)
    src = tmp / "manuscript_bcrt.docx"
    outdir = tmp / "doc"
    checkdir = tmp / "roundtrip"
    outdir.mkdir()
    checkdir.mkdir()
    doc.save(str(src))
    normalize_zip(src)

    soffice = shutil.which("soffice")
    if not soffice:
        raise SystemExit("LibreOffice 'soffice' is required to create the BCRT .doc file")
    conv = subprocess.run(
        [soffice, "--headless", "--convert-to", "doc", "--outdir", str(outdir), str(src)],
        capture_output=True, text=True)
    built_doc = outdir / "manuscript_bcrt.doc"
    if conv.returncode or not built_doc.is_file():
        raise SystemExit("DOC conversion failed: " + (conv.stderr or conv.stdout)[-500:])

    # A round-trip catches converters that produce a syntactically valid DOC but
    # silently drop embedded figures or tables.
    check = subprocess.run(
        [soffice, "--headless", "--convert-to", "docx", "--outdir", str(checkdir),
         str(built_doc)], capture_output=True, text=True)
    roundtrip = checkdir / "manuscript_bcrt.docx"
    if check.returncode or not roundtrip.is_file():
        raise SystemExit("DOC round-trip validation failed: " +
                         (check.stderr or check.stdout)[-500:])
    verify = Document(str(roundtrip))
    body = "\n".join(p.text for p in verify.paragraphs)
    if expected_title not in body:
        raise SystemExit("DOC round-trip lost the article title")
    if len(verify.inline_shapes) != nfig or len(verify.tables) != ntbl:
        raise SystemExit("DOC round-trip changed figures/tables: "
                         f"{len(verify.inline_shapes)}/{len(verify.tables)} "
                         f"vs expected {nfig}/{ntbl}")
    if not verify.sections or any("PAGE" not in section.footer._element.xml
                                  for section in verify.sections):
        raise SystemExit("DOC round-trip lost an automatic page-number field")
    shutil.copy2(built_doc, OUT)

log(f"WORDSRC-001: {OUT.name} with {nfig} figures + {ntbl} tables; "
    f"{OUT.stat().st_size / 1024:.0f} KB; DOC round-trip verified")
logf.close()
