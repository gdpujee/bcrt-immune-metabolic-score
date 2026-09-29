"""Build the formal submission-style article PDF (FINALPDF-001).

Publication-grade layout from submission_bcrt/manuscript_submission.md:
title block with authors/affiliations/corresponding line, keywords, abstract,
numbered IMRaD sections, Declarations, captioned figures and tables, then
numbered References at the end, with running header + page numbers (DejaVuSans).
Review/submission rendering only; no scientific content created or altered.
Outputs: submission_bcrt/manuscript_final.pdf, logs/final_pdf.log
"""
import re
import os
import importlib.util
import pandas as pd
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak, HRFlowable,
                                KeepTogether)
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# Reproducibility: reportlab stamps /CreationDate with the wall-clock time, so two
# identical builds of the same Markdown produced different bytes (verified: two
# consecutive runs gave 167228ac… then 05e55a98…). reportlab honours
# SOURCE_DATE_EPOCH for the embedded date, so pin it before any canvas is created.
# The document renders no date, so nothing visible changes.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1700000000")

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "submission_bcrt" / "manuscript_submission.md"
FIG = ROOT / "figures"
TAB = ROOT / "tables"
OUT = ROOT / "submission_bcrt" / "manuscript_final.pdf"
LOGS = ROOT / "logs"
_title_match = re.search(r"(?m)^## Title\s*\n([^\n]+)", MS.read_text())
if not _title_match:
    raise SystemExit("submission manuscript has no ## Title section")
article_title = _title_match.group(1).strip()
FONT = None
if importlib.util.find_spec("matplotlib"):
    import matplotlib
    FONT = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data",
                        "fonts", "ttf", "DejaVuSans.ttf")
if FONT is None or not os.path.isfile(FONT):
    _reportlab_fonts = Path(pdfmetrics.__file__).parents[1] / "fonts"
    FONT = next((p for p in ("/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
                             "/Library/Fonts/Arial Unicode.ttf",
                             str(_reportlab_fonts / "Vera.ttf"))
                 if os.path.isfile(p)), None)
if FONT is None:
    raise SystemExit("no Unicode TrueType font found (install matplotlib or provide a system font)")
pdfmetrics.registerFont(TTFont("DVS", FONT))
pdfmetrics.registerFontFamily("DVS", normal="DVS", bold="DVS",
                              italic="DVS", boldItalic="DVS")
logf = open(LOGS / "final_pdf.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5",
       "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-"}


def esc(t):
    return (t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def rich(t):
    t = esc(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)

    # Scientific notation and superscripts
    t = re.sub(r"([0-9.]+)&times;10([⁻¹²³⁴⁵⁶⁷⁸⁹⁰]+)",
               lambda m: r"%s&times;10<sup>%s</sup>" % (m.group(1), "".join(SUP.get(c, c) for c in m.group(2))), t)
    t = re.sub(r"([0-9.]+)×10([⁻¹²³⁴⁵⁶⁷⁸⁹⁰]+)",
               lambda m: r"%s&times;10<sup>%s</sup>" % (m.group(1), "".join(SUP.get(c, c) for c in m.group(2))), t)
    t = re.sub(r"×10([⁻¹²³⁴⁵⁶⁷⁸⁹⁰]+)",
               lambda m: r"&times;10<sup>%s</sup>" % "".join(SUP.get(c, c) for c in m.group(1)), t)
    t = re.sub(r"\b([0-9.]+)e-([0-9]+)\b", r"\1&times;10<sup>-\2</sup>", t)

    # Greek & math entities
    t = t.replace("χ²", "&chi;<sup>2</sup>")
    t = t.replace("χ", "&chi;")
    t = t.replace("ΔC", "&Delta;C")
    t = t.replace("Δ", "&Delta;")
    t = t.replace("×", "&times;")
    t = t.replace("±", "&plusmn;")
    t = t.replace("≥", "&ge;")
    t = t.replace("≤", "&le;")
    t = t.replace("→", "&rarr;")
    t = t.replace("–", "&ndash;").replace("—", "&mdash;")
    return t


title_st = ParagraphStyle("title", fontName="DVS", fontSize=16, leading=20,
                          alignment=TA_CENTER, spaceAfter=8)
auth_st = ParagraphStyle("auth", fontName="DVS", fontSize=11, leading=14,
                         alignment=TA_CENTER, spaceAfter=2)
aff_st = ParagraphStyle("aff", fontName="DVS", fontSize=9, leading=12,
                        alignment=TA_CENTER, spaceAfter=6)
abs_st = ParagraphStyle("abs", fontName="DVS", fontSize=9.5, leading=13.5,
                        alignment=TA_JUSTIFY, spaceAfter=6)
body = ParagraphStyle("body", fontName="DVS", fontSize=10, leading=14.5,
                      alignment=TA_JUSTIFY, spaceAfter=6)
h1 = ParagraphStyle("h1", fontName="DVS", fontSize=13, leading=17,
                    spaceBefore=12, spaceAfter=6)
h3 = ParagraphStyle("h3", fontName="DVS", fontSize=11, leading=14,
                    spaceBefore=10, spaceAfter=4)
cap = ParagraphStyle("cap", fontName="DVS", fontSize=8.5, leading=11.5,
                     spaceAfter=8)
ref_st = ParagraphStyle("ref", fontName="DVS", fontSize=8.5, leading=11.5,
                        leftIndent=14, firstLineIndent=-14, spaceAfter=4)

tbl_cell = ParagraphStyle("tbl_cell", fontName="DVS", fontSize=6.5, leading=8.5,
                          alignment=TA_CENTER)
tbl_cell_left = ParagraphStyle("tbl_cell_left", fontName="DVS", fontSize=6.5, leading=8.5,
                               alignment=TA_LEFT)
tbl_hdr = ParagraphStyle("tbl_hdr", fontName="DVS", fontSize=6.5, leading=8.5,
                         alignment=TA_CENTER)
tbl_note = ParagraphStyle("tbl_note", fontName="DVS", fontSize=7, leading=9.5,
                          textColor=colors.dimgrey, spaceAfter=8)


def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("DVS", 8)
    canvas.setFillColor(colors.grey)
    canvas.drawString(20 * mm, 282 * mm,
                      "Manuscript prepared for submission to Breast Cancer Research and Treatment")
    canvas.drawRightString(190 * mm, 15 * mm, "Page %d" % doc.page)
    canvas.restoreState()


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


def make_table_flowables(tbl_id, tsv_file):
    # keep_default_na=False: literal "NA" cells (e.g. GSE45827 has no survival)
    # must survive as "NA", and empty cells must render as "NA" rather than the
    # string "nan" that pandas emits for coerced NaN.
    d = pd.read_csv(TAB / tsv_file, sep="\t", keep_default_na=False,
                    dtype=str).replace({"": "NA"})
    flowables = []

    if tbl_id == "Table 1":
        col_widths = [30 * mm, 14 * mm, 15 * mm, 15 * mm, 20 * mm, 22 * mm, 13 * mm, 41 * mm]
        hdr = [Paragraph(f"<b>{rich(c)}</b>", tbl_hdr) for c in d.columns]
        rows = []
        for r in d.values:
            rows.append([Paragraph(rich(r[j]), tbl_cell_left if j in (0, 5, 7) else tbl_cell)
                         for j in range(len(r))])
        t = Table([hdr] + rows, colWidths=col_widths, repeatRows=1)
        t.setStyle(TableStyle([
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            # reportlab's default 6pt side padding left only ~6 mm of text width in
            # the 9-13 mm numeric columns, so headers and numbers wrapped mid-token
            # ("deat hs", "Slo pe", "1909/ 71").  2pt keeps every token intact.
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ]))
        flowables.append(t)
        flowables.append(Spacer(1, 6))

    elif tbl_id == "Table 2":
        col_widths = [26 * mm, 24 * mm, 24 * mm, 24 * mm, 24 * mm, 24 * mm, 24 * mm]
        hdr = [Paragraph(f"<b>{rich(c)}</b>", tbl_hdr) for c in d.columns]
        rows = []
        for r in d.values:
            rows.append([Paragraph(rich(r[j]), tbl_cell_left if j == 0 else tbl_cell)
                         for j in range(len(r))])
        t = Table([hdr] + rows, colWidths=col_widths, repeatRows=1)
        t.setStyle(TableStyle([
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            # reportlab's default 6pt side padding left only ~6 mm of text width in
            # the 9-13 mm numeric columns, so headers and numbers wrapped mid-token
            # ("deat hs", "Slo pe", "1909/ 71").  2pt keeps every token intact.
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ]))
        flowables.append(t)
        flowables.append(Spacer(1, 6))

    elif tbl_id == "Table 3":
        notes = []
        if "Note" in d.columns:
            for _, row in d.iterrows():
                if row["Note"]:
                    note_str = row["Note"]
                    note_str = note_str.replace("LumA/B null", "LumA/B no sig association")
                    notes.append(f"<b>{rich(row['Cohort'])}:</b> {rich(note_str)}")
            d_disp = d.drop(columns=["Note"])
        else:
            d_disp = d

        # Give the raw-score coefficient a usable label width. The old 9 mm
        # column split "Raw-score Cox coef." inside words in the rendered PDF.
        col_widths = [24 * mm, 9 * mm, 11 * mm, 20 * mm, 22 * mm, 20 * mm, 18 * mm, 14 * mm, 16 * mm, 15 * mm]
        hdr = [Paragraph(f"<b>{rich(c)}</b>", tbl_hdr) for c in d_disp.columns]
        rows = []
        for r in d_disp.values:
            rows.append([Paragraph(rich(r[j]), tbl_cell_left if j == 0 else tbl_cell)
                         for j in range(len(r))])
        t = Table([hdr] + rows, colWidths=col_widths, repeatRows=1)
        t.setStyle(TableStyle([
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            # reportlab's default 6pt side padding left only ~6 mm of text width in
            # the 9-13 mm numeric columns, so headers and numbers wrapped mid-token
            # ("deat hs", "Slo pe", "1909/ 71").  2pt keeps every token intact.
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ]))
        flowables.append(t)
        if notes:
            flowables.append(Spacer(1, 3))
            flowables.append(Paragraph("<b>Table 3 Notes:</b> " + " ".join(notes), tbl_note))
        flowables.append(Spacer(1, 4))

    elif tbl_id == "Table 4":
        col_widths = [65 * mm, 18 * mm, 32 * mm, 55 * mm]
        hdr = [Paragraph(f"<b>{rich(c)}</b>", tbl_hdr) for c in d.columns]
        rows = []
        for r in d.values:
            rows.append([Paragraph(rich(r[j]), tbl_cell_left if j in (0, 3) else tbl_cell)
                         for j in range(len(r))])
        t = Table([hdr] + rows, colWidths=col_widths, repeatRows=1)
        t.setStyle(TableStyle([
            ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
            ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            # reportlab's default 6pt side padding left only ~6 mm of text width in
            # the 9-13 mm numeric columns, so headers and numbers wrapped mid-token
            # ("deat hs", "Slo pe", "1909/ 71").  2pt keeps every token intact.
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ]))
        flowables.append(t)
        flowables.append(Spacer(1, 6))

    return flowables


story = []
lines = [ln for ln in MS.read_text().split("\n")
         if not re.match(r"^- \*\*Fig\. S\d+\*\*", ln)]
secnum, nfig, ntbl = 0, 0, 0
table_heading = None
NUMERED = {"Introduction", "Methods", "Results", "Discussion", "Conclusion"}
i = 0
guard = 0
while i < len(lines):
    guard += 1
    assert guard < 100000, "loop guard tripped"
    ln = lines[i].rstrip()
    if ln.startswith("# ") and not ln.startswith("## "):
        i += 1
        continue
    elif ln.startswith("## Title"):
        j = i + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        story.append(Paragraph(rich(lines[j].strip()), title_st))
        story.append(HRFlowable(width="100%", thickness=0.6,
                                color=colors.grey, spaceAfter=8))
        i = j
    elif ln.startswith("## Authors"):
        j = i + 1
        buf = []
        while j < len(lines) and not lines[j].startswith("## "):
            if lines[j].strip():
                buf.append(lines[j].strip())
            j += 1
        for k, b in enumerate(buf):
            story.append(Paragraph(rich(b), auth_st if k == 0 else aff_st))
        i = j - 1
    elif ln.startswith("## Keywords"):
        j = i + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        story.append(Paragraph("<b>Keywords:</b> " + rich(lines[j].strip()), body))
        i = j
    elif ln.startswith("## Abstract"):
        j = i + 1
        buf = []
        while j < len(lines) and not lines[j].startswith("## "):
            if lines[j].strip():
                buf.append(lines[j].strip())
            j += 1
        story.append(Paragraph("<b>Abstract</b><br/>" + "<br/><br/>".join(rich(b) for b in buf), abs_st))
        i = j - 1
    elif ln.startswith("## References"):
        story.append(Paragraph("References", h1))
        j = i + 1
        while j < len(lines) and not lines[j].startswith("## "):
            m = re.match(r"\[(\d+)\] (.*)", lines[j].strip())
            if m:
                story.append(Paragraph(f"[{m.group(1)}] {rich(m.group(2))}", ref_st))
            j += 1
        i = j - 1
    elif ln.startswith("## Figure legends"):
        # A blank Markdown line before this explicit page break becomes a tiny
        # Spacer flowable. If it cannot fit at the foot of the preceding page,
        # ReportLab can put it alone on a new page and then honor PageBreak(),
        # leaving a visually empty sheet before the legends.
        while story and isinstance(story[-1], Spacer):
            story.pop()
        story.append(PageBreak())
        story.append(Paragraph("Figure legends", h1))
        j = i + 1
        while j < len(lines) and not lines[j].startswith("## "):
            m = re.match(r"- \*\*(Fig\.? ?S?\d+)\*\* (.*)", lines[j].strip())
            if m:
                fig_id = m.group(1).replace("Fig ", "Fig. ")
                fig_cap = rich(m.group(2))
                num = re.sub(r"^Fig\.?\s*", "", fig_id)
                caption = Paragraph(f"<b>{fig_id}</b> {fig_cap}", cap)
                
                fig_file = FIG_MAP.get(fig_id)
                if not fig_file:
                    f_search = re.search(r"Files: ([A-Za-z0-9_]+\.png)", lines[j])
                    if f_search:
                        fig_file = f_search.group(1)
                
                if not fig_file or not (FIG / fig_file).exists():
                    raise SystemExit(f"no image for {fig_id} (FIG_MAP -> {fig_file!r}); a "
                                     "legend must not render as a caption with no figure")
                figure = Image(str(FIG / fig_file), width=150 * mm,
                               height=100 * mm, kind="proportional")
                story.append(KeepTogether([caption, figure, Spacer(1, 8)]))
                nfig += 1
            j += 1
        i = j - 1
    elif ln.startswith("## Tables"):
        while story and isinstance(story[-1], Spacer):
            story.pop()
        story.append(PageBreak())
        table_heading = Paragraph("Tables", h1)
        j = i + 1
        while j < len(lines) and not lines[j].startswith("## "):
            m = re.match(r"- \*\*(Table \d+)\.\*\* (.*)", lines[j].strip())
            if m:
                tbl_id = m.group(1)
                tbl_cap = rich(m.group(2))
                caption = Paragraph(f"<b>{tbl_id}.</b> {tbl_cap}", cap)
                
                tbl_file = TBL_MAP.get(tbl_id)
                if not tbl_file:
                    f_search = re.search(r"File: tables/([A-Za-z0-9_]+\.tsv)", lines[j])
                    if f_search:
                        tbl_file = f_search.group(1)
                
                if tbl_file and (TAB / tbl_file).exists():
                    t_flowables = make_table_flowables(tbl_id, tbl_file)
                    group = [caption, *t_flowables]
                    if table_heading is not None:
                        group.insert(0, table_heading)
                        table_heading = None
                    story.append(KeepTogether(group))
                    ntbl += 1
            j += 1
        i = j - 1
    elif ln.startswith("### "):
        story.append(Paragraph(rich(ln[4:]), h3))
    elif ln.startswith("## "):
        name = ln[3:].strip()
        if name in NUMERED:
            secnum += 1
            story.append(Paragraph(f"{secnum} {rich(name)}", h1))
        else:
            story.append(Paragraph(rich(name), h1))
    elif ln.startswith("- "):
        story.append(Paragraph("• " + rich(ln[2:]), body))
    elif ln.strip() == "":
        story.append(Spacer(1, 4))
    else:
        story.append(Paragraph(rich(ln), body))
    i += 1

doc = SimpleDocTemplate(str(OUT), pagesize=A4, topMargin=25 * mm,
                        bottomMargin=20 * mm, leftMargin=20 * mm,
                        rightMargin=20 * mm, title=article_title,
                        author="Danhua He, Qiang Li")
n_legends = sum(1 for l in lines if l.startswith("- **Fig"))
if nfig != n_legends:
    raise SystemExit(f"{nfig} image(s) embedded for {n_legends} figure legend(s)")
doc.build(story, onFirstPage=header_footer, onLaterPages=header_footer)
log(f"FINALPDF-001: {OUT.name} with {nfig} figures + {ntbl} tables; "
    f"{OUT.stat().st_size / 1024:.0f} KB")
logf.close()
