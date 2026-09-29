"""Build a single review PDF from the submission manuscript (REVIEWPDF-001).

Renders submission_bcrt/manuscript_submission.md with DejaVuSans (covers →,
×, –, χ², Δ, β), embeds the 12 main-figure PNGs at their legends, and renders the
4 submission tables from tables/*.tsv. Review rendering only; no scientific
content is created or altered here.
Outputs: submission_bcrt/manuscript_review.pdf, logs/review_pdf.log
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
                                Table, TableStyle, PageBreak, KeepTogether)
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# Reproducibility: reportlab stamps /CreationDate with the wall-clock time, so two
# identical builds of the same Markdown produced different bytes. reportlab honours
# SOURCE_DATE_EPOCH for the embedded date, so pin it before any canvas is created.
# The document renders no date, so nothing visible changes.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1700000000")

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "submission_bcrt" / "manuscript_submission.md"
FIG = ROOT / "figures"
TAB = ROOT / "tables"
OUT = ROOT / "submission_bcrt" / "manuscript_review.pdf"
LOGS = ROOT / "logs"
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

logf = open(LOGS / "review_pdf.log", "w")


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
    t = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1 (\2)", t)

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


body = ParagraphStyle("body", fontName="DVS", fontSize=10, leading=14,
                      alignment=TA_JUSTIFY, spaceAfter=6)
h1 = ParagraphStyle("h1", fontName="DVS", fontSize=16, leading=20,
                    spaceBefore=14, spaceAfter=8)
h2 = ParagraphStyle("h2", fontName="DVS", fontSize=13, leading=17,
                    spaceBefore=12, spaceAfter=6)
h3 = ParagraphStyle("h3", fontName="DVS", fontSize=11, leading=14,
                    spaceBefore=10, spaceAfter=4)
cap = ParagraphStyle("cap", fontName="DVS", fontSize=8.5, leading=11.5,
                     spaceAfter=8, textColor=colors.grey)

tbl_cell = ParagraphStyle("tbl_cell", fontName="DVS", fontSize=6.5, leading=8.5,
                          alignment=TA_CENTER)
tbl_cell_left = ParagraphStyle("tbl_cell_left", fontName="DVS", fontSize=6.5, leading=8.5,
                               alignment=TA_LEFT)
tbl_hdr = ParagraphStyle("tbl_hdr", fontName="DVS", fontSize=6.5, leading=8.5,
                         alignment=TA_CENTER)
tbl_note = ParagraphStyle("tbl_note", fontName="DVS", fontSize=7, leading=9.5,
                          textColor=colors.dimgrey, spaceAfter=8)


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
    # string "nan" that pandas emits for coerced NaN.  (The final-PDF builder was
    # fixed first; this one kept the defect until it was caught by the gate, which
    # is why the gate now scans every deliverable.)
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
                if row["Note"] and row["Note"] != "nan":
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
i, nfig, ntbl = 0, 0, 0
table_heading = None
while i < len(lines):
    ln = lines[i].rstrip()
    if ln.startswith("# "):
        story.append(Paragraph(rich(ln[2:]), h1))
    elif ln.startswith("### "):
        story.append(Paragraph(rich(ln[4:]), h3))
    elif ln.startswith("## "):
        if ln[3:].strip() == "Tables":
            table_heading = Paragraph("Tables", h2)
        else:
            story.append(Paragraph(rich(ln[3:]), h2))
    elif ln.startswith("- **Fig"):
        m = re.match(r"- \*\*(Fig\.? ?S?\d+)\*\* (.*)", ln)
        if m:
            fig_id = m.group(1).replace("Fig ", "Fig. ")
            fig_cap = rich(m.group(2))
            num = re.sub(r"^Fig\.?\s*", "", fig_id)
            caption = Paragraph(f"<b>{fig_id}</b> {fig_cap}", cap)
            
            fig_file = FIG_MAP.get(fig_id)
            if not fig_file:
                f_search = re.search(r"Files: ([A-Za-z0-9_]+\.png)", ln)
                if f_search:
                    fig_file = f_search.group(1)
            
            if not fig_file or not (FIG / fig_file).exists():
                raise SystemExit(f"no image for {fig_id} (FIG_MAP -> {fig_file!r}); a "
                                 "legend must not render as a caption with no figure")
            figure = Image(str(FIG / fig_file), width=150 * mm,
                           height=100 * mm, kind="proportional")
            story.append(KeepTogether([caption, figure, Spacer(1, 6)]))
            nfig += 1
    elif ln.startswith("- **Table"):
        m = re.match(r"- \*\*(Table \d+)\.\*\* (.*)", ln)
        if m:
            tbl_id = m.group(1)
            tbl_cap = rich(m.group(2))
            caption = Paragraph(f"<b>{tbl_id}.</b> {tbl_cap}", cap)
            
            tbl_file = TBL_MAP.get(tbl_id)
            if not tbl_file:
                f_search = re.search(r"File: tables/([A-Za-z0-9_]+\.tsv)", ln)
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
    elif ln.startswith("- "):
        story.append(Paragraph("• " + rich(ln[2:]), body))
    elif ln.strip() == "":
        story.append(Spacer(1, 4))
    else:
        story.append(Paragraph(rich(ln), body))
    i += 1

doc = SimpleDocTemplate(str(OUT), pagesize=A4, topMargin=20 * mm,
                        bottomMargin=20 * mm, leftMargin=20 * mm,
                        rightMargin=20 * mm)
n_legends = sum(1 for l in lines if l.startswith("- **Fig"))
if nfig != n_legends:
    raise SystemExit(f"{nfig} image(s) embedded for {n_legends} figure legend(s)")
doc.build(story)
log(f"REVIEWPDF-001: {OUT.name} with {nfig} figures + {ntbl} tables; "
    f"{OUT.stat().st_size / 1024:.0f} KB")
logf.close()
