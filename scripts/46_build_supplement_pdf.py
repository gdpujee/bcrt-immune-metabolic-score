"""Build the Supplementary Information PDF (SUPPPDF-001).

BCRT requires the supplementary material as a PDF, and requires every
supplementary file to carry the article title, the journal name, the author names
and the corresponding author's affiliation and e-mail address.  The previous
supplement was an internal engineering document (project codename, revision tag,
evidence-ledger reference, run ids, hashes, script paths), which is not something
a reader or an editor should receive.

This renders submission_bcrt/supplement.md — the source for the supporting methods,
Table S1 and figures — and appends the generated REMARK checklist from
submission_bcrt/REMARK_checklist.md into one standalone PDF with:

  * a header block carrying title / journal / authors / corresponding author,
  * the S1-S7 sections, with markdown tables and bullets,
  * the supplementary figures embedded beside their captions,
  * the completed 20-item REMARK checklist as the final section.

Byte-reproducibility: reportlab stamps the wall-clock time into the document
metadata, so SOURCE_DATE_EPOCH is pinned before any canvas exists (same technique
as scripts/25_build_final_pdf.py).  The build asserts that the four supplementary
figures are present and that the header block is complete, because a supplement
that silently loses its figures or its authorship line is exactly the defect this
script exists to prevent.

Outputs: submission_bcrt/ESM_1.pdf, logs/supplement_pdf.log
"""
import os
import re
from pathlib import Path

import importlib.util
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (HRFlowable, Image, KeepTogether,
                                Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

os.environ.setdefault("SOURCE_DATE_EPOCH", "1700000000")

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / "submission_bcrt"
MD = SUB / "supplement.md"
REMARK = SUB / "REMARK_checklist.md"
OUT = SUB / "ESM_1.pdf"
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
logf = open(LOGS / "supplement_pdf.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


# --------------------------------------------------------------------- text
SUP = {"\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3", "\u2074": "4",
       "\u2075": "5", "\u2076": "6", "\u2077": "7", "\u2078": "8", "\u2079": "9",
       "\u207b": "-"}


def esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def rich(t):
    """Inline markdown -> reportlab mini-HTML.

    Order matters: the escaped-asterisk placeholder is installed before the bold
    and italic passes so that a literal `\\*` (the corresponding-author marker)
    cannot be mistaken for the start of an italic run.
    """
    t = esc(t)
    t = t.replace("\\*", "\x00STAR\x00")
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<i>\1</i>", t)
    t = re.sub(r"`([^`\n]+)`", r"\1", t)

    # scientific notation and superscripts
    t = re.sub(r"([0-9.]+)(?:&times;|\u00d7)10([\u207b\u2070-\u2079]+)",
               lambda m: r"%s&times;10<sup>%s</sup>"
               % (m.group(1), "".join(SUP.get(c, c) for c in m.group(2))), t)
    # e-notation as compact unicode superscripts: an HTML <sup> exponent is a
    # line-break opportunity, which let "1.36x10-18" render as "1.36x10-1 / 8" in
    # the narrow LRT-p column (review round 6 item 4).
    _sup = "0123456789"
    _supmap = str.maketrans(_sup, "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")
    t = re.sub(r"\b([0-9.]+)e-([0-9]+)\b",
               lambda m: f"{m.group(1)}\u00d710\u207b{m.group(2).translate(_supmap)}", t)

    for a, b in (("\u03c7\u00b2", "&chi;<sup>2</sup>"), ("\u03c7", "&chi;"),
                 ("\u0394C", "&Delta;C"), ("\u0394", "&Delta;"),
                 ("\u00d7", "&times;"), ("\u00b1", "&plusmn;"),
                 ("\u2265", "&ge;"), ("\u2264", "&le;"),
                 ("\u2192", "&rarr;"), ("\u00b7", "&middot;"),
                 ("\u2212", "-"), ("\u2013", "&ndash;"), ("\u2014", "&mdash;")):
        t = t.replace(a, b)
    return t.replace("\x00STAR\x00", "*")


title_st = ParagraphStyle("title", fontName="DVS", fontSize=17, leading=21,
                          alignment=TA_CENTER, spaceAfter=10)
arttitle_st = ParagraphStyle("arttitle", fontName="DVS", fontSize=12, leading=16,
                             alignment=TA_CENTER, spaceAfter=6)
journal_st = ParagraphStyle("journal", fontName="DVS", fontSize=11, leading=15,
                            alignment=TA_CENTER, spaceAfter=6)
auth_st = ParagraphStyle("auth", fontName="DVS", fontSize=11, leading=15,
                         alignment=TA_CENTER, spaceAfter=4)
aff_st = ParagraphStyle("aff", fontName="DVS", fontSize=9, leading=12.5,
                        alignment=TA_CENTER, spaceAfter=3)
body = ParagraphStyle("body", fontName="DVS", fontSize=9.5, leading=13.5,
                      alignment=TA_JUSTIFY, spaceAfter=6)
bullet = ParagraphStyle("bullet", fontName="DVS", fontSize=9.5, leading=13.5,
                        alignment=TA_JUSTIFY, spaceAfter=4, leftIndent=10,
                        bulletIndent=0)
sub_bullet = ParagraphStyle("sub_bullet", fontName="DVS", fontSize=9, leading=12.5,
                            alignment=TA_LEFT, spaceAfter=2, leftIndent=24,
                            bulletIndent=12)
h1 = ParagraphStyle("h1", fontName="DVS", fontSize=13, leading=17,
                    spaceBefore=14, spaceAfter=6)
h3 = ParagraphStyle("h3", fontName="DVS", fontSize=10.5, leading=14,
                    spaceBefore=10, spaceAfter=4)
cap = ParagraphStyle("cap", fontName="DVS", fontSize=9, leading=12.5,
                     alignment=TA_JUSTIFY, spaceAfter=6)
cell = ParagraphStyle("cell", fontName="DVS", fontSize=6.5, leading=8.5,
                      alignment=TA_CENTER)
cell_l = ParagraphStyle("cell_l", fontName="DVS", fontSize=6.5, leading=8.5,
                        alignment=TA_LEFT)
cell_h = ParagraphStyle("cell_h", fontName="DVS", fontSize=6.5, leading=8.5,
                        alignment=TA_CENTER)

AVAIL = 170 * mm


def column_widths(rows):
    """Proportional column widths, floored so a short numeric column still fits.

    A purely proportional split starves columns whose content is short but whose
    header is not ("deaths", "Cohort SD"), and the header then wraps mid-word —
    the same legibility defect the review raised for the Word tables.  Floor 9,
    cap 30: wide prose columns stop dominating, narrow ones stay readable.
    """
    w = []
    for c in range(len(rows[0])):
        longest = max(len(r[c]) for r in rows)
        w.append(min(max(longest, 11), 30))
    tot = sum(w)
    return [AVAIL * x / tot for x in w]


def make_table(rows):
    hdr = [Paragraph(f"<b>{rich(c)}</b>", cell_h) for c in rows[0]]
    body_rows = [[Paragraph(rich(r[j]), cell_l if j == 0 else cell)
                  for j in range(len(r))] for r in rows[1:]]
    t = Table([hdr] + body_rows, colWidths=column_widths(rows), repeatRows=1)
    t.setStyle(TableStyle([
        ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
        ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 2.5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2.5),
    ]))
    return t


def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("DVS", 8)
    canvas.setFillColor(colors.grey)
    canvas.drawString(20 * mm, 282 * mm,
                      "Supplementary Information \u2014 Breast Cancer Research "
                      "and Treatment")
    canvas.drawRightString(190 * mm, 15 * mm, "Page %d" % doc.page)
    canvas.restoreState()


# --------------------------------------------------------------------- parse
text = MD.read_text()
_title_match = re.search(r"(?m)^\*\*Article title:\*\*\s*(.+)$", text)
if not _title_match:
    raise SystemExit("supplement source has no article title")
article_title = _title_match.group(1).strip()

# A supplement that names the project codename, a revision tag, an internal
# ledger or a result path is the failure mode this script was written for, so the
# source is checked here rather than only in the submission gate.
TRACES = [r"scripts/\d", r"results/raw", r"logs/", r"\bbio-dsh\b",
          r"EVIDENCE_LEDGER", r"DECISIONS\.md", r"MASTER_CHECKLIST", r"STATE\.md",
          r"\bv3\.\d\b", r"\bW4\b", r"\bsha256\b", r"metadata/"]
for pat in TRACES:
    hit = re.search(pat, text)
    if hit:
        raise SystemExit(f"supplement.md carries an internal trace {hit.group(0)!r} "
                         "and must not be rendered as-is")

# The header block is a hard requirement (title / journal / authors / corresponding
# author), so its parts are asserted before anything is drawn.
REQUIRED = {
    "article title": article_title,
    "journal name": "Breast Cancer Research and Treatment",
    "first author": "Danhua He",
    "corresponding author": "Qiang Li",
    "corresponding e-mail": "gdpujee@gmail.com",
}
for what, needle in REQUIRED.items():
    if needle not in text:
        raise SystemExit(f"supplement header is missing the {what} ({needle!r})")

IMG_RE = re.compile(r"^!\[(?P<alt>[^\]]*)\]\((?P<path>[^)]+)\)$")
TABLE_SEP = re.compile(r"^\|[\s\-:|]+\|$")

story = []
lines = text.split("\n")
i = 0
n_img = 0
n_tbl = 0
seen_section = False
guard = 0
while i < len(lines):
    guard += 1
    assert guard < 100000, "loop guard tripped"
    s = lines[i].strip()
    if not s:
        story.append(Spacer(1, 3))
        i += 1
        continue
    if s.startswith("<!--"):
        i += 1
        continue

    m = IMG_RE.match(s)
    if m:
        path = (MD.parent / m.group("path")).resolve()
        if not path.exists():
            raise SystemExit(f"supplement references a missing image: {m.group('path')}")
        story.append(Image(str(path), width=150 * mm, height=98 * mm,
                           kind="proportional"))
        story.append(Spacer(1, 10))
        n_img += 1
        i += 1
        continue

    if s.startswith("|") and i + 1 < len(lines) and TABLE_SEP.match(lines[i + 1].strip()):
        rows = []
        while i < len(lines) and lines[i].strip().startswith("|"):
            raw = lines[i].strip().strip("|")
            if not TABLE_SEP.match(lines[i].strip()):
                rows.append([c.strip() for c in raw.split("|")])
            i += 1
        story.append(Spacer(1, 2))
        story.append(make_table(rows))
        story.append(Spacer(1, 8))
        n_tbl += 1
        continue

    if s.startswith("# ") and not s.startswith("## "):
        story.append(Paragraph(rich(s[2:]), title_st))
        i += 1
        continue
    if s.startswith("### "):
        story.append(Paragraph(rich(s[4:]), h3))
        i += 1
        continue
    if s.startswith("## "):
        if not seen_section:
            story.append(HRFlowable(width="100%", thickness=0.6, color=colors.grey,
                                    spaceBefore=6, spaceAfter=2))
            seen_section = True
        story.append(Paragraph(rich(s[3:]), h1))
        i += 1
        continue

    m = re.match(r"\*\*(Article title|Journal|Authors):\*\*\s*(.*)", s)
    if m:
        style = {"Article title": arttitle_st, "Journal": journal_st,
                 "Authors": auth_st}[m.group(1)]
        story.append(Paragraph(rich(m.group(2)), style))
        i += 1
        continue

    if s.startswith("  - "):
        story.append(Paragraph(rich(s[4:]), sub_bullet,
                               bulletText="\u2013"))
        i += 1
        continue
    if s.startswith("- "):
        story.append(Paragraph(rich(s[2:]), bullet, bulletText="\u2022"))
        i += 1
        continue

    # a plain paragraph; in the header block it is an affiliation or the
    # corresponding-author line, so it is centred rather than justified
    story.append(Paragraph(rich(s), body if seen_section else aff_st))
    i += 1

# The figures are the reason the supplement exists; losing one silently is the
# defect this asserts against.  Compare against the source, not a literal, so
# adding a supplementary figure cannot leave the check behind.
n_src_img = sum(1 for ln in lines if IMG_RE.match(ln.strip()))
if n_img != n_src_img:
    raise SystemExit(f"{n_img} image(s) embedded for {n_src_img} reference(s)")
if n_img < 4:
    raise SystemExit(f"only {n_img} supplementary figure(s) rendered; expected >= 4")

# The main text states that the completed REMARK checklist is supplementary
# material. Include the generated checklist in this PDF so the delivered item is
# a standard SI file rather than a Markdown-only appendix.
remark_rows = []
for ln in REMARK.read_text().splitlines():
    s = ln.strip()
    if TABLE_SEP.match(s):
        continue
    if s.startswith("|"):
        cells = [c.strip() for c in s.strip("|").split("|")]
        if len(cells) == 4:
            remark_rows.append(cells)
    elif remark_rows:
        break
if len(remark_rows) != 21:
    raise SystemExit("REMARK checklist must have a header and 20 rows; "
                     f"found {len(remark_rows)}")

story.append(Paragraph("REMARK checklist", h1))
story.append(Paragraph(
    "Completed checklist for tumor-marker prognostic studies. Item locations and "
    "evidence refer to the manuscript, tables and figures in this submission.", body))
remark_hdr = [Paragraph(f"<b>{rich(c)}</b>", cell_h) for c in remark_rows[0]]
remark_body = [[Paragraph(rich(c), cell_l) for c in row] for row in remark_rows[1:]]
remark_table = Table(
    [remark_hdr] + remark_body,
    colWidths=[12 * mm, 35 * mm, 38 * mm, 85 * mm], repeatRows=1)
remark_table.setStyle(TableStyle([
    ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
    ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("TOPPADDING", (0, 0), (-1, -1), 2.5),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ("LEFTPADDING", (0, 0), (-1, -1), 2.5),
    ("RIGHTPADDING", (0, 0), (-1, -1), 2.5),
]))
story.append(remark_table)
n_tbl += 1

# A trailing Spacer after the final figure pushed an empty page (review round 6
# item 4: blank page 7); BCRT publishes supplements as received, so drop it.
while story and type(story[-1]).__name__ == "Spacer":
    story.pop()

doc = SimpleDocTemplate(
    str(OUT), pagesize=A4, topMargin=25 * mm, bottomMargin=20 * mm,
    leftMargin=20 * mm, rightMargin=20 * mm,
    title=f"Supplementary Information - {article_title}",
    author="Danhua He, Qiang Li")
doc.build(story, onFirstPage=header_footer, onLaterPages=header_footer)
log(f"SUPPPDF-001: {OUT.name} with {n_img} figures and {n_tbl} tables; "
    f"{OUT.stat().st_size / 1024:.0f} KB")

# Verify the delivered artifact, not the intent: the header block must survive
# into page 1 of the rendered PDF.
from pypdf import PdfReader

pages = PdfReader(str(OUT)).pages
p1 = re.sub(r"\s+", " ", pages[0].extract_text() or "")
for what, needle in (("article title", article_title),
                     ("journal name", "Breast Cancer Research and Treatment"),
                     ("author names", "Danhua He"),
                     ("corresponding e-mail", "gdpujee@gmail.com")):
    if needle not in p1:
        raise SystemExit(f"page 1 of {OUT.name} does not carry the {what} "
                         f"({needle!r}); BCRT requires it on every supplementary file")
log(f"SUPPPDF-001: page 1 verified to carry title, journal, authors and "
    f"corresponding e-mail; {len(pages)} page(s)")
logf.close()
