"""Build BCRT LaTeX source from the submission manuscript (LATEXSRC-001).

Standard article class (compiles with any LaTeX engine; BCRT accepts LaTeX
source at submission). Includes full text, 4 tables (tabular from TSVs) and
14 figure environments referencing the committed PDF twins. Review rendering
only; no scientific content is created or altered here.
Outputs: manuscript/latex/main.tex, submission_bcrt/manuscript_bcrt.tex, logs/latex_src.log
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "submission_bcrt" / "manuscript_submission.md"
TAB = ROOT / "tables"
OUT_LATEX = ROOT / "manuscript" / "latex" / "main.tex"
OUT_BCRT = ROOT / "submission_bcrt" / "manuscript_bcrt.tex"
FIG = ROOT / "figures"
LOGS = ROOT / "logs"
logf = open(LOGS / "latex_src.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5",
       "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-"}


def tex(t):
    # LaTeX math replacements
    t = t.replace("χ²", r"$\chi^2$")
    t = t.replace("χ", r"$\chi$")
    t = t.replace("ΔC", r"$\Delta C$")
    t = t.replace("Δ", r"$\Delta$")
    t = t.replace("→", r"$\rightarrow$")
    t = t.replace("≥", r"$\geq$")
    t = t.replace("≤", r"$\leq$")
    t = t.replace("−", r"$-$")
    t = t.replace("∩", r"$\cap$")
    t = t.replace("∪", r"$\cup$")
    
    # Scientific notation: e.g. 8.3×10⁻¹⁴ -> $8.3\times 10^{-14}$
    t = re.sub(r"([0-9.]+)×10([⁻¹²³⁴⁵⁶⁷⁸⁹⁰]+)",
               lambda m: r"$%s\times 10^{%s}$" % (m.group(1), "".join(SUP.get(c, c) for c in m.group(2))), t)
    t = re.sub(r"×10([⁻¹²³⁴⁵⁶⁷⁸⁹⁰]+)",
               lambda m: r"$\times 10^{%s}$" % "".join(SUP.get(c, c) for c in m.group(1)), t)
    t = re.sub(r"\b([0-9.]+)e-([0-9]+)\b", r"$\1\\times 10^{-\2}$", t)
    
    t = t.replace("×", r"$\times$")
    t = t.replace("–", "--").replace("—", "---")
    t = t.replace("/", r"/\allowbreak ")
    t = re.sub(r"\|log2FC\|([><]=?)([0-9.]+)", r"$|\\mathrm{log}_2\\mathrm{FC}| \1 \2$", t)
    
    # Standard LaTeX escaping
    for a, b in [("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_")]:
        t = t.replace(a, b)
    t = re.sub(r"\*\*(.+?)\*\*", r"{\\bf \1}", t)

    # Greek letters and the middle dot go LAST, after the escaping loop: the LaTeX
    # they introduce contains "_", and escaping first would turn $\beta_1$ into
    # $\beta\_1$ and break math mode.  They were previously unmapped, so the engine
    # silently dropped them ("Missing character: There is no β (U+03B2) in font
    # lmroman10-regular") and the PDF shipped with the symbol missing.
    t = t.replace("β1", r"$\beta_1$").replace("β2", r"$\beta_2$")
    t = t.replace("β", r"$\beta$").replace("α", r"$\alpha$")
    t = t.replace("·", r"$\cdot$")
    return t


def section_lines(name):
    out, on = [], False
    for ln in MS.read_text().split("\n"):
        if ln.strip() == "## " + name:
            on = True
            continue
        if on and ln.startswith("## "):
            break
        if on and ln.strip():
            out.append(ln.strip())
    return out


FIG_MAP = {
    "Fig. 1": "Fig1_flow.pdf",
    "Fig. 2": "Fig2_volcano_42568.pdf",
    "Fig. 3": "Fig3_coef_forest_train.pdf",
    "Fig. 4": "train_KM.pdf",
    "Fig. 5": "Fig3_ROC_train.pdf",
    "Fig. 6": "valid_KM.pdf",
    "Fig. 7": "Fig4_ROC_valid.pdf",
    "Fig. 8": "Fig5_pathway_valid.pdf",
    "Fig. 9": "Fig5_checkpoints_valid.pdf",
    "Fig. 10": "RNASEQ_KM.pdf",
    "Fig. 11": "Fig_PAM50_forest.pdf",
    "Fig. 12": "Fig12_crosscohort_forest.pdf",
    "Fig. S1": "QC_followup_expr.pdf",
    "Fig. S2": "risk_distributions.pdf",
    "Fig. S3": "PH_schoenfeld.pdf",
    "Fig. S4": "HR_t_curves.pdf",
}

TBL_MAP = {
    "Table 1": "Tab1_cohorts_v3.tsv",
    "Table 2": "Tab2_coefficients.tsv",
    "Table 3": "Tab3_performance_v3.tsv",
    "Table 4": "Tab4_incremental_value.tsv",
}

TITLE_TEX = tex(" ".join(section_lines("Title")))
out = ["\\documentclass[11pt,a4paper]{article}",
       "\\usepackage[utf8]{inputenc}",
       "\\usepackage[english]{babel}",
       "\\usepackage{microtype}",
       "\\usepackage{graphicx,booktabs,hyperref,caption,amsmath,amssymb,enumitem}",
       "\\usepackage[margin=25mm]{geometry}",
       "\\title{%s}" % TITLE_TEX,
       "\\author{Danhua He$^{1}$, Qiang Li$^{2,*}$\\\\[2ex]"
       "\\small $^{1}$Guangdong Provincial Hospital of Chinese Medicine, Guangzhou, Guangdong, China\\\\"
       "\\small $^{2}$Foshan, Guangdong, China (Unaffiliated / Independent Researcher)\\\\"
       "\\small $^{*}$Corresponding author: Qiang Li (Email: gdpujee@gmail.com)}",
       "\\date{}",
       "\\begin{document}",
       "\\sloppy",
       "\\maketitle", ""]
nfig = ntbl = 0
skip_sec = False
for ln in MS.read_text().split("\n"):
    ln = ln.rstrip()
    if ln.startswith("# "):
        continue
    if ln.startswith("### "):
        out.append("\\subsection*{%s}" % tex(ln[4:]))
    elif ln.startswith("## "):
        skip_sec = ln.strip() in ("## Title", "## Authors and affiliations")
        if skip_sec:
            continue
        if ln.strip() == "## References":
            # Keep preceding figure/table floats ahead of the final bibliography.
            out.append("\\clearpage")
        out.append("\\section*{%s}" % tex(ln[3:]))
    elif skip_sec:
        continue
    elif ln.startswith("- **Fig"):
        m = re.match(r"- \*\*(Fig\.? ?S?\d+)\.\*\* (.*)", ln)
        if m:
            fig_id = m.group(1).replace("Fig ", "Fig. ")
            fig_cap = tex(m.group(2))
            fig_file = FIG_MAP.get(fig_id)
            if not fig_file:
                f_search = re.search(r"Files: ([A-Za-z0-9_]+)\.png", ln)
                if f_search:
                    fig_file = f_search.group(1) + ".pdf"
            
            out += ["\\begin{figure}[htbp]", "\\centering"]
            if not fig_file or not (FIG / fig_file).exists():
                raise SystemExit(f"no image for {fig_id} (FIG_MAP -> {fig_file!r}); a "
                                 "legend must not render as a caption with no figure")
            out.append("\\includegraphics[width=0.92\\textwidth]{../../figures/%s}" % fig_file)
            nfig += 1
            out += ["\\caption{%s}" % fig_cap, "\\end{figure}"]
    elif ln.startswith("- **Table"):
        m = re.match(r"- \*\*(Table \d+)\.\*\* (.*)", ln)
        if m:
            tbl_id = m.group(1)
            tbl_cap = tex(m.group(2))
            tbl_file = TBL_MAP.get(tbl_id)
            if not tbl_file:
                f_search = re.search(r"File: tables/([A-Za-z0-9_]+\.tsv)", ln)
                if f_search:
                    tbl_file = f_search.group(1)
            
            if tbl_file and (TAB / tbl_file).exists():
                import pandas as pd
                d = pd.read_csv(TAB / tbl_file, sep="\t", keep_default_na=False,
                                dtype=str).replace({"": "NA"})
                out.append("\\begin{table}[htbp]")
                out.append("\\centering")
                out.append("\\caption{%s}" % tbl_cap)
                
                if tbl_id == "Table 3":
                    # Table 3: Extract Note column into footnotes so table fits beautifully within textwidth
                    notes = []
                    if "Note" in d.columns:
                        for _, row in d.iterrows():
                            if row["Note"]:
                                notes.append(f"\\item \\textbf{{{tex(row['Cohort'])}}}: {tex(row['Note'])}")
                        d_disp = d.drop(columns=["Note"])
                    else:
                        d_disp = d
                    cols = "l" * len(d_disp.columns)
                    out.append("\\resizebox{\\textwidth}{!}{")
                    out.append("\\begin{tabular}{%s}\n\\toprule" % cols)
                    out.append(" & ".join(tex(c) for c in d_disp.columns) + " \\\\ \\midrule")
                    for _, r in d_disp.iterrows():
                        out.append(" & ".join(tex(x) for x in r.tolist()) + " \\\\")
                    out.append("\\bottomrule\n\\end{tabular}}")
                    if notes:
                        out.append("\\vspace{2mm}")
                        out.append("{\\footnotesize\\begin{itemize}[leftmargin=*]\\setlength\\itemsep{0.2em}")
                        out.extend(notes)
                        out.append("\\end{itemize}}")
                else:
                    cols = "l" * len(d.columns)
                    out.append("\\resizebox{\\textwidth}{!}{")
                    out.append("\\begin{tabular}{%s}\n\\toprule" % cols)
                    out.append(" & ".join(tex(c) for c in d.columns) + " \\\\ \\midrule")
                    for _, r in d.iterrows():
                        out.append(" & ".join(tex(x) for x in r.tolist()) + " \\\\")
                    out.append("\\bottomrule\n\\end{tabular}}")
                
                out.append("\\end{table}")
                ntbl += 1
    elif ln.startswith("- "):
        out.append("\\begin{itemize}\\item %s\\end{itemize}" % tex(ln[2:]))
    elif ln.strip() == "":
        out.append("")
    else:
        out.append(tex(ln) + "\n")
out.append("\\end{document}")

n_legends = sum(1 for l in MS.read_text().split("\n") if l.startswith("- **Fig"))
if nfig != n_legends:
    raise SystemExit(f"{nfig} image(s) embedded for {n_legends} figure legend(s)")
latex_content = "\n".join(out)

# The TeX engine reports an unmapped character as a warning, not an error, and then
# ships the PDF with the glyph missing.  Fail here instead: every symbol the source
# uses must have been translated to a LaTeX command by tex().
_odd = sorted({c for c in latex_content if ord(c) > 127})
if _odd:
    raise SystemExit("non-ASCII characters survived into the LaTeX source "
                     + ", ".join(f"{c!r} (U+{ord(c):04X})" for c in _odd)
                     + " — add them to tex() so the engine does not drop them")

OUT_LATEX.write_text(latex_content)
OUT_BCRT.write_text(latex_content)
log(f"LATEXSRC-001: {OUT_LATEX.name} and {OUT_BCRT.name} written with {nfig} figures + {ntbl} tables")
logf.close()
