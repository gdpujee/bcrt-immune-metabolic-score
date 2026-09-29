#!/usr/bin/env python3
"""Outcome-blind scoring-scale sensitivities for the locked 14-gene model.

Compares the recorded frozen-reference score with two explicitly retrospective
alternatives: cohort-wise gene z-scores and within-cohort inverse-normal ranks.
The latter two use each validation cohort's expression distribution and are not
single-sample deployment procedures. No outcome data are used to transform or
weight expression values.

Outputs: results/raw/score_scale_sensitivity.json, logs/score_scale_sensitivity.log
"""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
from statsmodels.duration.hazard_regression import PHReg

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
PROC = ROOT / "data/processed"
META = ROOT / "metadata"
DER = ROOT / "results/derived"
RES = ROOT / "results/raw"
LOGS = ROOT / "logs"


def log(message):
    print(message, flush=True)
    with open(LOGS / "score_scale_sensitivity.log", "a") as fh:
        fh.write(message + "\n")


def c_index(time, event, risk):
    comparable = concordant = 0.0
    for i in np.flatnonzero(event == 1):
        later = time > time[i]
        comparable += int(later.sum())
        concordant += int(np.sum(risk[later] < risk[i]))
        concordant += 0.5 * int(np.sum(risk[later] == risk[i]))
    return float(concordant / comparable) if comparable else 0.5


def align_cohort(name, expr, clinical, sample_col, time_col, event_col,
                 time_divisor=1.0):
    """Return genes x patients expression and aligned survival arrays."""
    clinical = clinical.copy()
    clinical[time_col] = pd.to_numeric(clinical[time_col], errors="coerce")
    clinical[event_col] = pd.to_numeric(clinical[event_col], errors="coerce")
    clinical = clinical.dropna(subset=[sample_col, time_col, event_col])
    clinical = clinical.drop_duplicates(sample_col, keep="first").set_index(sample_col)
    common = [sample for sample in expr.columns if sample in clinical.index]
    if not common:
        raise ValueError(f"{name}: no shared expression/clinical samples")
    x = expr.loc[:, common].astype(float)
    c = clinical.loc[common]
    t = c[time_col].to_numpy(dtype=float) / time_divisor
    e = c[event_col].to_numpy(dtype=int)
    if not set(np.unique(e)).issubset({0, 1}):
        raise ValueError(f"{name}: event field is not binary")
    return x, t, e


def load_scanb(genes):
    path = RAW / "GSE96058_gene_expression.csv.gz"
    with gzip.open(path, "rt", newline="") as fh:
        header = next(fh).rstrip("\n\r").split(",")
        titles = [x.strip('"') for x in header[1:]]
        keep = [i for i, title in enumerate(titles) if "repl" not in title]
        titles = [titles[i] for i in keep]
        rows = {}
        for line in fh:
            gene = line.split(",", 1)[0].strip().strip('"')
            if gene not in genes:
                continue
            fields = line.rstrip("\n\r").split(",")
            rows[gene] = [float(fields[i + 1]) for i in keep]
    missing = [g for g in genes if g not in rows]
    if missing:
        raise ValueError(f"SCAN-B missing locked genes: {missing}")
    expr = pd.DataFrame([rows[g] for g in genes], index=genes, columns=titles)
    clin = pd.read_csv(META / "GSE96058_clinical_raw.tsv", sep="\t")
    clin = clin[~clin.title.astype(str).str.contains("repl", na=False)]
    clin["OS_years"] = pd.to_numeric(clin.overall_survival_days, errors="coerce") / 365.25
    clin["OS_event"] = pd.to_numeric(clin.overall_survival_event, errors="coerce")
    return align_cohort("SCAN-B", expr, clin, "title", "OS_years", "OS_event")


def load_metabric(genes):
    expr = pd.read_csv(RAW / "METABRIC_mrna_14g.tsv", sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    missing = [g for g in genes if g not in expr.index]
    if missing:
        raise ValueError(f"METABRIC missing locked genes: {missing}")
    clin = pd.read_csv(META / "METABRIC_clinical_dl.tsv", sep="\t")
    clin["OS_years"] = pd.to_numeric(clin.OS_MONTHS, errors="coerce") / 12.0
    clin["OS_event"] = clin.OS_STATUS.astype(str).str.contains("DECEASED", na=False).astype(int)
    return align_cohort("METABRIC", expr.loc[genes], clin, "patientId", "OS_years", "OS_event")


def summarize(name, x, time, event, genes, beta, means, sds, cutoff):
    if x.index.tolist() != genes:
        x = x.loc[genes]
    if not np.isfinite(x.to_numpy(dtype=float)).all():
        raise ValueError(f"{name}: non-finite expression values")
    values = x.T
    frozen = ((values - means) / sds).to_numpy() @ beta
    cohort_z = ((values - values.mean(axis=0)) / values.std(axis=0, ddof=0)).to_numpy() @ beta
    ranks = rankdata(values.to_numpy(), axis=0, method="average")
    rank_normal = norm.ppf((ranks - 0.5) / len(values)) @ beta

    scores = {
        "frozen_reference": frozen,
        "cohort_gene_z": cohort_z,
        "cohort_inverse_normal_rank": rank_normal,
    }
    output = {"n": int(len(time)), "events": int(event.sum()), "scores": {}}
    for method, score in scores.items():
        sd = float(np.std(score, ddof=0))
        if not np.isfinite(sd) or sd == 0:
            raise ValueError(f"{name}/{method}: score has zero or invalid SD")
        fit = PHReg(time, (score / sd).reshape(-1, 1), event).fit(disp=0)
        coef, se = float(fit.params[0]), float(fit.bse[0])
        output["scores"][method] = {
            "score_sd": sd,
            "HR_per_SD": float(np.exp(coef)),
            "CI95": [float(np.exp(coef - 1.96 * se)), float(np.exp(coef + 1.96 * se))],
            "p": float(fit.pvalues[0]),
            "C_index": c_index(time, event, score),
        }
    output["score_correlations_vs_frozen"] = {
        method: float(np.corrcoef(frozen, score)[0, 1])
        for method, score in scores.items() if method != "frozen_reference"
    }
    above = frozen > cutoff
    output["frozen_cutoff"] = {
        "cutoff": float(cutoff), "high": int(above.sum()), "low": int((~above).sum())
    }
    if name == "METABRIC":
        binary_fit = PHReg(time, above.astype(int).reshape(-1, 1), event).fit(disp=0)
        coef, se = float(binary_fit.params[0]), float(binary_fit.bse[0])
        output["frozen_cutoff"].update({
            "HR": float(np.exp(coef)),
            "CI95": [float(np.exp(coef - 1.96 * se)), float(np.exp(coef + 1.96 * se))],
            "p": float(binary_fit.pvalues[0]),
        })
    return output


def main():
    (LOGS).mkdir(parents=True, exist_ok=True)
    (RES).mkdir(parents=True, exist_ok=True)
    (LOGS / "score_scale_sensitivity.log").write_text("")
    locked = json.load(open(DER / "locked_model.json"))
    genes = locked["genes"]
    beta = np.asarray(locked["coefs"], dtype=float)
    means = pd.Series(locked["scaling"]["means"]).loc[genes]
    sds = pd.Series(locked["scaling"]["sds"]).loc[genes]
    cutoff = float(locked["cutoff"])

    e206 = pd.read_csv(PROC / "GSE20685_expr_gene.tsv", sep="\t", index_col=0)
    c206 = pd.read_csv(META / "GSE20685_clinical_curated.tsv", sep="\t")
    c206 = c206[(c206.tissue == "primary breast cancer") & c206.OS_years.notna()
                & c206.OS_event.notna()].copy()
    x206 = e206.loc[genes]
    x206, t206, d206 = align_cohort(
        "GSE20685", x206, c206, "GSM", "OS_years", "OS_event")

    xscan, tscan, dscan = load_scanb(genes)
    xmeta, tmeta, dmeta = load_metabric(genes)

    cohorts = {
        "GSE20685": (x206, t206, d206),
        "SCAN-B": (xscan, tscan, dscan),
        "METABRIC": (xmeta, tmeta, dmeta),
    }
    results = {
        "analysis": "post hoc, outcome-blind scoring-scale sensitivity",
        "transformations": {
            "frozen_reference": "recorded derivation-cohort gene means/SDs, fixed coefficients",
            "cohort_gene_z": "within-cohort gene-wise mean/SD, fixed coefficients",
            "cohort_inverse_normal_rank": "within-cohort inverse-normal gene ranks, fixed coefficients",
        },
        "interpretation": "The latter two transformations use the full validation-cohort expression distribution and are retrospective sensitivities, not single-sample deployment procedures.",
        "cohorts": {},
    }
    for name, (x, t, d) in cohorts.items():
        results["cohorts"][name] = summarize(name, x, t, d, genes, beta, means, sds, cutoff)
        detail = results["cohorts"][name]
        log(f"{name}: n={detail['n']} events={detail['events']}")
        for method, res in detail["scores"].items():
            log(f"  {method}: HR/SD={res['HR_per_SD']:.3f} "
                f"[{res['CI95'][0]:.3f},{res['CI95'][1]:.3f}], "
                f"p={res['p']:.3g}, C={res['C_index']:.3f}")
        log(f"  frozen cutoff split={detail['frozen_cutoff']['high']}/"
            f"{detail['frozen_cutoff']['low']}")
    outpath = RES / "score_scale_sensitivity.json"
    outpath.write_text(json.dumps(results, indent=2) + "\n")
    log(f"WROTE {outpath.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
