"""Bootstrap selection stability for the 14-gene LASSO panel (BOOTSTAB-001).

The default run resamples the 104 derivation tumors (seed 42), repeats the
univariable Cox screen over the frozen 297-gene candidate pool, and fits the
L1 Cox model at the alpha recorded in the locked model. This estimates gene
selection stability conditional on that alpha; it does not repeat alpha tuning
or provide full-pipeline internal validation. Set BOOT_ALPHA_MODE=reselect to
repeat the 20-alpha, five-fold tuning grid in every resample, which is much
slower and separately budgeted.

Selection frequencies use only resamples that complete the selection. Samples
with too few events, too small a screen, fit failures, or timeouts are counted
as skipped, not as non-selections. A per-resample SIGALRM guard and a total
wall-clock budget bound the run. The JSON reports requested and completed B,
skipped counts, the locked alpha, and the actual denominator. The bootstrap is
post hoc, leaves the recorded model unchanged, and is not an optimism correction.

Outputs: results/raw/bootstrap_stability.json, figures/BOOT_stability.png(.pdf),
logs/bootstrap.log
"""

import json
import os
import signal
import sys
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.duration.hazard_regression import PHReg
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/raw"
DER = ROOT / "results/derived"
FIG = ROOT / "figures"
LOGS = ROOT / "logs"
B = int(sys.argv[1]) if len(sys.argv) > 1 else 500
SEED = 42
rng = np.random.default_rng(SEED)
ALPHAS = np.logspace(-3, 1, 20)
ALPHA_MODE = os.environ.get("BOOT_ALPHA_MODE", "locked")
if ALPHA_MODE not in {"locked", "reselect"}:
    raise SystemExit("BOOT_ALPHA_MODE must be 'locked' or 'reselect'")
RESAMPLE_GUARD_S = int(os.environ.get("BOOT_GUARD_S", "45" if ALPHA_MODE == "locked" else "300"))
BUDGET_S = int(os.environ.get("BOOT_BUDGET_S", "5400"))
PROGRESS_EVERY = int(os.environ.get("BOOT_PROGRESS_EVERY", "25"))

logf = open(LOGS / "bootstrap.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


class ResampleTimeout(Exception):
    pass


def _on_alarm(signum, frame):
    raise ResampleTimeout()


signal.signal(signal.SIGALRM, _on_alarm)


def guard_selftest():
    """Prove the per-resample time guard actually escapes the inner handlers.

    The first two runs of this script stalled for tens of minutes without ever
    logging a timeout or a progress mark, because the inner `except Exception`
    clauses wrapping `fit_regularized` swallowed ResampleTimeout — a timeout was
    silently recorded as a 0.5 fold score and the resample continued. A guard
    that cannot fire is worse than no guard, since it looks like protection.
    This exercises the same nesting, so a regression fails loudly in one second
    instead of costing an hour of wall clock.
    """
    swallowed = {"inner": False}

    def body():
        try:
            try:
                time.sleep(5)
            except ResampleTimeout:
                raise
            except Exception:
                swallowed["inner"] = True
        except ResampleTimeout:
            return "escaped"
        return "swallowed"

    signal.setitimer(signal.ITIMER_REAL, 1)
    try:
        outcome = body()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    if outcome != "escaped" or swallowed["inner"]:
        raise SystemExit("GUARD SELFTEST FAILED: ResampleTimeout does not escape "
                         "the inner handlers, so the time guard cannot fire")
    log("guard selftest: OK (ResampleTimeout escapes the inner handlers)")


guard_selftest()


def cindex(T, E, r):
    order = np.argsort(T, kind="stable")
    T, E, r = T[order], E[order], r[order]
    num = den = 0.0
    for i in np.where(E == 1)[0]:
        later = np.arange(i + 1, len(T))
        later = later[T[later] > T[i]]
        if len(later) == 0:
            continue
        den += len(later)
        num += np.sum(r[later] < r[i]) + 0.5 * np.sum(r[later] == r[i])
    return num / den if den else 0.5


e = pd.read_csv(ROOT / "data/processed/GSE42568_expr_gene.tsv", sep="\t", index_col=0)
clin = pd.read_csv(ROOT / "metadata/GSE42568_clinical_curated.tsv", sep="\t")
gsm = clin[(clin.tissue == "breast cancer") & clin.OS_time_days.notna()
           & clin.OS_event.notna()]["GSM"].tolist()
T0 = clin.set_index("GSM").loc[gsm, "OS_time_days"].values.astype(float) / 365.25
E0 = clin.set_index("GSM").loc[gsm, "OS_event"].values.astype(int)
pool = pd.read_csv(RES / "candidate_pool.tsv", sep="\t")["symbol"].tolist()
missing_pool = [g for g in pool if g not in e.index]
if missing_pool:
    raise SystemExit(f"candidate genes missing from expression matrix: {missing_pool[:10]}")
# Keep the matrix columns in the same order as `pool`. The expression file is
# ordered by mean expression, whereas candidate_pool.tsv is alphabetical; using
# e[gsm].T and then indexing columns by enumerate(pool) silently assigned
# bootstrap p-values and coefficients to the wrong genes.
Z0 = e.loc[pool, gsm].T
if Z0.columns.tolist() != pool:
    raise SystemExit("candidate-pool gene order does not match bootstrap matrix columns")
locked = json.load(open(DER / "locked_model.json"))
locked_genes = locked["genes"]
locked_alpha = float(locked["alpha"])
freq = {g: 0 for g in pool}
coefs = {g: [] for g in locked_genes}
signs = {g: [] for g in locked_genes}
nsel = []
skipped = {"few_events": 0, "small_screen": 0, "fit_failed": 0, "timeout": 0}

log(f"BOOTSTAB-001 B_requested={B} seed={SEED} alpha_mode={ALPHA_MODE} "
    f"locked_alpha={locked_alpha:.8g} guard={RESAMPLE_GUARD_S}s "
    f"budget={BUDGET_S}s pool={len(pool)} n={len(T0)} events={int(E0.sum())}")

t_start = time.time()
completed = 0
for b in range(B):
    if time.time() - t_start > BUDGET_S:
        log(f"wall-clock budget {BUDGET_S}s reached after {completed} completed "
            f"resample(s) at iteration {b + 1}/{B}; stopping")
        break
    t_iter = time.time()
    signal.setitimer(signal.ITIMER_REAL, RESAMPLE_GUARD_S)
    try:
        idx = rng.integers(0, len(T0), len(T0))
        T, E, Zb = T0[idx], E0[idx], Z0.values[idx]
        if E.sum() < 5:
            skipped["few_events"] += 1
            continue
        pv = {}
        for j, g in enumerate(pool):
            try:
                f = PHReg(T, Zb[:, j].reshape(-1, 1), E).fit(disp=0)
                pv[g] = float(f.pvalues[0])
            except ResampleTimeout:
                raise
            except Exception:
                pv[g] = 1.0
        sel = [g for g in pool if pv[g] < 0.01]
        if len(sel) < 2:
            skipped["small_screen"] += 1
            continue
        cols = [pool.index(g) for g in sel]
        Zs = (Zb[:, cols] - Zb[:, cols].mean(0)) / (Zb[:, cols].std(0) + 1e-9)
        kf = np.array_split(rng.permutation(len(T)), 5)
        if ALPHA_MODE == "reselect":
            best, bestc = None, -1
            for a in ALPHAS:
                cs = []
                for k in range(5):
                    te = kf[k]
                    tr = np.concatenate([kf[j] for j in range(5) if j != k])
                    try:
                        f = PHReg(T[tr], Zs[tr], E[tr]).fit_regularized(
                            method="elastic_net", alpha=a, L1_wt=1.0, maxiter=200)
                        cs.append(cindex(T[te], E[te], Zs[te] @ np.asarray(f.params).ravel()))
                    except ResampleTimeout:
                        # Must propagate: ResampleTimeout derives from Exception, so a
                        # bare `except Exception` here would swallow the time guard and
                        # silently record a 0.5 fold score instead of skipping the
                        # resample. That is exactly how the first two runs stalled for
                        # tens of minutes without ever logging a timeout.
                        raise
                    except Exception:
                        cs.append(0.5)
                if np.mean(cs) > bestc:
                    bestc, best = np.mean(cs), a
        else:
            # Post hoc conditional-stability sensitivity: hold the recorded alpha
            # fixed to the value selected in the original derivation cohort. This
            # avoids claiming that a small, budget-limited resampling run re-tuned
            # the complete selection pipeline in every bootstrap sample.
            best = locked_alpha
        fr = PHReg(T, Zs, E).fit_regularized(method="elastic_net", alpha=best,
                                             L1_wt=1.0, maxiter=500)
        params = np.asarray(fr.params).ravel()
        nz = [sel[j] for j, v in enumerate(params) if abs(v) > 1e-8]
    except ResampleTimeout:
        skipped["timeout"] += 1
        log(f"iteration {b + 1}: exceeded the {RESAMPLE_GUARD_S}s per-resample "
            f"guard, skipped")
        continue
    except Exception as exc:
        skipped["fit_failed"] += 1
        log(f"iteration {b + 1}: {type(exc).__name__}: {exc}; skipped")
        continue
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    completed += 1
    nsel.append(len(nz))
    for g in nz:
        freq[g] += 1
    for j, g in enumerate(sel):
        if g in coefs:
            v = float(params[j])
            if abs(v) > 1e-8:
                coefs[g].append(v)
                signs[g].append(int(np.sign(v)))
    if completed % PROGRESS_EVERY == 0:
        rate = (time.time() - t_start) / completed
        log(f"completed {completed}/{B} resample(s) in {time.time() - t_start:.0f}s "
            f"({rate:.1f}s each, {skipped})")

n_total = len(nsel)
if n_total == 0:
    raise SystemExit("no resample completed the selection — nothing to report")
log(f"completed {n_total}/{B} resample(s); denominator for all rates is "
    f"{n_total}; skipped={skipped}")
out = {"B_requested": B, "B_used": n_total, "alpha_mode": ALPHA_MODE,
       "locked_alpha": locked_alpha, "skipped": skipped,
       "resample_guard_s": RESAMPLE_GUARD_S, "budget_s": BUDGET_S,
       "median_nsel": float(np.median(nsel)),
       "genes": {}}
for g in locked_genes:
    out["genes"][g] = {"freq": round(freq[g] / n_total, 3),
                       "n_selected": freq[g],
                       "coef_median": round(float(np.median(coefs[g])), 4) if coefs[g] else None,
                       "sign_consistent": bool(len(set(signs[g])) == 1) if signs[g] else None}
    log(f"{g}: freq={freq[g] / n_total:.3f} n_coef={len(coefs[g])} "
        f"sign_consistent={out['genes'][g]['sign_consistent']}")
json.dump(out, open(RES / "bootstrap_stability.json", "w"), indent=2)
fig, ax = plt.subplots(figsize=(9, 4))
gs = sorted(locked_genes, key=lambda g: -freq[g] / n_total)
ax.bar(gs, [freq[g] / n_total for g in gs])
ax.axhline(0.5, ls="--", color="grey", lw=1)
ax.set_ylabel(f"selection frequency (B={n_total})")
ax.set_ylim(0, 1)
plt.xticks(rotation=30, ha="right", fontsize=8)
plt.tight_layout()
plt.savefig(FIG / "BOOT_stability.png", dpi=150)
plt.savefig(FIG / "BOOT_stability.pdf")
log("WROTE bootstrap_stability.json (BOOTSTAB-001)")
