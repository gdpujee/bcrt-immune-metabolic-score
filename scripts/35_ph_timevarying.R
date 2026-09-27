#!/usr/bin/env Rscript
# Formal time-varying (non-proportional-hazards) analysis (PH-TVC-001).
#
# The external review asked for a formal treatment of the non-PH finding rather
# than only split-time sensitivity: a time-varying-coefficient Cox model
#   h(t | risk) = h0(t) * exp(beta1 * risk + beta2 * risk * log(t))
# so that HR(t) = exp(beta1 + beta2 * log t) is estimated with a proper
# covariance matrix.  Cohorts: GSE20685 (microarray), SCAN-B (RNA-seq),
# METABRIC (Illumina).  risk is the locked linear predictor scaled by each
# cohort's own SD (the same per-SD unit used in the primary analysis).
#
# Outputs: results/raw/ph_timevarying.json, figures/HR_t_curves.png/.pdf,
#          logs/ph_timevarying.log
library(survival)
library(jsonlite)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(script_arg) == 1)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
read_tsv <- function(path) read.delim(file.path(root, path), check.names = FALSE)
logcon <- file(file.path(root, "logs/ph_timevarying.log"), open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, logcon); flush(logcon) }

fit_tvc <- function(time, event, risk, label) {
  keep <- is.finite(time) & time > 0 & !is.na(risk)
  if (any(!keep)) say(sprintf("%s: excluded %d records with non-positive/NA follow-up time", label, sum(!keep)))
  time <- time[keep]; event <- event[keep]; risk <- risk[keep]
  stopifnot(all(is.finite(time)), all(time > 0), all(event %in% c(0, 1)))
  d <- data.frame(time = time, event = event, risk = risk)
  fit <- coxph(Surv(time, event) ~ risk + tt(risk),
               tt = function(x, t, ...) x * log(t), data = d, ties = "breslow", x = TRUE)
  fit0 <- coxph(Surv(time, event) ~ risk, data = d, ties = "breslow")
  b <- coef(fit); V <- vcov(fit)
  grid <- c(1, 2, 3, 5, 7, 10)
  hr <- exp(b[1] + b[2] * log(grid))
  se <- sqrt(V[1, 1] + log(grid)^2 * V[2, 2] + 2 * log(grid) * V[1, 2])
  # LRT of the time-varying term: model WITH tt vs the proportional-hazards model
  lrt <- 2 * (fit$loglik[2] - fit0$loglik[2])
  p_lrt <- pchisq(lrt, df = 1, lower.tail = FALSE)
  say(sprintf("%s: n=%d deaths=%d beta_risk=%.4f (p=%.3g) beta_risk:logt=%.4f (p=%.3g) LRT chi2=%.2f p=%.3g",
              label, nrow(d), sum(event), b[1], summary(fit)$coefficients[1, 5],
              b[2], summary(fit)$coefficients[2, 5], lrt, p_lrt))
  list(label = label, n = nrow(d), deaths = as.integer(sum(event)),
       beta_risk = unname(b[1]), beta_risk_p = unname(summary(fit)$coefficients[1, 5]),
       beta_tt = unname(b[2]), beta_tt_p = unname(summary(fit)$coefficients[2, 5]),
       LRT_chisq = unname(lrt), LRT_df = 1L, LRT_p = unname(p_lrt),
       hr_at = setNames(as.list(round(hr, 4)), paste0("t", grid)),
       ci_lo = setNames(as.list(round(exp(b[1] + b[2] * log(grid) - 1.96 * se), 4)), paste0("t", grid)),
       ci_hi = setNames(as.list(round(exp(b[1] + b[2] * log(grid) + 1.96 * se), 4)), paste0("t", grid)),
       tvc = list(fit = fit, grid = grid))
}

# --- GSE20685 (microarray) ---
g <- read_tsv("metadata/GSE20685_clinical_curated.tsv")
gr <- read_tsv("results/raw/validation_risk_GSE20685.tsv")
g <- merge(g, gr[, c("GSM", "risk")], by = "GSM")
G <- fit_tvc(g$follow_up_duration_years, g$event_death, g$risk / sd(g$risk), "GSE20685")

# --- SCAN-B (RNA-seq) ---
s <- read_tsv("results/raw/rnaseq_risk_GSE96058.tsv")
S <- fit_tvc(s$OS_years, s$OS_event, s$risk / sd(s$risk), "SCANB")

# --- METABRIC (Illumina) ---
mc <- read_tsv("metadata/METABRIC_clinical_dl.tsv")
mc <- mc[!duplicated(mc$patientId), ]
mr <- read_tsv("results/raw/metabric_risk.tsv")
mc <- merge(mc, mr, by = "patientId")
M <- fit_tvc(as.numeric(mc$OS_MONTHS) / 12,
             as.integer(grepl("DECEASED", mc$OS_STATUS)),
             mc$risk / sd(mc$risk), "METABRIC")

# --- HR(t) figure ---
tt <- seq(0.5, 12, length.out = 200)
curve_ci <- function(fit) {
  b <- coef(fit); V <- vcov(fit)
  lp <- b[1] + b[2] * log(tt)
  se <- sqrt(V[1, 1] + log(tt)^2 * V[2, 2] + 2 * log(tt) * V[1, 2])
  list(hr = exp(lp), lo = exp(lp - 1.96 * se), hi = exp(lp + 1.96 * se))
}
series <- list(
  list(nm = "GSE20685 (microarray, n=327)", col = "#1f77b4", ci = curve_ci(G$tvc$fit)),
  list(nm = "SCAN-B (RNA-seq, n=3273)", col = "#2ca02c", ci = curve_ci(S$tvc$fit)),
  list(nm = sprintf("METABRIC (Illumina, n=%d)", M$n), col = "#d62728", ci = curve_ci(M$tvc$fit))
)
draw <- function() {
  par(mar = c(4.2, 4.4, 1.2, 1.0), cex = 0.9)
  plot(NA, xlim = c(0.5, 12), ylim = c(0.55, 2.6), log = "xy",
       xlab = "Follow-up time (years)", ylab = "Hazard ratio per SD of locked score",
       axes = FALSE)
  axis(1, at = c(0.5, 1, 2, 3, 5, 7, 10, 12), labels = c(0.5, 1, 2, 3, 5, 7, 10, 12))
  axis(2, at = c(0.6, 0.8, 1, 1.5, 2, 2.5), labels = c(0.6, 0.8, 1, 1.5, 2, 2.5))
  box()
  abline(h = 1, lty = 3, col = "grey40")
  for (s in series) {
    polygon(c(tt, rev(tt)), c(s$ci$lo, rev(s$ci$hi)),
            col = adjustcolor(s$col, alpha.f = 0.15), border = NA)
    lines(tt, s$ci$hr, col = s$col, lwd = 2)
  }
  legend("topright", legend = vapply(series, function(s) s$nm, ""),
         col = vapply(series, function(s) s$col, ""), lwd = 2, bty = "n", cex = 0.8)
}
png(file.path(root, "figures/HR_t_curves.png"), width = 1500, height = 1050, res = 150)
draw(); dev.off()
pdf(file.path(root, "figures/HR_t_curves.pdf"), width = 7.2, height = 5.0)
draw(); dev.off()

out <- list(
  method = "survival::coxph with tt(risk) = risk*log(t) (time-varying coefficient, Breslow ties)",
  survival_version = as.character(packageVersion("survival")),
  unit = "per SD of locked linear predictor, within-cohort SD",
  GSE20685 = G[c("label", "n", "deaths", "beta_risk", "beta_risk_p", "beta_tt", "beta_tt_p", "LRT_chisq", "LRT_df", "LRT_p", "hr_at", "ci_lo", "ci_hi")],
  SCANB = S[c("label", "n", "deaths", "beta_risk", "beta_risk_p", "beta_tt", "beta_tt_p", "LRT_chisq", "LRT_df", "LRT_p", "hr_at", "ci_lo", "ci_hi")],
  METABRIC = M[c("label", "n", "deaths", "beta_risk", "beta_risk_p", "beta_tt", "beta_tt_p", "LRT_chisq", "LRT_df", "LRT_p", "hr_at", "ci_lo", "ci_hi")]
)
write_json(out, file.path(root, "results/raw/ph_timevarying.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 10)
say("WROTE results/raw/ph_timevarying.json + figures/HR_t_curves.png/.pdf (PH-TVC-001)")
close(logcon)
