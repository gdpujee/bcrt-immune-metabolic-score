#!/usr/bin/env Rscript
# Exact Grambsch-Therneau rank-time PH checks, independent of the Python screen.
library(survival)
library(jsonlite)

# Rscript --file=... includes the prefix in the option; resolve from the script path.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(script_arg) == 1)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
read_tsv <- function(path) read.delim(file.path(root, path), check.names = FALSE)
impute_numeric <- function(d, cols) {
  for (nm in cols) {
    v <- suppressWarnings(as.numeric(d[[nm]]))
    v[is.na(v)] <- median(v, na.rm = TRUE)
    d[[nm]] <- v
  }
  d[cols]
}
check <- function(time, event, x) {
  stopifnot(length(time) == length(event), nrow(x) == length(time), all(is.finite(time)), all(event %in% c(0, 1)))
  fit <- coxph(Surv(time, event) ~ ., data = x, ties = "breslow", x = TRUE)
  tab <- cox.zph(fit, transform = "rank")$table
  vars <- lapply(seq_len(ncol(x)), function(i) list(chisq = unname(tab[i, "chisq"]), p = unname(tab[i, "p"])))
  names(vars) <- colnames(x)
  list(n = length(time), events = sum(event), vars = vars,
       global = list(chisq = unname(tab[nrow(tab), "chisq"]),
                     df = unname(tab[nrow(tab), "df"]), p = unname(tab[nrow(tab), "p"])))
}
split_hr <- function(time, event, risk, cut) {
  # Keep everyone at risk until the cut; administratively censor survivors.
  early <- coxph(Surv(pmin(time, cut), as.integer(event == 1 & time <= cut)) ~ risk,
                 ties = "breslow")
  keep <- time > cut
  late <- coxph(Surv(rep(cut, sum(keep)), time[keep], event[keep]) ~ risk[keep],
                ties = "breslow")
  item <- function(fit, n, deaths) {
    ci <- exp(confint(fit))[1, ]
    list(n = n, deaths = deaths, HR = unname(exp(coef(fit))[1]),
         CI = unname(as.numeric(ci)), p = unname(summary(fit)$coefficients[1, "Pr(>|z|)"]))
  }
  list(early = item(early, length(time), sum(event == 1 & time <= cut)),
       late = item(late, sum(keep), sum(event[keep])))
}

g <- read_tsv("metadata/GSE20685_clinical_curated.tsv")
gr <- read_tsv("results/raw/validation_risk_GSE20685.tsv")
g <- merge(g, gr[, c("GSM", "risk")], by = "GSM")
g_risk <- g$risk / sd(g$risk)
g_cont <- check(g$follow_up_duration_years, g$event_death, data.frame(risk = g_risk))
g_adj <- check(g$follow_up_duration_years, g$event_death,
               impute_numeric(g, c("risk", "age_at_diagnosis", "t_stage", "n_stage")))
g_split5 <- split_hr(g$follow_up_duration_years, g$event_death, g_risk, 5)

s <- read_tsv("results/raw/rnaseq_risk_GSE96058.tsv")
s_cont <- check(s$OS_years, s$OS_event, data.frame(risk = s$risk / sd(s$risk)))
sc <- read_tsv("metadata/GSE96058_clinical_raw.tsv")
sc <- sc[!grepl("repl", sc$title), ]
sc <- merge(sc, s[, c("sample", "risk")], by.x = "title", by.y = "sample")
s_adj <- check(sc$overall_survival_days / 365.25, sc$overall_survival_event,
               impute_numeric(sc, c("risk", "age_at_diagnosis", "er_status", "her2_status")))

mc <- read_tsv("metadata/METABRIC_clinical_dl.tsv")
mc <- mc[!duplicated(mc$patientId), ]
mr <- read_tsv("results/raw/metabric_risk.tsv")
mc <- merge(mc, mr, by = "patientId")
m_cont <- check(as.numeric(mc$OS_MONTHS) / 12, as.integer(grepl("DECEASED", mc$OS_STATUS)),
                data.frame(risk = mc$risk / sd(mc$risk)))
mt <- as.numeric(mc$OS_MONTHS) / 12
me <- as.integer(grepl("DECEASED", mc$OS_STATUS))
mx <- mc$risk / sd(mc$risk)
m_splits <- lapply(c(3, 5, 7), function(cut) split_hr(mt, me, mx, cut))
names(m_splits) <- c("3", "5", "7")

out <- list(method = "survival::coxph ties=breslow; survival::cox.zph transform=rank",
            survival_version = as.character(packageVersion("survival")),
            GSE20685_cont = g_cont, GSE20685_adj = g_adj,
            SCANB_cont = s_cont, SCANB_adj = s_adj, METABRIC_cont = m_cont,
            GSE20685_split5y = g_split5, METABRIC_timesplit = m_splits)
write_json(out, file.path(root, "results/raw/ph_diagnostics_exact.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
for (nm in c("GSE20685_cont", "GSE20685_adj", "SCANB_cont", "SCANB_adj", "METABRIC_cont")) {
  x <- out[[nm]]
  cat(nm, "n=", x$n, "events=", x$events,
      "risk_p=", format(x$vars$risk$p, digits = 5),
      "global_p=", format(x$global$p, digits = 5), "\n")
}
