# Exact ranked-time Grambsch-Therneau test for the documented METABRIC
# complete-case adjusted model (MBADJ-ZPH-001); selection timing is not established:
# GSE20685 and SCAN-B already report adjusted-model cox.zph; METABRIC must too.
# The frame is built by scripts/52_metabric_adj_frame.py, which asserts it
# reproduces the published complete-case n/deaths (1815/1041) before fitting.
library(survival)
library(jsonlite)

log_lines <- c()
say <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines <<- c(log_lines, m) }

d <- read.delim("metadata/METABRIC_adjusted_cc_frame.tsv", stringsAsFactors = FALSE)
stopifnot(nrow(d) == 1815, sum(d$OS_event) == 1041)

x <- d[, c("OS_years", "OS_event", "risk_sd", "age", "nodal_burden", "grade", "er", "her2", "tumor_size")]
fit <- coxph(Surv(OS_years, OS_event) ~ ., data = x, ties = "breslow", x = TRUE)
z <- cox.zph(fit, transform = "rank")
tab <- z$table

covs <- c("risk_sd", "age", "nodal_burden", "grade", "er", "her2", "tumor_size")
vars <- lapply(seq_along(covs), function(i) list(
  chisq = unname(tab[i, "chisq"]), p = unname(tab[i, "p"])))
names(vars) <- covs
out <- list(
  method = "survival::cox.zph (rank transform, Breslow ties) on the documented complete-case adjusted model; covariate-set selection timing not established",
  survival_version = as.character(packageVersion("survival")),
  unit = "risk_sd = per SD of the locked linear predictor, within-cohort SD",
  n = nrow(d), events = sum(d$OS_event),
  vars = vars,
  global = list(chisq = unname(tab[nrow(tab), "chisq"]),
                df = unname(tab[nrow(tab), "df"]),
                p = unname(tab[nrow(tab), "p"])))
say("METABRIC adjusted cox.zph: risk chisq=%.4f p=%.4g; global chisq=%.3f df=%d p=%.4g",
    tab["risk_sd", "chisq"], tab["risk_sd", "p"],
    tab[nrow(tab), "chisq"], tab[nrow(tab), "df"], tab[nrow(tab), "p"])
for (v in names(vars)) say("  %s: chisq=%.4f p=%.4g", v, vars[[v]]$chisq, vars[[v]]$p)

write(toJSON(out, auto_unbox = TRUE, digits = 10, pretty = TRUE),
      "results/raw/metabric_adjusted_zph.json")
writeLines(log_lines, "logs/metabric_adj_zph.log")
say("WROTE results/raw/metabric_adjusted_zph.json (MBADJ-ZPH-001)")
