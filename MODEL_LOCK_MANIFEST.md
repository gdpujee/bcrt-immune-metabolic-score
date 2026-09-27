# Model Lock Manifest: 14-Gene Immune–Metabolic Prognostic Score

**Document Version:** 1.0  
**Repository:** `gdpujee/bcrt-immune-metabolic-score`  
**Purpose:** Audit trail documenting the frozen model specification, safeguards against post-hoc tuning, and the authors' chronology of external validation.

---

## 1. Model Specification and Architecture

- **Model Identifier:** 14-gene immune-metabolic signature (`v1.0 TRAIN-001`)
- **Derivation Cohort:** GEO GSE42568 (Affymetrix HG-U133_Plus_2 GPL570; $n=104$ breast tumors, 35 deaths; 17 adjacent normal controls)
- **Candidate Pool Construction:** 
  - Differential expression: Mann–Whitney U test, Benjamini–Hochberg FDR $< 0.05$, relaxed $|\log_2\text{FC}| > 0.5$ (3,763 genes at strict $|\log_2\text{FC}| > 1$)
  - Pathway intersection: 11 pre-specified KEGG metabolic and immune pathways $\to$ 297 candidate genes
  - Univariable Cox pre-filter: $p < 0.01 \to 41$ genes
- **Feature Selection:** LASSO-penalized Cox regression with 5-fold cross-validation ($\alpha = 0.02976$) yielding 14 genes.
- **Multivariable Refitting:** Coefficients refit via unpenalized multivariable Cox regression in derivation cohort and strictly locked.

---

## 2. Locked Parameters and Scoring Formula

The individual prognostic risk score is calculated as a standardized linear combination:

$$\text{Score} = \sum_{i=1}^{14} \beta_i \left( \frac{X_i - \mu_i}{\sigma_i} \right)$$

where $X_i$ is the $\log_2$ gene expression level, $\mu_i$ is the derivation-cohort mean, $\sigma_i$ is the derivation-cohort standard deviation, and $\beta_i$ is the locked multivariable Cox regression coefficient.

### Complete Parameter Registry

| Gene Symbol | Functional Group | Locked Coefficient ($\beta_i$) | Training Mean ($\mu_i$) | Training SD ($\sigma_i$) |
| :--- | :--- | :--- | :--- | :--- |
| **FBP1** | Metabolic | -0.35971662232914275 | 7.306464180826923 | 1.8702794199539863 |
| **GRK6** | Immune | 0.15927842017141208 | 6.615874547894230 | 0.6886933262015541 |
| **IKBKB** | Immune | -0.09047055528507283 | 8.527793790769230 | 1.1025878816264900 |
| **RAP1B** | Immune | -0.42630761245035115 | 9.747929494701928 | 0.7152122838863929 |
| **TNFRSF19** | Immune | -0.19951373950008247 | 4.9963760304807705 | 1.3432194772168748 |
| **ALDH2** | Metabolic | -0.44352487720438130 | 8.686792377490380 | 1.3302826365700170 |
| **CXCL14** | Immune | -0.13224840870206375 | 9.117919255701922 | 2.5671127707074533 |
| **DLG1** | Immune | -0.04899088758068194 | 8.512231564673078 | 0.5555295846150139 |
| **IDNK** | Metabolic | -0.21195787820553647 | 4.768315584057688 | 0.8227711160479159 |
| **TNFRSF12A** | Immune | 0.14443631662641723 | 7.0575305361538465 | 1.2122786410673785 |
| **PRPS2** | Metabolic | -0.46153378939885180 | 7.560836554778846 | 0.9822798174685683 |
| **TKT** | Metabolic | 0.41755920456427910 | 8.585334860980772 | 0.9360152373121320 |
| **UQCRHL** | Metabolic | 0.52941717709690640 | 10.477828195528842 | 0.8731032606269039 |
| **PLCG1** | Immune | 0.57804157955870340 | 7.849202317663456 | 0.6637494992764560 |

### Locked Classification Cutoff
- **Risk Cutoff:** `-0.2880220748019792` (exact derivation median)
- **High Risk:** $\text{Score} > -0.2880220748019792$
- **Low Risk:** $\text{Score} \le -0.2880220748019792$

---

## 3. Cryptographic Verification

- **Artifact File:** `results/derived/locked_model.json`
- **SHA-256 Checksum:** `e8f2e00c6a9056998dcaa6d5b7f888e6db3119ea252a5a52089215b2e5917baf`
- **Verification Command:**
  ```bash
  shasum -a 256 results/derived/locked_model.json
  ```

---

## 4. Sequential Validation Timeline & Git Audit Trail

| Milestone | Git Commit Hash | Timestamp (UTC+8) | Description & Scope |
| :--- | :--- | :--- | :--- |
| **Model Lock** | `1b4c2c0` | 2026-09-20 01:07:21 | Model parameters, standardization parameters, and cutoff frozen. Initial validation on same-platform GSE20685 ($n=327$). |
| **SCAN-B RNA-seq Validation** | `88972bc` | 2026-09-20 06:38:50 | Applied locked model without re-training to SCAN-B GSE96058 ($n=3,273$). Identified PAM50 heterogeneity. |
| **METABRIC Extension** | `d50651b` | 2026-09-26 07:31:10 | Applied locked model without re-training to METABRIC ($n=1,980$). Tested PAM50 interaction replication. |

**Public verifiability status (updated 2026-09-27).** The three hashes above
come from the development repository's full history (347 commits) and are kept
as the **original local audit-trail timestamps**. That full development history
remains local to protect pre-publication drafting files, so those specific commit
hashes are not resolvable on the public remote and are not claimed to be.

What IS publicly checkable: the open code-and-data release at
`https://github.com/gdpujee/bcrt-immune-metabolic-score` (first commit `b1ff2b1`,
tagged `v1.0.1`) contains `results/derived/locked_model.json` with exactly the
SHA-256 recorded above (`e8f2e00c…917baf`, independently verifiable by anyone
who clones that commit). The public release thus verifies the frozen parameter
artifact and checksum; the original development timeline is documented by the
authors' local history.

---

## 5. Anti-Leakage & Governance Declaration

1. **No Post-Hoc Tuning:** No features were added or deleted, no coefficients refit, and no cutoffs readjusted in any validation cohort.
2. **Outcome-Blind Scoring:** Standardized score calculations in validation cohorts used only gene expression data and were executed blinded to clinical follow-up and survival endpoints.
3. **Full Disclosure of Limitations:** Failures of absolute cutoff transport (SCAN-B 4.77% high; METABRIC 96.41% high) and non-replication of PAM50 heterogeneity in METABRIC ($p=0.42$) are transparently reported without selective omission.
