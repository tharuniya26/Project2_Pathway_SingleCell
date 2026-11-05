# Final validation report

- **Project:** Pathway-Aware Single-Cell Analysis of Type 2 Diabetes
- **Validation date:** 5 November 2026
- **Author:** Tharuniya Anpalakan

## Overall status

**PASS with documented scientific limitations.** The canonical donor-level pipeline executes on the available real datasets, produces the required three-way representation comparison, performs genuinely independent external validation, and generates machine-readable results and figures. No commit or push was performed.

## Execution record

| Component | Status | Evidence |
|---|---|---|
| Workspace and Git audit | Pass | Existing history/results inspected; pre-existing manuscript edits preserved |
| GSE221156 metadata validation | Pass | 54 GEO library records, 48 donor library IDs, all 34 primary donors represented |
| Discovery pseudobulk integrity | Pass | 17 ND + 17 T2D; exact count/metadata alignment; unique genes/donors; non-negative integer counts |
| Donor-level QC | Pass | 63,350 beta cells represented; donor QC and cohort summaries generated |
| Adjusted differential expression | Pass | 17,408 genes tested; 927 FDR-significant; 520 also pass effect threshold |
| KEGG enrichment | Pass | Tested-gene background used; separate directional BH correction; four significant pathway-direction tests |
| Internal model comparison | Pass | Identical 10x repeated five-fold donor splits for genes, pathways, and random controls |
| Matched random controls | Pass | 50 size- and expression-decile-matched representations completed |
| External cohort processing | Pass | 270 metadata-annotated beta cells aggregated to 6 ND + 4 T2D donor pseudobulks |
| Independent validation | Pass | Fixed discovery models applied without retraining; gene/pathway/random comparisons generated |
| Uncertainty | Pass | 2,000 stratified donor bootstraps plus randomization distributions |
| Figures and tables | Pass | Required non-empty outputs generated under `results/figures`, `results/tables`, and `results/metrics` |
| Final PDF report | Pass | Built directly from `PROJECT_REPORT.md`; typography, figures, tables, references, page numbering, metadata, and extracted text inspected |
| Reproducible environment | Pass | Python 3.12 requirements and Conda environment supplied; pandas compatibility pin tested |
| Automated validation | Pass | `src/validate_project.py` checks inputs, outputs, metrics, cohorts, random replicate counts, and local Markdown links |

## Commands executed successfully

```powershell
.\.venv\Scripts\python.exe src\run_adjusted_differential_expression.py --n-cpus 1 --quiet
.\.venv\Scripts\python.exe src\run_comparative_analysis.py --cv-repeats 10 --random-replicates 50 --bootstrap-replicates 2000
.\.venv\Scripts\python.exe -m py_compile src\run_adjusted_differential_expression.py src\run_comparative_analysis.py src\validate_project.py
.\.venv\Scripts\python.exe src\build_report_pdf.py
.\.venv\Scripts\python.exe src\validate_project.py
.\.venv\Scripts\python.exe -m pip check
```

The first PyDESeq2 attempt exposed an incompatibility with pandas 3.0 during dispersion fitting. Pinning pandas 2.3.3 fixed the underlying dtype error; the complete 17,408-gene adjusted model then finished successfully.

## Validated numerical findings

- Adjusted DE: 927 genes at FDR <0.05; 321 higher and 199 lower in T2D after the additional |log2FC| >= log2(1.5) threshold.
- Internal ROC AUC: genes 0.927, pathways 0.872, matched random sets 0.905 mean.
- External ROC AUC: genes 1.000, pathways 0.833, matched random sets 0.916 mean.
- Pathway-versus-random superiority tests: internal *p*=0.980; external *p*=0.922. The central pathway-superiority hypothesis is not supported.
- Significant adjusted KEGG results: neuroactive ligand-receptor and ECM-receptor interaction among upregulated genes; maturity-onset diabetes of the young and a distinct neuroactive subset among downregulated genes.

## Leakage and pseudoreplication audit

- Biological donor, not cell, is the unit for DE and prediction.
- All scaling is fitted inside internal training folds.
- The same donor splits are used across fixed representations.
- No disease-label-driven gene or pathway selection is performed.
- External labels are used only for final scoring, never feature construction, fitting, threshold selection, or calibration.
- Expression filtering and random-set expression matching are label-blind but use the complete discovery feature matrix; this cohort-wide unsupervised step is disclosed as a limitation rather than described as fold-specific.

## Files generated or updated for review

### Canonical code and environment

- `src/run_adjusted_differential_expression.py`
- `src/run_comparative_analysis.py`
- `src/validate_project.py`
- `src/build_report_pdf.py`
- `run_all.ps1`
- `requirements.txt`
- `environment.yml`

### Documentation

- `README.md`
- `PROJECT_REPORT.md`
- `FINAL_VALIDATION_REPORT.md`
- `PROGRESS_LOG.md`
- `Tharuniya_Anpalakan_Project2_Research_Report.pdf`

### Principal new results

- `results/tables/beta_T2D_vs_ND_pydeseq2_adjusted.csv`
- `results/tables/adjusted_KEGG_overrepresentation.csv`
- `results/tables/discovery_cohort_summary.csv`
- `results/tables/discovery_donor_qc.csv`
- `results/tables/discovery_pathway_definitions.csv`
- `results/tables/external_matched_pathway_definitions.csv`
- all files under `results/metrics/`
- new adjusted-QC, enrichment, random-control, and comparative figures under `results/figures/`

## Items not treated as canonical

The three notebooks are preserved as the historical exploratory record. They contain duplicated/out-of-order cells and environment-specific online calls; they were inspected but are not the clean execution interface. Their scientific stages have been reimplemented and executed in deterministic scripts. The earlier IEEE manuscript draft under `paper/` predates the matched-random and comparative external analyses; it must not be presented as the final account of this repository. `PROJECT_REPORT.md` is the authoritative report.

## Remaining scientific limitations

1. Only 34 discovery and ten external donors are available.
2. External AUC=1.0 for genes is based on 24 cross-class pairs and has high unmeasured between-cohort uncertainty.
3. BMI and beta-cell yield differ between discovery groups; residual confounding is possible.
4. Discovery and external protocols differ (10x versus Smart-seq2), and probabilities are not cross-study calibrated.
5. Mean pathway scores ignore topology and opposing gene effects.
6. KEGG annotations overlap and are incomplete.
7. No prospective, clinical, or multi-cohort validation was performed.

## Git safety and publication readiness

Large raw/processed data and caches remain ignored. The release scope is limited to the canonical scripts, environments, validated result artifacts, reports, and final PDF. Pre-existing working-tree edits under `paper/`, the modified historical `results/figures/external_validation_probabilities.png`, `paper/page3.png`, and `src/polish_external_validation_figure.py` are unrelated to the canonical release and are intentionally excluded from staging. No history rewrite or force-push is required.
