# Progress log

## 5 November 2026

- Audited Git history, uncommitted manuscript work, datasets, notebooks, scripts, results, and environment.
- Preserved pre-existing manuscript/figure changes and excluded them from the canonical release commits.
- Verified GSE221156 metadata against GEO and its associated publication.
- Verified E-MTAB-5061 as a genuinely independent six-ND/four-T2D donor cohort.
- Installed PyDESeq2 0.5.3 and diagnosed its pandas 3.0 incompatibility; pinned pandas 2.3.3.
- Ran covariate-adjusted donor-level DE on 17,408 genes.
- Implemented and executed the missing 50-replicate matched-random gene-set control.
- Added fixed-fold internal comparison and comparative independent validation for genes, pathways, and random sets.
- Added donor bootstrap intervals, randomization diagnostics, QC/enrichment/model figures, reproducibility files, and automated validation.
- Replaced planned-workflow README language with completed findings and wrote the scientific report.
- Finalized the report against saved metrics and generated a visually inspected, 15-page A4 PDF directly from `PROJECT_REPORT.md`.
- Ran strengthened automated validation: 23 required artifacts, exact headline AUCs, all donor/model invariants, 50 random controls per cohort, the PDF signature/size, and 14 local Markdown links passed; `pip check` reported no broken requirements.
- Audited the release manifest, Git status, document links, PDF metadata/text/layout, credential patterns, and file sizes before staging.
- Published the validated release to the existing GitHub `main` branch using normal fast-forward pushes.
