# Pathway-Aware Single-Cell Analysis of Type 2 Diabetes

**Research project by Tharuniya Anpalakan · Human pancreatic islet beta cells · Computational transcriptomics**

## Research question

Do KEGG biological pathway representations improve donor-level type 2 diabetes (T2D) classification and transfer across cohorts compared with individual genes and size-matched random gene sets?

## Study design

- Discovery cohort: [GSE221156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156), human pancreatic-islet single-cell RNA-seq; donor-level beta-cell pseudobulk analysis.
- Independent external cohort: E-MTAB-5061, with cohort-specific preprocessing and aligned biological features.
- Comparisons: gene-level features, KEGG pathway features and matched random gene groups.
- Evaluation: internal cross-validation, external ROC AUC, and empirical tests of pathway superiority against random controls.

## Completed findings

**The pathway-superiority hypothesis was not supported.** Pathway features retained predictive information, but did not outperform genes or matched random groups in the completed comparison.

| Feature representation | Internal ROC AUC | External ROC AUC |
| --- | ---: | ---: |
| Individual genes | 0.927 | 1.000 |
| KEGG pathways | 0.872 | 0.833 |
| Matched random sets (mean) | 0.905 | 0.916 |

Pathway performance ranked at the **2nd percentile** of matched random controls internally and the **16th percentile** externally. Empirical tests for pathway superiority were not significant (**p = 0.980** internal; **p = 0.922** external).

**Interpretation:** Biologically named pathway features do not necessarily confer a predictive advantage. The result is a negative hypothesis test, not a failure of the research. Small cohorts, study differences and feature-construction choices limit generalisation. In particular, the external gene-level ROC AUC of 1.000 must not be interpreted as clinical validation.

## Research outputs and reproducibility

- [Existing IEEE-format manuscript PDF](paper/ieee_pathway_singlecell.pdf)
- [Editable manuscript source](paper/ieee_pathway_singlecell.tex)
- [Results tables](results/tables/) and [figures](results/figures/)
- [Analysis code](src/) and [notebooks](notebooks/)

**Version-consistency notice (5 November 2026):** The summary above reflects the later completed evaluation reported by the project author, including matched-random-set controls. The currently committed PDF and TeX manuscript document an earlier analysis and **must be regenerated and checked against the final local results before being described as the final paper**. The later result files have not been independently verified within this GitHub repository.

## Tools

Python; Scanpy; AnnData; pandas; NumPy; scikit-learn; matplotlib; Git and GitHub.

## Research status

Completed comparative analysis, with manuscript and repository artifact reconciliation pending. No peer-reviewed publication or clinical deployment is claimed. Raw and processed single-cell datasets are not redistributed in this repository.
