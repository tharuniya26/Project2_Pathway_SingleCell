"""Fail-fast validation of required inputs, results, and scientific invariants."""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def require(relative: str) -> Path:
    path = ROOT / relative
    if not path.exists() or path.stat().st_size == 0:
        raise AssertionError(f"Required non-empty file is missing: {relative}")
    return path


def main() -> None:
    required = [
        "data/processed/beta_counts.csv",
        "data/processed/beta_design.csv",
        "results/tables/beta_T2D_vs_ND_pydeseq2_adjusted.csv",
        "results/tables/adjusted_KEGG_overrepresentation.csv",
        "results/metrics/data_integrity_checks.csv",
        "results/metrics/internal_model_performance.csv",
        "results/metrics/internal_random_gene_set_performance.csv",
        "results/metrics/external_model_performance.csv",
        "results/metrics/external_random_gene_set_performance.csv",
        "results/metrics/pathway_vs_random_diagnostics.csv",
        "results/metrics/measured_covariate_associations.csv",
        "results/metrics/technical_covariate_model_performance.csv",
        "results/figures/discovery_pseudobulk_PCA.png",
        "results/figures/discovery_beta_cell_counts.png",
        "results/figures/internal_random_set_controls.png",
        "results/figures/comparative_model_performance.png",
        "results/figures/external_random_set_controls.png",
        "results/figures/external_comparative_probabilities.png",
        "results/figures/adjusted_T2D_vs_ND_volcano.png",
        "PROJECT_REPORT.md",
        "Tharuniya_Anpalakan_Project2_Research_Report.pdf",
        "README.md",
        "FINAL_VALIDATION_REPORT.md",
    ]
    for relative in required:
        require(relative)

    counts = pd.read_csv(ROOT / "data/processed/beta_counts.csv", index_col=0)
    design = pd.read_csv(ROOT / "data/processed/beta_design.csv", index_col=0)
    primary = design.loc[design.Condition.isin(["ND", "T2D"])]
    assert primary.Condition.value_counts().to_dict() == {"ND": 17, "T2D": 17}
    assert int(design.Ncells.sum()) == 99_029
    assert int(primary.Ncells.sum()) == 63_350
    assert list(counts.loc[:, primary.index].columns) == list(primary.index)
    assert counts.index.is_unique and primary.index.is_unique
    assert np.isfinite(counts.to_numpy()).all() and (counts.to_numpy() >= 0).all()

    integrity = pd.read_csv(ROOT / "results/metrics/data_integrity_checks.csv")
    assert integrity.Passed.astype(str).str.lower().eq("true").all()

    internal = pd.read_csv(ROOT / "results/metrics/internal_model_performance.csv")
    external = pd.read_csv(ROOT / "results/metrics/external_model_performance.csv")
    expected_models = {"Gene features", "Biological pathways", "Matched random gene sets"}
    assert set(internal.Model) == expected_models
    assert set(external.Model) == expected_models
    for frame in [internal, external]:
        for metric in ["Accuracy", "Balanced_accuracy", "ROC_AUC"]:
            assert frame[metric].between(0, 1).all()

    expected_auc = {
        "Gene features": (0.927336, 1.000000),
        "Biological pathways": (0.871972, 0.833333),
        "Matched random gene sets": (0.904983, 0.915833),
    }
    for model, (expected_internal, expected_external) in expected_auc.items():
        observed_internal = internal.loc[internal.Model.eq(model), "ROC_AUC"].item()
        observed_external = external.loc[external.Model.eq(model), "ROC_AUC"].item()
        assert np.isclose(observed_internal, expected_internal, atol=5e-7)
        assert np.isclose(observed_external, expected_external, atol=5e-7)

    random_internal = pd.read_csv(ROOT / "results/metrics/internal_random_gene_set_performance.csv")
    random_external = pd.read_csv(ROOT / "results/metrics/external_random_gene_set_performance.csv")
    assert len(random_internal) >= 50 and len(random_external) >= 50

    adjusted_de = pd.read_csv(ROOT / "results/tables/beta_T2D_vs_ND_pydeseq2_adjusted.csv")
    assert len(adjusted_de) == 17_408
    assert adjusted_de.gene.is_unique
    assert adjusted_de.padj.dropna().between(0, 1).all()

    external_predictions = pd.read_csv(
        ROOT / "results/metrics/external_comparative_predictions.csv"
    )
    assert len(external_predictions) == 10
    assert external_predictions.Condition.value_counts().to_dict() == {"ND": 6, "T2D": 4}
    external_design = pd.read_csv(ROOT / "data/processed/external_beta_design.csv")
    assert int(external_design.Ncells.sum()) == 270

    report_text = (ROOT / "PROJECT_REPORT.md").read_text(encoding="utf-8")
    readme_text = (ROOT / "README.md").read_text(encoding="utf-8")
    for rounded_value in ["0.927", "0.872", "0.905", "1.000", "0.833", "0.916"]:
        assert rounded_value in report_text and rounded_value in readme_text
    assert "central hypothesis was not supported" in report_text.lower()
    assert "did not outperform" in readme_text.lower()

    pdf = ROOT / "Tharuniya_Anpalakan_Project2_Research_Report.pdf"
    assert pdf.stat().st_size > 500_000
    assert pdf.read_bytes()[:5] == b"%PDF-"

    # Verify repository-relative Markdown links and embedded figures.
    markdown_link = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
    checked_links = 0
    for markdown_file in [ROOT / "README.md", ROOT / "PROJECT_REPORT.md", ROOT / "FINAL_VALIDATION_REPORT.md"]:
        for target in markdown_link.findall(markdown_file.read_text(encoding="utf-8")):
            clean_target = target.split("#", 1)[0]
            if not clean_target or clean_target.startswith(("http://", "https://", "mailto:")):
                continue
            linked = (markdown_file.parent / clean_target).resolve()
            if not linked.exists():
                raise AssertionError(f"Broken local link in {markdown_file.name}: {target}")
            checked_links += 1

    print(f"PASS: {len(required)} required artifacts found")
    print("PASS: discovery donor alignment, classes, 99,029 beta cells, counts, and identifiers")
    print("PASS: three-way internal and independent external comparisons")
    print("PASS: exact internal/external headline AUCs and negative hypothesis conclusion")
    print("PASS: at least 50 matched-random controls per cohort")
    print("PASS: adjusted differential-expression table")
    print("PASS: independent cohort contains 270 beta cells from 6 ND and 4 T2D donors")
    print("PASS: final PDF report signature and expected size")
    print(f"PASS: {checked_links} repository-relative Markdown links")


if __name__ == "__main__":
    main()
