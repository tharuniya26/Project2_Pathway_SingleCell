"""Run donor-level, covariate-adjusted differential expression with PyDESeq2.

The input matrix contains raw beta-cell pseudobulk counts (genes x donors).  The
model follows the covariate set reported by the GSE221156 study for glycaemic-
state comparisons: sequencing chemistry, sex, race, and BMI.  Donor is the unit
of replication; individual cells are never treated as independent samples.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cpus", type=int, default=2)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts_path = ROOT / "data" / "processed" / "beta_counts.csv"
    design_path = ROOT / "data" / "processed" / "beta_design.csv"
    output_dir = ROOT / "results" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(counts_path, index_col=0)
    metadata = pd.read_csv(design_path, index_col=0)
    metadata = metadata.loc[metadata["Condition"].isin(["ND", "T2D"])].copy()
    counts = counts.loc[:, metadata.index]

    if list(counts.columns) != list(metadata.index):
        raise ValueError("Count columns and donor metadata rows are not aligned")
    if counts.index.has_duplicates:
        raise ValueError("Development pseudobulk matrix has duplicate gene symbols")
    if (counts.to_numpy() < 0).any():
        raise ValueError("Counts must be non-negative")

    # Match the repository's pre-specified expression filter.
    keep = (counts >= 10).sum(axis=1) >= 5
    counts = counts.loc[keep]
    counts_for_deseq = counts.T.round().astype(np.int64)

    model_metadata = metadata[["Chemistry", "Sex", "Race", "BMI", "Condition"]].copy()
    model_metadata["Condition"] = pd.Categorical(
        model_metadata["Condition"], categories=["ND", "T2D"]
    )
    for column in ["Chemistry", "Sex", "Race"]:
        model_metadata[column] = model_metadata[column].astype("category")
    model_metadata["BMI"] = pd.to_numeric(model_metadata["BMI"], errors="raise")

    print(f"Analysing {counts_for_deseq.shape[0]} donors and "
          f"{counts_for_deseq.shape[1]} genes")
    print(model_metadata["Condition"].value_counts().sort_index())

    dds = DeseqDataSet(
        counts=counts_for_deseq,
        metadata=model_metadata,
        design="~ Chemistry + Sex + Race + BMI + Condition",
        refit_cooks=True,
        n_cpus=args.n_cpus,
        quiet=args.quiet,
    )
    dds.deseq2()

    stats = DeseqStats(
        dds,
        contrast=["Condition", "T2D", "ND"],
        alpha=0.05,
        n_cpus=args.n_cpus,
        quiet=args.quiet,
    )
    stats.summary()
    results = stats.results_df.copy()
    results.index.name = "gene"
    results = results.reset_index().sort_values("padj", na_position="last")

    output_path = output_dir / "beta_T2D_vs_ND_pydeseq2_adjusted.csv"
    results.to_csv(output_path, index=False)

    significant = results[results["padj"].lt(0.05)].copy()
    publication_threshold = significant[
        significant["log2FoldChange"].abs().ge(np.log2(1.5))
    ]
    summary = pd.DataFrame(
        [
            {"Measure": "Donors", "Value": len(model_metadata)},
            {"Measure": "ND donors", "Value": int((model_metadata.Condition == "ND").sum())},
            {"Measure": "T2D donors", "Value": int((model_metadata.Condition == "T2D").sum())},
            {"Measure": "Genes tested", "Value": len(results)},
            {"Measure": "FDR < 0.05", "Value": len(significant)},
            {
                "Measure": "FDR < 0.05 and abs(log2FC) >= log2(1.5)",
                "Value": len(publication_threshold),
            },
            {
                "Measure": "Up at FDR < 0.05 and effect threshold",
                "Value": int((publication_threshold.log2FoldChange > 0).sum()),
            },
            {
                "Measure": "Down at FDR < 0.05 and effect threshold",
                "Value": int((publication_threshold.log2FoldChange < 0).sum()),
            },
        ]
    )
    summary.to_csv(output_dir / "adjusted_de_summary.csv", index=False)

    print(summary.to_string(index=False))
    print(f"Saved {output_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
