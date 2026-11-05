"""Reproducible donor-level representation comparison and external validation.

This script addresses the project's central question using the same donors and
folds for all internal comparisons.  It evaluates individual genes, KEGG pathway
means, and expression/size-matched random gene sets.  The independent
E-MTAB-5061 cohort is processed at donor level and is never used for fitting.
"""

from __future__ import annotations

import argparse
import gzip
import importlib
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
METRICS = ROOT / "results" / "metrics"
FIGURES = ROOT / "results" / "figures"
CACHE = ROOT / "data" / "cache"
PROCESSED = ROOT / "data" / "processed"
EXTERNAL_RAW = ROOT / "data" / "raw" / "E-MTAB-5061"
RANDOM_SEED = 20261105


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cv-repeats", type=int, default=10)
    parser.add_argument("--random-replicates", type=int, default=50)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--skip-external", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for directory in [TABLES, METRICS, FIGURES, CACHE, PROCESSED]:
        directory.mkdir(parents=True, exist_ok=True)


def load_discovery() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(PROCESSED / "beta_counts.csv", index_col=0)
    design = pd.read_csv(PROCESSED / "beta_design.csv", index_col=0)
    design = design.loc[design["Condition"].isin(["ND", "T2D"])].copy()
    counts = counts.loc[:, design.index]

    checks = {
        "donor_alignment": list(counts.columns) == list(design.index),
        "unique_donor_ids": design.index.is_unique,
        "unique_gene_symbols": counts.index.is_unique,
        "counts_nonnegative": bool((counts.to_numpy() >= 0).all()),
        "counts_integer_valued": bool(np.allclose(counts.to_numpy(), counts.to_numpy().round())),
        "two_balanced_classes": design["Condition"].value_counts().to_dict()
        == {"ND": 17, "T2D": 17},
        "all_beta_pseudobulks_have_cells": bool((design["Ncells"] > 0).all()),
    }
    if not all(checks.values()):
        raise ValueError(f"Discovery-data validation failed: {checks}")

    keep = (counts >= 10).sum(axis=1) >= 5
    filtered = counts.loc[keep].copy()
    library_sizes = filtered.sum(axis=0)
    log_cpm = np.log2(filtered.divide(library_sizes, axis=1) * 1_000_000 + 1)

    pd.DataFrame(
        [{"Check": key, "Passed": value} for key, value in checks.items()]
    ).to_csv(METRICS / "data_integrity_checks.csv", index=False)

    donor_qc = design.copy()
    donor_qc["Filtered_library_size"] = library_sizes
    donor_qc["Detected_filtered_genes"] = (filtered > 0).sum(axis=0)
    donor_qc.index.name = "Donor"
    donor_qc.to_csv(TABLES / "discovery_donor_qc.csv")

    cohort_summary = (
        design.groupby("Condition", observed=True)
        .agg(
            Donors=("Condition", "size"),
            Beta_cells=("Ncells", "sum"),
            Median_beta_cells=("Ncells", "median"),
            Min_beta_cells=("Ncells", "min"),
            Max_beta_cells=("Ncells", "max"),
            Mean_age=("Age", "mean"),
            Mean_BMI=("BMI", "mean"),
            Mean_HbA1c=("HbA1c", "mean"),
        )
        .reset_index()
    )
    cohort_summary.to_csv(TABLES / "discovery_cohort_summary.csv", index=False)
    return filtered, design, log_cpm


def validate_geo_metadata(design: pd.DataFrame) -> None:
    path = ROOT / "data" / "metadata" / "GSE221156_series_matrix.txt.gz"
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    geo_line = next(line for line in lines if line.startswith("!Sample_geo_accession"))
    title_line = next(line for line in lines if line.startswith("!Sample_title"))
    geo_ids = [x.strip('"') for x in geo_line.rstrip().split("\t")[1:]]
    titles = [x.strip('"') for x in title_line.rstrip().split("\t")[1:]]
    records = [{"GSM": gsm, "sample_title": title} for gsm, title in zip(geo_ids, titles)]
    for line in lines:
        if line.startswith("!Sample_characteristics_ch1"):
            values = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
            for i, value in enumerate(values):
                if ": " in value:
                    key, val = value.split(": ", 1)
                    records[i][key.strip().lower().replace(" ", "_")] = val.strip()
    metadata = pd.DataFrame(records)
    metadata.to_csv(TABLES / "GSE221156_geo_sample_metadata.csv", index=False)

    expected = set(design.index)
    observed_titles = set(metadata["sample_title"])
    summary = pd.DataFrame(
        [
            {"Measure": "GEO library records", "Value": len(metadata)},
            {"Measure": "Unique sample titles", "Value": metadata.sample_title.nunique()},
            {"Measure": "Unique library IDs", "Value": metadata.library_id.nunique()},
            {"Measure": "Primary-analysis donors", "Value": len(design)},
            {"Measure": "Primary donors represented in GEO titles", "Value": len(expected & observed_titles)},
        ]
    )
    summary.to_csv(METRICS / "GSE221156_metadata_validation.csv", index=False)


def load_kegg() -> dict[str, list[str]]:
    """Load KEGG 2021 Human through Enrichr with a workspace-local cache."""
    cache_dir = CACHE / "gseapy"
    cache_dir.mkdir(parents=True, exist_ok=True)
    utils = importlib.import_module("gseapy.utils")
    enrichr = importlib.import_module("gseapy.enrichr")
    parser = importlib.import_module("gseapy.parser")
    for module in [utils, enrichr, parser]:
        module.DEFAULT_CACHE_PATH = str(cache_dir)
    sets = parser.get_library(name="KEGG_2021_Human", organism="human")
    return {name: sorted(set(map(str, genes))) for name, genes in sets.items()}


def filter_pathways(
    gene_sets: dict[str, list[str]], universe: Iterable[str], min_size: int = 10, max_size: int = 200
) -> dict[str, list[str]]:
    universe_set = set(map(str, universe))
    result = {}
    for name, genes in gene_sets.items():
        overlap = sorted(universe_set.intersection(genes))
        if min_size <= len(overlap) <= max_size:
            result[name] = overlap
    return result


def score_sets(log_cpm: pd.DataFrame, sets: dict[str, list[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        {name: log_cpm.loc[genes].mean(axis=0) for name, genes in sets.items()},
        index=log_cpm.columns,
    )


def expression_bins(log_cpm: pd.DataFrame, n_bins: int = 10) -> pd.Series:
    means = log_cpm.mean(axis=1).rank(method="first")
    return pd.qcut(means, q=n_bins, labels=False, duplicates="drop")


def matched_random_sets(
    template_sets: dict[str, list[str]], log_cpm: pd.DataFrame, rng: np.random.Generator
) -> dict[str, list[str]]:
    """Match every pathway's size and mean-expression-decile composition."""
    bins = expression_bins(log_cpm)
    pools = {int(b): bins.index[bins.eq(b)].to_numpy() for b in sorted(bins.unique())}
    random_sets: dict[str, list[str]] = {}
    for name, genes in template_sets.items():
        template_bins = bins.loc[genes].value_counts().sort_index()
        sampled: list[str] = []
        for bin_id, n_genes in template_bins.items():
            sampled.extend(rng.choice(pools[int(bin_id)], size=int(n_genes), replace=False).tolist())
        random_sets[name] = sorted(sampled)
    return random_sets


def model() -> Pipeline:
    return Pipeline(
        [
            ("scale", StandardScaler()),
            ("logistic", LogisticRegression(solver="liblinear", C=1.0, max_iter=5000, random_state=RANDOM_SEED)),
        ]
    )


def metric_dict(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    prediction = (probability >= 0.5).astype(int)
    return {
        "Accuracy": accuracy_score(y, prediction),
        "Balanced_accuracy": balanced_accuracy_score(y, prediction),
        "ROC_AUC": roc_auc_score(y, probability),
    }


def cross_validated_predictions(
    X: pd.DataFrame, y: np.ndarray, splits: list[tuple[np.ndarray, np.ndarray]]
) -> tuple[np.ndarray, pd.DataFrame]:
    probability_sum = np.zeros(len(y), dtype=float)
    prediction_counts = np.zeros(len(y), dtype=int)
    fold_rows = []
    for fold_id, (train, test) in enumerate(splits, start=1):
        fitted = model().fit(X.iloc[train], y[train])
        probability = fitted.predict_proba(X.iloc[test])[:, 1]
        probability_sum[test] += probability
        prediction_counts[test] += 1
        row = {"Fold": fold_id, **metric_dict(y[test], probability)}
        fold_rows.append(row)
    if (prediction_counts == 0).any():
        raise RuntimeError("At least one donor received no out-of-fold prediction")
    return probability_sum / prediction_counts, pd.DataFrame(fold_rows)


def stratified_bootstrap(
    y: np.ndarray, probability: np.ndarray, n_boot: int, seed: int
) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(seed)
    class_indices = [np.flatnonzero(y == label) for label in [0, 1]]
    records = []
    for _ in range(n_boot):
        index = np.concatenate([rng.choice(idx, size=len(idx), replace=True) for idx in class_indices])
        records.append(metric_dict(y[index], probability[index]))
    frame = pd.DataFrame(records)
    return {column: tuple(frame[column].quantile([0.025, 0.975])) for column in frame}


def internal_comparison(
    log_cpm: pd.DataFrame,
    design: pd.DataFrame,
    pathways: dict[str, list[str]],
    cv_repeats: int,
    random_replicates: int,
    bootstrap_replicates: int,
) -> tuple[pd.DataFrame, list[dict[str, list[str]]]]:
    y = design["Condition"].map({"ND": 0, "T2D": 1}).to_numpy(dtype=int)
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=cv_repeats, random_state=RANDOM_SEED)
    splits = list(cv.split(np.zeros(len(y)), y))
    representations = {
        "Gene features": log_cpm.T,
        "Biological pathways": score_sets(log_cpm, pathways),
    }
    prediction_table = pd.DataFrame({"Donor": design.index, "Condition": design.Condition, "True_label": y})
    summary_rows = []
    for i, (name, X) in enumerate(representations.items()):
        probability, fold_metrics = cross_validated_predictions(X, y, splits)
        fold_metrics.insert(0, "Model", name)
        fold_metrics.to_csv(METRICS / f"internal_{i}_fold_metrics.csv", index=False)
        metrics = metric_dict(y, probability)
        ci = stratified_bootstrap(y, probability, bootstrap_replicates, RANDOM_SEED + i)
        prediction_table[name] = probability
        summary_rows.append(
            {
                "Model": name,
                "Features": X.shape[1],
                **metrics,
                **{f"{key}_CI_low": value[0] for key, value in ci.items()},
                **{f"{key}_CI_high": value[1] for key, value in ci.items()},
                "CV_split_ROC_AUC_mean": fold_metrics.ROC_AUC.mean(),
                "CV_split_ROC_AUC_SD": fold_metrics.ROC_AUC.std(),
            }
        )

    rng = np.random.default_rng(RANDOM_SEED)
    random_definitions = []
    random_rows = []
    for replicate in range(1, random_replicates + 1):
        random_sets = matched_random_sets(pathways, log_cpm, rng)
        random_definitions.append(random_sets)
        X_random = score_sets(log_cpm, random_sets)
        probability, _ = cross_validated_predictions(X_random, y, splits)
        random_rows.append({"Replicate": replicate, "Features": X_random.shape[1], **metric_dict(y, probability)})
        if replicate % 10 == 0 or replicate == random_replicates:
            print(f"Internal random controls: {replicate}/{random_replicates}", flush=True)

    random_frame = pd.DataFrame(random_rows)
    random_frame.to_csv(METRICS / "internal_random_gene_set_performance.csv", index=False)
    prediction_table.to_csv(METRICS / "internal_out_of_fold_predictions.csv", index=False)

    for metric in ["Accuracy", "Balanced_accuracy", "ROC_AUC"]:
        random_values = random_frame[metric]
        summary_rows.append(
            {
                "Model": "Matched random gene sets" if metric == "Accuracy" else f"__random_{metric}",
                "Features": int(random_frame.Features.iloc[0]),
                metric: random_values.mean(),
                f"{metric}_CI_low": random_values.quantile(0.025),
                f"{metric}_CI_high": random_values.quantile(0.975),
            }
        )
    # Merge the three temporary random rows into one model row.
    summary = pd.DataFrame(summary_rows)
    random_summary = {"Model": "Matched random gene sets", "Features": len(pathways)}
    temporary = summary[summary.Model.str.startswith(("Matched random", "__random"), na=False)]
    for column in summary.columns:
        values = temporary[column].dropna()
        if column not in random_summary and len(values):
            random_summary[column] = values.iloc[0]
    summary = pd.concat(
        [summary[~summary.Model.str.startswith(("Matched random", "__random"), na=False)], pd.DataFrame([random_summary])],
        ignore_index=True,
    )
    summary.to_csv(METRICS / "internal_model_performance.csv", index=False)
    return summary, random_definitions


def technical_covariate_assessment(design: pd.DataFrame, cv_repeats: int) -> None:
    """Quantify measured donor/technical imbalance without using HbA1c as a predictor."""
    rows = []
    chemistry_table = pd.crosstab(design.Condition, design.Chemistry)
    sex_table = pd.crosstab(design.Condition, design.Sex)
    for name, table in [("Chemistry", chemistry_table), ("Sex", sex_table)]:
        result = fisher_exact(table)
        rows.append({"Variable": name, "Test": "Fisher exact", "Statistic": result.statistic, "P_value": result.pvalue})
    for variable in ["Age", "BMI", "Ncells"]:
        result = mannwhitneyu(
            design.loc[design.Condition.eq("ND"), variable],
            design.loc[design.Condition.eq("T2D"), variable],
            alternative="two-sided",
        )
        rows.append({"Variable": variable, "Test": "Mann-Whitney U", "Statistic": result.statistic, "P_value": result.pvalue})

    y = design.Condition.map({"ND": 0, "T2D": 1}).to_numpy(dtype=int)
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=cv_repeats, random_state=RANDOM_SEED)
    splits = list(cv.split(np.zeros(len(y)), y))
    chemistry = pd.get_dummies(design[["Chemistry"]], drop_first=True, dtype=float)
    covariates = pd.concat(
        [
            design[["Age", "BMI"]].astype(float),
            pd.get_dummies(design[["Chemistry", "Sex", "Race"]], drop_first=True, dtype=float),
        ],
        axis=1,
    )
    performance = []
    for name, matrix in [("Chemistry only", chemistry), ("Age/BMI/sex/race/chemistry", covariates)]:
        probability, _ = cross_validated_predictions(matrix, y, splits)
        performance.append({"Model": name, "Features": matrix.shape[1], **metric_dict(y, probability)})
    pd.DataFrame(rows).to_csv(METRICS / "measured_covariate_associations.csv", index=False)
    pd.DataFrame(performance).to_csv(METRICS / "technical_covariate_model_performance.csv", index=False)


def build_external_pseudobulk() -> tuple[pd.DataFrame, pd.DataFrame]:
    cache_counts = PROCESSED / "external_beta_pseudobulk_counts.csv"
    cache_design = PROCESSED / "external_beta_design.csv"
    if cache_counts.exists() and cache_design.exists():
        return pd.read_csv(cache_counts, index_col=0), pd.read_csv(cache_design, index_col=0)

    expression_path = EXTERNAL_RAW / "pancreas_refseq_rpkms_counts_3514sc.txt"
    sdrf = pd.read_csv(EXTERNAL_RAW / "E-MTAB-5061.sdrf.txt", sep="\t", low_memory=False)
    with expression_path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    cell_names = header[1:]
    beta_meta = sdrf.loc[sdrf["Characteristics [inferred cell type]"].eq("beta cell")].copy()
    beta_meta = beta_meta.set_index("Source Name")
    beta_names = [cell for cell in cell_names if cell in beta_meta.index]
    beta_meta = beta_meta.loc[beta_names]
    positions = [cell_names.index(cell) for cell in beta_names]
    count_positions = [2 + len(cell_names) + position for position in positions]
    donors = beta_meta["Characteristics [individual]"].astype(str)
    donor_names = list(dict.fromkeys(donors.tolist()))
    donor_to_positions = {donor: np.flatnonzero(donors.to_numpy() == donor) for donor in donor_names}

    genes: list[str] = []
    rows: list[list[float]] = []
    with expression_path.open("r", encoding="utf-8", errors="replace") as handle:
        next(handle)
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[0])
            cell_counts = np.fromiter((float(fields[j]) for j in count_positions), dtype=float, count=len(count_positions))
            rows.append([cell_counts[idx].sum() for idx in donor_to_positions.values()])
            if line_number % 10_000 == 0:
                print(f"Streamed {line_number} external genes", flush=True)
    pseudobulk = pd.DataFrame(rows, index=genes, columns=donor_names)
    pseudobulk = pseudobulk.groupby(level=0, sort=False).sum()

    design = (
        beta_meta[["Characteristics [individual]", "Characteristics [disease]"]]
        .drop_duplicates()
        .set_index("Characteristics [individual]")
        .loc[donor_names]
    )
    design["Condition"] = design["Characteristics [disease]"].map(
        {"normal": "ND", "type II diabetes mellitus": "T2D"}
    )
    design["Ncells"] = donors.value_counts().reindex(design.index)
    pseudobulk.to_csv(cache_counts)
    design.to_csv(cache_design)
    return pseudobulk, design


def normalize_external(counts: pd.DataFrame) -> pd.DataFrame:
    keep = (counts >= 10).sum(axis=1) >= 5
    filtered = counts.loc[keep]
    return np.log2(filtered.divide(filtered.sum(axis=0), axis=1) * 1_000_000 + 1)


def external_comparison(
    train_log: pd.DataFrame,
    train_design: pd.DataFrame,
    pathways: dict[str, list[str]],
    random_replicates: int,
    bootstrap_replicates: int,
) -> pd.DataFrame:
    external_counts, external_design = build_external_pseudobulk()
    external_log = normalize_external(external_counts)
    y_train = train_design.Condition.map({"ND": 0, "T2D": 1}).to_numpy(dtype=int)
    y_external = external_design.Condition.map({"ND": 0, "T2D": 1}).to_numpy(dtype=int)
    common = sorted(set(train_log.index) & set(external_log.index))
    matched_pathways = {
        name: sorted(set(genes) & set(common)) for name, genes in pathways.items()
        if len(set(genes) & set(common)) >= 10
    }
    matrices = {
        "Gene features": (train_log.loc[common].T, external_log.loc[common].T),
        "Biological pathways": (score_sets(train_log, matched_pathways), score_sets(external_log, matched_pathways)),
    }
    rows = []
    predictions = pd.DataFrame(
        {"Donor": external_design.index, "Condition": external_design.Condition, "True_label": y_external}
    )
    for i, (name, (X_train, X_external)) in enumerate(matrices.items()):
        fitted = model().fit(X_train, y_train)
        probability = fitted.predict_proba(X_external)[:, 1]
        metrics = metric_dict(y_external, probability)
        ci = stratified_bootstrap(y_external, probability, bootstrap_replicates, RANDOM_SEED + 100 + i)
        rows.append(
            {
                "Model": name,
                "Features": X_train.shape[1],
                "Shared_genes": len(common),
                **metrics,
                **{f"{key}_CI_low": value[0] for key, value in ci.items()},
                **{f"{key}_CI_high": value[1] for key, value in ci.items()},
            }
        )
        predictions[name] = probability

    # External random sets are generated from the label-free common-gene universe
    # and match the retained biological pathways in size and expression strata.
    common_train_log = train_log.loc[common]
    rng = np.random.default_rng(RANDOM_SEED + 500)
    random_rows = []
    for replicate in range(1, random_replicates + 1):
        random_sets = matched_random_sets(matched_pathways, common_train_log, rng)
        X_train = score_sets(train_log, random_sets)
        X_external = score_sets(external_log, random_sets)
        fitted = model().fit(X_train, y_train)
        probability = fitted.predict_proba(X_external)[:, 1]
        random_rows.append({"Replicate": replicate, "Features": len(random_sets), **metric_dict(y_external, probability)})
        if replicate % 10 == 0 or replicate == random_replicates:
            print(f"External random controls: {replicate}/{random_replicates}", flush=True)
    random_frame = pd.DataFrame(random_rows)
    random_frame.to_csv(METRICS / "external_random_gene_set_performance.csv", index=False)
    rows.append(
        {
            "Model": "Matched random gene sets",
            "Features": len(matched_pathways),
            "Shared_genes": len(common),
            **{metric: random_frame[metric].mean() for metric in ["Accuracy", "Balanced_accuracy", "ROC_AUC"]},
            **{f"{metric}_CI_low": random_frame[metric].quantile(0.025) for metric in ["Accuracy", "Balanced_accuracy", "ROC_AUC"]},
            **{f"{metric}_CI_high": random_frame[metric].quantile(0.975) for metric in ["Accuracy", "Balanced_accuracy", "ROC_AUC"]},
        }
    )
    summary = pd.DataFrame(rows)
    summary.to_csv(METRICS / "external_model_performance.csv", index=False)
    predictions.to_csv(METRICS / "external_comparative_predictions.csv", index=False)
    pd.DataFrame(
        [{"Pathway": name, "N_shared_genes": len(genes), "Genes": ";".join(genes)} for name, genes in matched_pathways.items()]
    ).to_csv(TABLES / "external_matched_pathway_definitions.csv", index=False)
    return summary


def pathway_enrichment(pathways: dict[str, list[str]], background: set[str]) -> pd.DataFrame:
    adjusted_path = TABLES / "beta_T2D_vs_ND_pydeseq2_adjusted.csv"
    source = adjusted_path if adjusted_path.exists() else TABLES / "beta_T2D_vs_ND_DESeq2.csv"
    de = pd.read_csv(source)
    significant = de.loc[de.padj.lt(0.05) & de.log2FoldChange.abs().ge(np.log2(1.5))].copy()
    records = []
    for direction, selected in {
        "Up in T2D": set(significant.loc[significant.log2FoldChange > 0, "gene"]) & background,
        "Down in T2D": set(significant.loc[significant.log2FoldChange < 0, "gene"]) & background,
    }.items():
        for name, genes in pathways.items():
            gene_set = set(genes) & background
            overlap = selected & gene_set
            p_value = hypergeom.sf(len(overlap) - 1, len(background), len(gene_set), len(selected))
            records.append(
                {
                    "Direction": direction,
                    "Pathway": name,
                    "Selected_genes": len(selected),
                    "Pathway_genes": len(gene_set),
                    "Overlap": len(overlap),
                    "P_value": p_value,
                    "Genes": ";".join(sorted(overlap)),
                    "DE_source": source.name,
                }
            )
    result = pd.DataFrame(records)
    result["FDR"] = np.nan
    for direction in result.Direction.unique():
        mask = result.Direction.eq(direction)
        result.loc[mask, "FDR"] = multipletests(result.loc[mask, "P_value"], method="fdr_bh")[1]
    result = result.sort_values(["Direction", "FDR"])
    result.to_csv(TABLES / "adjusted_KEGG_overrepresentation.csv", index=False)
    return result


def plot_results(
    log_cpm: pd.DataFrame,
    design: pd.DataFrame,
    internal: pd.DataFrame,
    external: pd.DataFrame | None,
    enrichment: pd.DataFrame,
) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    colors = {"ND": "#3575A8", "T2D": "#C45A3C"}

    pca_model = PCA(n_components=2)
    pca = pca_model.fit_transform(StandardScaler().fit_transform(log_cpm.T))
    pca_frame = pd.DataFrame(pca, columns=["PC1", "PC2"], index=design.index).join(design)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    sns.scatterplot(data=pca_frame, x="PC1", y="PC2", hue="Condition", style="Chemistry", palette=colors, s=75, ax=ax)
    ax.set_xlabel(f"PC1 ({pca_model.explained_variance_ratio_[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca_model.explained_variance_ratio_[1] * 100:.1f}% variance)")
    ax.set_title("Discovery beta-cell pseudobulk profiles")
    fig.tight_layout()
    fig.savefig(FIGURES / "discovery_pseudobulk_PCA.png", dpi=300)
    plt.close(fig)

    donor_plot = design.reset_index(names="Donor").sort_values(["Condition", "Ncells"], ascending=[True, False])
    fig, ax = plt.subplots(figsize=(9, 4.8))
    sns.barplot(data=donor_plot, x="Donor", y="Ncells", hue="Condition", dodge=False, palette=colors, ax=ax)
    ax.tick_params(axis="x", rotation=70, labelsize=7)
    ax.set_ylabel("Beta cells in donor pseudobulk")
    ax.set_title("Beta-cell contribution by discovery donor")
    fig.tight_layout()
    fig.savefig(FIGURES / "discovery_beta_cell_counts.png", dpi=300)
    plt.close(fig)

    random_internal = pd.read_csv(METRICS / "internal_random_gene_set_performance.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    line_colors = {"Gene features": "#2F6B9A", "Biological pathways": "#C45A3C"}
    for ax, metric, title in zip(axes, ["Balanced_accuracy", "ROC_AUC"], ["Balanced accuracy", "ROC AUC"]):
        sns.histplot(random_internal[metric], bins=12, color="#9A9A9A", ax=ax)
        for _, row in internal.loc[internal.Model.ne("Matched random gene sets")].iterrows():
            ax.axvline(row[metric], label=row.Model, linewidth=2.5, color=line_colors[row.Model])
        ax.set_xlabel(title)
        ax.set_title("Matched random-set distribution")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "internal_random_set_controls.png", dpi=300)
    plt.close(fig)

    plot_frames = [internal.assign(Cohort="Discovery (repeated CV)")]
    if external is not None:
        plot_frames.append(external.assign(Cohort="Independent E-MTAB-5061"))
    comparison = pd.concat(plot_frames, ignore_index=True)
    comparison["Model"] = comparison["Model"].replace(
        {"Gene features": "Genes", "Biological pathways": "KEGG pathways", "Matched random gene sets": "Matched random sets"}
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharey=True)
    for ax, metric, title in zip(axes, ["Balanced_accuracy", "ROC_AUC"], ["Balanced accuracy", "ROC AUC"]):
        sns.barplot(data=comparison, x="Model", y=metric, hue="Cohort", errorbar=None, ax=ax)
        ax.set_ylim(0, 1)
        ax.set_xlabel("")
        ax.set_ylabel("Score" if ax is axes[0] else "")
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=15)
    axes[0].legend_.remove()
    axes[1].legend(loc="lower right", fontsize=8)
    fig.suptitle("Representation performance across discovery and independent cohorts", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIGURES / "comparative_model_performance.png", dpi=300)
    plt.close(fig)

    if external is not None:
        random_external = pd.read_csv(METRICS / "external_random_gene_set_performance.csv")
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        for ax, metric, title in zip(axes, ["Balanced_accuracy", "ROC_AUC"], ["Balanced accuracy", "ROC AUC"]):
            sns.histplot(random_external[metric], bins=10, color="#9A9A9A", ax=ax)
            for _, row in external.loc[external.Model.ne("Matched random gene sets")].iterrows():
                ax.axvline(row[metric], label=row.Model, linewidth=2.5, color=line_colors[row.Model])
            ax.set_xlabel(title)
            ax.set_title("External matched-random distribution")
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGURES / "external_random_set_controls.png", dpi=300)
        plt.close(fig)

        external_predictions = pd.read_csv(METRICS / "external_comparative_predictions.csv")
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharey=True)
        for ax, model_name in zip(axes, ["Gene features", "Biological pathways"]):
            ordered = external_predictions.sort_values(["Condition", model_name]).reset_index(drop=True)
            palette = ordered.Condition.map(colors)
            ax.scatter(np.arange(len(ordered)), ordered[model_name], c=palette, s=55)
            ax.axhline(0.5, color="#555555", linestyle="--", linewidth=1)
            ax.set_xticks(np.arange(len(ordered)), ordered.Donor, rotation=45, ha="right")
            ax.set_title(model_name)
            ax.set_xlabel("Independent donor")
        axes[0].set_ylabel("Predicted probability of T2D")
        fig.suptitle("Fixed discovery models applied to E-MTAB-5061")
        fig.tight_layout()
        fig.savefig(FIGURES / "external_comparative_probabilities.png", dpi=300)
        plt.close(fig)

    top = enrichment.loc[enrichment.FDR.lt(0.05)].groupby("Direction", observed=True).head(10).copy()
    if len(top):
        top["minus_log10_FDR"] = -np.log10(top.FDR.clip(lower=1e-300))
        top["Label"] = top.Pathway.str.replace(r" - Homo sapiens.*$", "", regex=True)
        fig, ax = plt.subplots(figsize=(9, max(4, 0.3 * len(top))))
        sns.barplot(data=top, y="Label", x="minus_log10_FDR", hue="Direction", ax=ax)
        ax.set_xlabel("-log10(FDR)")
        ax.set_ylabel("")
        ax.set_title("Covariate-adjusted KEGG over-representation")
        fig.tight_layout()
        fig.savefig(FIGURES / "adjusted_KEGG_enrichment.png", dpi=300)
        plt.close(fig)

    adjusted_de = pd.read_csv(TABLES / "beta_T2D_vs_ND_pydeseq2_adjusted.csv")
    adjusted_de = adjusted_de.dropna(subset=["padj", "log2FoldChange"]).copy()
    adjusted_de["minus_log10_FDR"] = -np.log10(adjusted_de.padj.clip(lower=1e-300))
    adjusted_de["Status"] = "Not significant"
    effect = adjusted_de.log2FoldChange.abs().ge(np.log2(1.5))
    adjusted_de.loc[adjusted_de.padj.lt(0.05) & effect & adjusted_de.log2FoldChange.gt(0), "Status"] = "Up in T2D"
    adjusted_de.loc[adjusted_de.padj.lt(0.05) & effect & adjusted_de.log2FoldChange.lt(0), "Status"] = "Down in T2D"
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    volcano_colors = {"Not significant": "#B7B7B7", "Up in T2D": "#C45A3C", "Down in T2D": "#3575A8"}
    for status in ["Not significant", "Down in T2D", "Up in T2D"]:
        subset = adjusted_de.loc[adjusted_de.Status.eq(status)]
        ax.scatter(subset.log2FoldChange, subset.minus_log10_FDR, s=10, alpha=0.65, color=volcano_colors[status], label=status)
    ax.axhline(-np.log10(0.05), color="#555555", linestyle="--", linewidth=1)
    ax.axvline(np.log2(1.5), color="#777777", linestyle=":", linewidth=1)
    ax.axvline(-np.log2(1.5), color="#777777", linestyle=":", linewidth=1)
    for _, row in adjusted_de.sort_values("padj").head(12).iterrows():
        ax.annotate(row.gene, (row.log2FoldChange, row.minus_log10_FDR), fontsize=7, xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("log2 fold change (T2D vs ND)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Covariate-adjusted donor-level differential expression")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "adjusted_T2D_vs_ND_volcano.png", dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    ensure_dirs()
    print("Loading and validating discovery pseudobulk data", flush=True)
    _, design, log_cpm = load_discovery()
    validate_geo_metadata(design)
    print("Loading KEGG 2021 Human gene sets", flush=True)
    pathways = filter_pathways(load_kegg(), log_cpm.index)
    if len(pathways) != 291:
        raise ValueError(f"Expected 291 retained pathways, observed {len(pathways)}")
    pd.DataFrame(
        [{"Pathway": name, "N_genes": len(genes), "Genes": ";".join(genes)} for name, genes in pathways.items()]
    ).to_csv(TABLES / "discovery_pathway_definitions.csv", index=False)

    print("Running leakage-free internal representation comparison", flush=True)
    internal, _ = internal_comparison(
        log_cpm, design, pathways, args.cv_repeats, args.random_replicates, args.bootstrap_replicates
    )
    technical_covariate_assessment(design, args.cv_repeats)
    external = None
    if not args.skip_external:
        print("Processing and evaluating independent E-MTAB-5061 cohort", flush=True)
        external = external_comparison(
            log_cpm, design, pathways, args.random_replicates, args.bootstrap_replicates
        )
    if external is not None:
        random_internal = pd.read_csv(METRICS / "internal_random_gene_set_performance.csv")
        random_external = pd.read_csv(METRICS / "external_random_gene_set_performance.csv")
        diagnostics = []
        for cohort, observed, random_frame in [
            ("Discovery", internal, random_internal),
            ("Independent E-MTAB-5061", external, random_external),
        ]:
            for metric in ["Balanced_accuracy", "ROC_AUC"]:
                pathway_value = observed.loc[observed.Model.eq("Biological pathways"), metric].iat[0]
                diagnostics.append(
                    {
                        "Cohort": cohort,
                        "Metric": metric,
                        "Pathway_value": pathway_value,
                        "Random_mean": random_frame[metric].mean(),
                        "Random_2.5_percentile": random_frame[metric].quantile(0.025),
                        "Random_97.5_percentile": random_frame[metric].quantile(0.975),
                        "Pathway_percentile_among_random": 100 * (random_frame[metric] <= pathway_value).mean(),
                        "One_sided_empirical_P_pathway_superior": (1 + (random_frame[metric] >= pathway_value).sum()) / (1 + len(random_frame)),
                    }
                )
        pd.DataFrame(diagnostics).to_csv(METRICS / "pathway_vs_random_diagnostics.csv", index=False)
    enrichment = pathway_enrichment(pathways, set(log_cpm.index))
    plot_results(log_cpm, design, internal, external, enrichment)
    print("\nInternal performance")
    print(internal.to_string(index=False))
    if external is not None:
        print("\nExternal performance")
        print(external.to_string(index=False))
    print("\nComparative analysis completed successfully.")


if __name__ == "__main__":
    main()
