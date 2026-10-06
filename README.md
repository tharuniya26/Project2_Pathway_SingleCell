# Project 2: Pathway-Aware Single-Cell Analysis of Type 2 Diabetes

## Project goal

This project investigates whether biologically informed pathway-level representations produce more reproducible type 2 diabetes signals across independent single-cell RNA-seq cohorts than conventional gene-level representations.

## Current biological system

Human pancreatic islet beta cells.

## Main dataset

GSE221156 — human pancreatic islet single-cell RNA sequencing dataset containing non-diabetic, prediabetic, and type 2 diabetic donors.

## Main research question

Do curated biological pathway representations provide a reproducible advantage over gene-level and matched random gene-group representations when transferring donor-level disease signals across independent cohorts?

## Planned workflow

1. Audit dataset and donor metadata.
2. Obtain beta-cell donor-level pseudobulk expression.
3. Compare non-diabetic and type 2 diabetic donors.
4. Build simple gene-level baseline models.
5. Create pathway-level representations.
6. Compare real biological pathways against matched random gene sets.
7. Perform external validation using an independent cohort.
8. Evaluate both predictive performance and biological signal stability.

## Tools

- Python 3.12
- Scanpy
- AnnData
- pandas
- NumPy
- scikit-learn
- matplotlib
- VS Code
- Git/GitHub
- Kaggle for additional compute if required

## Repository note

Large raw and processed datasets are not stored in GitHub. Dataset accessions and reproducible download/preprocessing instructions will be documented instead.