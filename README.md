# OA_Zonal_Pathomics

Reproducible analysis code for interpretable zonal pathomics of human osteoarthritic cartilage, integrating multistain histopathology, Visium spatial transcriptomics and Geo-seq.

## Overview

This repository accompanies a research manuscript investigating zonal alterations in human knee osteoarthritis (OA) through computational histopathology and spatial transcriptomic analyses.

The repository is organized into three analysis components:

- **`simtk/`** — Multistain cartilage histopathology, zonal feature extraction, statistical analyses and predictive modelling.
- **`visium/`** — Spatial transcriptomic analysis of a developing human knee reference (GSE260926), used to investigate anatomical zonal marker patterns. This dataset is not an adult OA cohort.
- **`geoseq/`** — Adult human knee OA Geo-seq analysis (GSE254844), including zonal gene expression, functional module scoring and molecular matrix indices.

## Repository structure

```text
OA_Zonal_Pathomics/
├── simtk/
├── visium/
├── geoseq/
├── docs/
├── tests/
├── provenance/
├── outputs/
├── manuscript_changes/
├── README.md
└── requirements.txt
```

The directory structure above describes the planned code release. Individual files and directories will be added progressively.

## Data and reproducibility

The analyses use study-specific histopathology data and publicly available spatial transcriptomic datasets. Input data must be obtained separately where they are not distributed with this repository.

The repository will provide analysis scripts, input specifications, figure-to-code mappings and execution guidance. Reproducing individual analyses requires the corresponding input files and software dependencies.

The release is being assembled and reviewed incrementally. Availability of a script does not, by itself, indicate that its numerical outputs have been independently reproduced.

## Analysis conventions

For the Geo-seq component, 122 quality-controlled microregions are retained from 124 original microregions, comprising 102 OA and 20 control microregions.

Geo-seq differential expression analyses use DESeq2 with `sfType = "poscounts"`. Exploratory expression-based module scores use gene-wise standardized log₂(CPM + 1) values.

The currently adopted hypertrophy-associated module comprises `SPP1`, `IBSP` and `COL10A1`.

## Documentation

Additional documentation, including input manifests, execution order, figure-to-code correspondence and known reproducibility limitations, will be added under `docs/`.

## Status

**Code release in preparation.** Files are being deployed individually and manually reviewed before the repository is considered ready for public reproducibility use.


## Validation and reproducibility

The repository provides automated functional tests, analysis scripts, and documentation of numerical reproduction using the study datasets.

### Automated code validation

GitHub Actions runs Python tests on repository updates. The latest documented test run completed successfully with **11 passing synthetic functional tests**, covering Figure 3 paired analysis, Figure 4 predictive modelling, and zonal-depth computations. A Python syntax check is also included for `simtk/figure4_models.py`.

These tests verify selected code behaviours using synthetic inputs. They do not independently reproduce the manuscript's numerical findings from the original study data.

### Figure 4 real-data numerical reproduction

The Figure 4 modelling workflow has also been run using the study's specimen-level score table and final zonal pathomics feature table.

The documented analysis uses **89 distinct medial-compartment patients** and includes:

- Fivefold out-of-fold ROC-AUC analyses for three models (traditional grading, brightfield zonal pathomics, and combined predictors), using **PLM_SZ ≥ 3** as the Figure 4c endpoint.
- An additional three-model out-of-fold analysis for **PLM_Total ≥ 5**, reported as supplementary numerical results without changing the Figure 4c–d panel definitions.
- Exploratory random-forest feature importance for **PLM_Total ≥ 5**, corresponding to Figure 4d.

The validation report documents the analysis scope, generated output files, and limitations:

**[Read the validation report](docs/VALIDATION_REPORT.md)**

### Scope and limitations

The real-data numerical reproduction is distinct from the GitHub Actions synthetic tests. It does not constitute independent validation of the upstream microscopy-image segmentation, histological feature extraction, specimen-score curation, or image-to-feature processing pipeline.

Study-level input data and generated outputs are not automatically distributed through this repository; their availability and permitted use must be assessed separately.
