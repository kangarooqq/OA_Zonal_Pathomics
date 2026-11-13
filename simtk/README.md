
# Multistain Histopathology (SIMTK)

## Overview

This directory contains computational histopathology analysis code associated with the OA_Zonal_Pathomics study of human knee osteoarthritis.

The analysis integrates quantitative features extracted from multiple histological stains to investigate cartilage zonal organization, tissue degeneration and relationships between histological measurements.

## Study scope

The histopathology component includes:

- 177 tissue specimens
- 791 histological images
- 87 within-patient paired sampling sites
- 89 patients included in predictive modelling

These figures describe different analysis units and should not be interpreted as interchangeable sample sizes.

## Analysis components

The intended workflows include:

1. Histological image preprocessing and quality control.
2. Cartilage region and zonal annotation.
3. Extraction of interpretable multistain histological features.
4. Analysis of zonal patterns and within-patient comparisons.
5. Patient-level predictive modelling and evaluation.
6. Generation of manuscript figures and supplementary outputs.

Individual scripts and their input requirements will be documented as they are added to this directory.

## Input data

The workflows may require study-specific histological images, annotations, patient metadata and derived feature tables.

These inputs are not assumed to be publicly redistributable.

Consult `../docs/INPUT_MANIFEST.csv` for the current input inventory.

Do not commit identifiable patient information or restricted source images to the public repository.

## Reproducibility

The repository is being deployed incrementally.

The presence of a script does not imply that it has been independently executed or that its numerical outputs have been reproduced.

Validation status is recorded in `../docs/VALIDATION_REPORT.md`.

## Outputs

Generated figures, tables, intermediate files and model artifacts should be stored outside the tracked source-code directories, using the output paths specified by the corresponding scripts.

## Manuscript correspondence

The final mapping between manuscript figures and analysis scripts is maintained in:

`../docs/FIGURE_TO_CODE.csv`

Mappings should be updated when scripts are added or revised.

## Status

Code deployment is complete. Scientific validation is ongoing.
