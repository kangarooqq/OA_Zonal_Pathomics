# Validation Report

## OA_Zonal_Pathomics

**Report date:** 2026-11-13  

## 1. Purpose

This report distinguishes code-level checks from numerical reproduction of the scientific results.

A successful syntax check, synthetic test or file deployment does not establish that manuscript figures and statistical results have been reproduced.

## 2. Checks reported during release preparation

### Python checks

During preparation of a local release candidate, four synthetic Python tests were reported as passing.

These checks provide limited evidence about selected code paths. They do not verify the scientific results or establish that the current GitHub repository passes the same tests.

**GitHub status:** The corresponding test files have not yet been deployed or rerun from the published repository.

### R checks

Some historical scripts depend on objects already present in an interactive R session, including `counts`, `meta` and `score_df4`.

## 3. Input dataset checks

### Histopathology

Study-specific multistain images, metadata and annotations are required.

**Status:** Full input-level and end-to-end validation pending.

## 4. Current analysis conventions

### Geo-seq differential expression

DESeq2 normalization uses:

`sfType = "poscounts"`

Experimental pseudocount-based modifications to raw DESeq2 count inputs are not part of the adopted analysis.

### Geo-seq module scoring

Exploratory module scores are calculated from gene-wise standardized log2(CPM + 1) expression.

Any numerical results previously associated with a different gene set require verification before reuse.

### Figure correspondence

Final Figure 5 contains panels a and b only.

Final Figure 6e is a conceptual schematic, not a gene-expression plot.

Historical exploratory figure labels must not be treated as final manuscript figure identifiers without review.
### Visium

Dataset: GSE260926.

Role: Developing human knee anatomical reference, not an adult OA cohort.

**Status:** Final Figure 5a–b and supplementary output reproduction pending independent verification.

### Geo-seq

Dataset: GSE254844.

Expected sample accounting:

- Original microregions: 124
- Excluded during quality control: 2
- Retained microregions: 122
- Control microregions: 20
- OA non-weight-bearing microregions: 52
- OA weight-bearing microregions: 50

**Status:** Sample accounting documented; fresh execution against the required inputs remains pending.

## 5. Requirements for declaring numerical reproduction

A manuscript result should be marked as numerically reproduced only when:

1. The input data and sample selection are documented.
2. The corresponding code executes in a recorded environment.
3. The expected outputs are generated.
4. Numerical estimates and sample sizes are compared with the manuscript.
5. Any discrepancies and accepted tolerances are recorded.
6. The validation evidence is retained.

## 6. Conclusion

The repository is undergoing incremental deployment and manual review.

Initial code organization and limited synthetic testing have been performed during release preparation. However, full execution with original research inputs, numerical agreement with manuscript results and independent end-to-end reproduction have not yet been established.

This report will be updated as additional scripts, tests and validation evidence become available.


## GitHub Actions validation update — 2026-11-13

**Workflow:** Python tests (`.github/workflows/python-tests.yml`)

**Run:** [Python tests #5](https://github.com/kangarooqq/OA_Zonal_Pathomics/actions/runs/43444401774)

**Commit:** `a56b84f`

**Status:** SUCCESS — 7 tests passed.

The GitHub Actions workflow ran under Python 3.11 on Linux and completed all seven synthetic functional tests successfully.

### Figure 3 paired analysis (4 tests)

- Complete 87-pair synthetic input produces the expected output structure.
- An incomplete 86-pair input is rejected for every outcome.
- Duplicate patient-region records are rejected.
- Paired Wilcoxon results agree with an independent calculation.

### Zonal depth computation (3 tests)

- Output structure and background handling are verified.
- Three-zone assignments and depth thresholds are verified.
- Micrometres-per-pixel scaling affects thickness without changing normalized depth.

### Scope and limitations

These tests verify selected software functionality using synthetic inputs. They do not establish numerical reproduction of Figure 3 from the actual 87 patient pairs, validate cartilage segmentation on real microscopy images, or reproduce other manuscript figures and statistical results.

**Scientific validation status: pending real-data reproduction.**

## Figure 4: Real-data numerical reproduction (2026-11-13)

### Input datasets and cohort verification

The Figure 4 analyses were rerun using the following two input datasets:

- `OA_Specimen_Scores_177.csv` — corrected specimen-level pathology scores and PLM outcomes.
- `OA_Zonal_Pathomics_Master_177_20261101.csv` — final specimen-level zonal pathomic features.

Each input contained 177 specimen records. The datasets were merged using `Patient`, `Region`, and `SpecimenID`. No duplicate specimen keys or conflicting overlapping fields were identified.

After restricting the merged dataset to medial specimens with non-missing `PLM_Total`, the analysis cohort comprised 89 records from 89 distinct patients. All five traditional grading predictors, 27 brightfield zonal pathomic predictors, and both PLM outcomes were available without missing values.

### Figure 4c: Superficial collagen injury

Endpoint: `PLM_SZ >= 3`.

The endpoint was positive in 48 of 89 patients (53.9%).

The original fivefold stratified cross-validation procedure was used, with fold-specific standardisation, logistic regression and random forest probability averaging, and 500 bootstrap resamples of out-of-fold predictions for confidence intervals.

| Model | OOF ROC-AUC | 95% CI |
|---|---:|---:|
| Traditional grading | 0.780 | 0.675–0.867 |
| Brightfield zonal pathomics | 0.727 | 0.609–0.830 |
| Combined | 0.803 | 0.698–0.884 |

These results reproduce the manuscript's reported Figure 4c AUCs and confidence intervals to three decimal places.

### Additional full-thickness collagen injury analysis

Endpoint: `PLM_Total >= 5`.

The endpoint was positive in 47 of 89 patients (52.8%). The same three model definitions and fivefold out-of-fold evaluation procedure were applied.

| Model | OOF ROC-AUC | 95% CI |
|---|---:|---:|
| Traditional grading | 0.799 | 0.698–0.881 |
| Brightfield zonal pathomics | 0.662 | 0.546–0.775 |
| Combined | 0.775 | 0.674–0.869 |

The traditional grading and combined model AUCs reproduce the corresponding values reported in the manuscript. This supplementary numerical analysis does not change the definition or visual content of Figure 4c or Figure 4d.

### Figure 4d: Random-forest feature importance

For `PLM_Total >= 5`, a random forest was fitted using all 89 patients and the 27 brightfield zonal pathomic features. The ten highest-ranked features were:

| Rank | Feature | Mean decrease in impurity |
|---|---|---:|
| 1 | `MZ_HE_Basophilia` | 0.070 |
| 2 | `SZ_SurfaceAngle_deg` | 0.068 |
| 3 | `DZ_SurfaceAngle_deg` | 0.060 |
| 4 | `SZ_HE_Eosinophilia` | 0.058 |
| 5 | `Grad_Angle_DZ_SZ` | 0.054 |
| 6 | `DZ_CloningRate` | 0.045 |
| 7 | `Surface_Roughness` | 0.044 |
| 8 | `SZ_CloningRate` | 0.043 |
| 9 | `DZ_ChondronArea_um2` | 0.041 |
| 10 | `DZ_SafO_PG` | 0.040 |

The ranking and displayed three-decimal importance values agree with the supplied final Figure 4d.

### Reproduction outputs

The analysis generated the following numerical source-data files:

- `Figure4c_OOF_AUC_Summary.csv`
- `Figure4c_PatientLevel_OOF_Predictions.csv`
- `Figure4_FullThickness_PLMTotal5_OOF_AUC_Summary.csv`
- `Figure4_FullThickness_PLMTotal5_PatientLevel_OOF_Predictions.csv`
- `Figure4d_All27_FeatureImportance.csv`

Figure 4c ROC and Figure 4d feature-importance visualisations were also generated in PDF, SVG, and TIFF formats.

The numerical outputs were generated in the real-data reproduction session; they are not automatically included in the public GitHub repository. The input datasets should not be assumed to be publicly distributable.

### Scope and limitations

This verification establishes numerical reproducibility of the Figure 4 modelling analyses from the supplied specimen-level CSV datasets under the documented Python implementation.

It does not independently validate the upstream image segmentation, histological feature extraction, specimen-level score curation, or raw microscopy-to-feature processing pipeline.

The Figure 4d impurity-based feature importance is exploratory and derived from a model fitted to the full analysis cohort. It should not be interpreted as cross-validated predictive importance or evidence of causal relationships.

These real-data numerical reproduction results are distinct from the GitHub Actions synthetic functional tests and Python syntax checks.
