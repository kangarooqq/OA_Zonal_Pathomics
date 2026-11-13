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
