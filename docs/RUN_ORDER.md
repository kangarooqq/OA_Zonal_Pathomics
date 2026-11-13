
# Analysis Run Order

## OA_Zonal_Pathomics

This document describes the intended execution sequence for the three analysis components of the OA_Zonal_Pathomics repository.

**Release status:** In preparation. Scripts and documentation are being deployed and reviewed incrementally. The commands and entry points will be finalized after the corresponding scripts and input manifests have been verified.

## 1. Prepare the software environment

Clone the repository:

```bash
git clone https://github.com/kangarooqq/OA_Zonal_Pathomics.git
cd OA_Zonal_Pathomics
```

Create a Python environment and install the initial dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, the virtual-environment activation command may differ.

R analyses require a separate R environment, including DESeq2 and the packages imported by the relevant scripts. R dependencies are not installed by `requirements.txt`.

The Python dependency list is preliminary and has not yet been validated as a fully reproducible environment.

## 2. Prepare input data

The repository contains three distinct analysis components.

### Histopathology (`simtk/`)

Requires study-specific multistain histopathology images, associated specimen metadata and the corresponding analysis inputs.

These data are not assumed to be publicly redistributable. Their availability and permitted use must be established separately.

### Visium (`visium/`)

Uses GSE260926, a developing human knee spatial transcriptomic dataset, as an anatomical reference for zonal gene-expression patterns.

This dataset must not be described as an adult osteoarthritis cohort.

### Geo-seq (`geoseq/`)

Uses GSE254844, an adult human knee cartilage Geo-seq dataset.

The analysis retains 122 quality-controlled microregions from 124 original microregions:

- 20 control microregions
- 52 OA non-weight-bearing (NWB) microregions
- 50 OA weight-bearing (WB) microregions

Input files must be checked for sample identifiers, gene identifiers, dimensions and correspondence with metadata before analysis.

## 3. Recommended component order

The three components can be prepared independently. The following order is recommended for organizing the manuscript-related outputs.

### Step A — Histopathology

1. Validate image and metadata inputs.
2. Run image preprocessing and zonal feature extraction.
3. Generate histopathological measurements and quality-control outputs.
4. Run statistical analyses and predictive modelling.
5. Export manuscript-related figures and source data.

**Execution status:** Exact script commands will be documented after the `simtk/` files are deployed and reviewed.

### Step B — Visium anatomical reference

1. Load the GSE260926 spatial transcriptomic data.
2. Apply the documented preprocessing and region-of-interest definitions.
3. Calculate and visualize anatomical zonal marker patterns.
4. Generate the final Figure 5a–b reference panels.
5. Generate Supplementary Figure S4 where supported by the corresponding scripts.

**Execution status:** Exact commands and required input paths will be confirmed during `visium/` deployment.

### Step C — Adult OA Geo-seq

1. Load and validate the cleaned GSE254844 counts and metadata.
2. Confirm the retained 122 microregions and sample alignment.
3. Perform zonal differential expression analyses using DESeq2 with `sfType = "poscounts"`.
4. Calculate exploratory functional module scores from gene-wise Z-standardized log2(CPM + 1) expression.
5. Calculate molecular proteoglycan retention and collagen damage indices using their documented gene sets.
6. Generate the supported manuscript figures, supplementary outputs and source-data tables.
7. Compare numerical outputs with the manuscript before declaring reproduction successful.

**Execution status:** Script-level commands and final figure mappings remain subject to review.

## 4. Outputs

Generated figures, tables and intermediate results should be written to local output directories rather than committed as unreviewed research artifacts.

The repository's `outputs/` directory is reserved for generated files. Its contents are excluded from version control except for a placeholder file.

## 5. Reproducibility checks

Before treating a result as reproduced:

1. Verify the input dataset and sample counts.
2. Record the software environment and relevant package versions.
3. Run the corresponding script without relying on undocumented interactive-session variables.
4. Check that the expected output files are generated.
5. Compare reported sample sizes, effect estimates, statistical results and figure panels against the manuscript.
6. Document any discrepancies rather than silently modifying results to match historical values.

## 6. Known limitations

- Some analyses require data that are not included in the repository.
- R and Python environments have not yet been fully locked and independently tested.
- Some historical exploratory scripts do not correspond to final manuscript figure panels.
- Figure-to-code mappings and numerical validation remain incomplete.
- A successful script execution does not by itself establish scientific reproducibility.

This document will be updated as each analysis file is deployed and manually reviewed.
