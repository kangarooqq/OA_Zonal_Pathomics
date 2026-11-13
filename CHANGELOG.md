
# Changelog

All notable changes to the OA_Zonal_Pathomics code repository will be documented in this file.

This repository accompanies research on zonal computational pathology and spatial transcriptomics in human knee osteoarthritis.

## [Unreleased] — 2026-11-13

### Added
- Initialized the public GitHub repository.
- Expanded `README.md` with the project overview, planned directory structure, data sources and reproducibility scope.
- Added `.gitignore` to exclude common local environments, caches, raw sequencing data and generated outputs.
- Added an initial `requirements.txt` listing Python dependencies.
- Began incremental, manually reviewed deployment of the research code and documentation.

### Analysis conventions
- The Geo-seq component uses GSE254844, with 122 retained microregions after quality control.
- Geo-seq differential expression analyses use DESeq2 with `sfType = "poscounts"`.
- The Visium component uses GSE260926 as a developing human knee anatomical reference, not an adult OA cohort.

### Pending
- Deploy and review the `simtk/`, `visium/` and `geoseq/` analysis scripts.
- Add input manifests, execution instructions and figure-to-code mappings.
- Verify numerical reproducibility using the required input datasets and software environments.
- Document unresolved differences between historical exploratory analyses and final manuscript figures.

## Release status

No versioned software release has been published yet.

The repository is under active preparation. Individual file commits should not be interpreted as completion of scientific reproducibility validation.
