# Usage: Rscript geoseq/run_all.R [analysis|all]
# Work directory must be the project root.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)>0) args[[1]] else "analysis"
source("geoseq/00_common.R")
if (mode %in% c("de","all")) source("geoseq/01_poscounts_de.R")
if (mode %in% c("analysis","all")) {
  source("geoseq/04_figure6.R")
  source("geoseq/05_supp_figure_s5abc.R")
  source("geoseq/06_figure6_source_data.R")
}
