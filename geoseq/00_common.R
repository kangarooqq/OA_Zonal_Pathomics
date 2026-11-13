# Common input contract for GSE254844 (adult OA four-zone Geo-seq).
# Metadata: sample rows, counts: gene rows, sample columns.
get_data_dir <- function() {
  d <- Sys.getenv("OA_GEOSEQ_DATA_DIR", unset = "")
  if (d == "") stop("Set OA_GEOSEQ_DATA_DIR to your authorised GSE254844 input directory")
  d
}
get_output_dir <- function() {
  d <- Sys.getenv("OA_OUTPUT_DIR", unset = "outputs/geoseq")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
DATA_DIR <- get_data_dir(); OUT_DIR <- get_output_dir()
metadata_file <- file.path(DATA_DIR, "GSE254844_Geoseq-metadata_CLEANED.csv.gz")
counts_file <- file.path(DATA_DIR, "GSE254844_processed_counts_CLEANED.txt.gz")
if (!all(file.exists(c(metadata_file,counts_file)))) stop("Required cleaned GEO inputs are missing; see README")
meta <- read.csv(gzfile(metadata_file), row.names = 1, check.names = FALSE)
counts <- read.delim(gzfile(counts_file), row.names = 1, check.names = FALSE)
if (!setequal(rownames(meta), colnames(counts))) stop("Counts/metadata samples do not match")
counts <- counts[, rownames(meta), drop = FALSE]
if (nrow(meta) != 122 || ncol(counts) != 122) stop("Expected 122 retained microregions")
if (!all(c("Condition", "Zone", "Patient_ID", "Trait") %in% names(meta))) stop("Missing cleaned metadata fields")
if (!all(colSums(counts)>0)) stop("Found zero-library sample")
cat("GSE254844: 122 aligned spatial microregions; donor count:", length(unique(meta$Patient_ID)), "\n")
