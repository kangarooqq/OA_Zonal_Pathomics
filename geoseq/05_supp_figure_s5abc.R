# Supplementary Figures S5a–S5d | Quality control and analysis plots
# Supplementary Figure S5a | Library-size quality control
# GSE254844 | 124 spatial samples
#
# Plotting and export settings:
# 1. Enlarged typography (Titles 13pt, Axes 12pt, Ticks 10pt) for high readability
# 2. Standard bold 'a' panel identifier
# 3. Mathematical expressions for scientific Y-axis log10 labels
# 4. Output to unified directory in PDF / SVG / TIFF (350 dpi, LZW)
# ==============================================================

library(ggplot2)

DATA_DIR <- get_data_dir()
OUT_DIR <- get_output_dir()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DPI <- 350

counts <- read.table(
  file.path(DATA_DIR, "GSE254844_processed_counts.txt.gz"),
  header = TRUE,
  row.names = 1,
  sep = "",
  check.names = FALSE
)

qc <- data.frame(
  Sample = colnames(counts),
  LibrarySize = colSums(counts),
  DetectedGenes = colSums(counts > 0)
)
stopifnot(nrow(qc) == 124)

qc$Status <- ifelse(
  qc$Sample %in% c("Control.AS6", "NWB.DZ9"),
  "Excluded (n = 2)",
  "Retained (n = 122)"
)

qc <- qc[order(qc$LibrarySize, decreasing = TRUE), ]
qc$Rank <- seq_len(nrow(qc))

write.csv(
  qc,
  file.path(OUT_DIR, "S5a_GSE254844_LibrarySize_QC_Audit.csv"),
  row.names = FALSE
)

p5a <- ggplot(qc, aes(x = Rank, y = LibrarySize)) +
  geom_col(
    aes(fill = Status),
    width = 0.85,
    color = NA
  ) +
  scale_fill_manual(
    values = c(
      "Retained (n = 122)" = "#20A39E",  # Academic Teal
      "Excluded (n = 2)"   = "#D65A8A"   # Vibrant Rose
    )
  ) +
  scale_y_log10(
    breaks = c(1e2, 1e3, 1e4, 1e5, 1e6),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_continuous(
    breaks = seq(0, 125, 25),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Library-size distribution",
    subtitle = "GSE254844 | 124 spatial samples",
    x = "Samples ranked by library size",
    y = expression(bold("Total counts (")*bolditalic(log)[10]*bold(" scale)")),
    fill = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(
      face = "bold", colour = "#252525", size = 13.5, margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      colour = "#555555", size = 10.5, margin = margin(b = 10)
    ),
    axis.title = element_text(
      face = "bold", colour = "#252525", size = 12.0
    ),
    axis.text = element_text(
      colour = "#252525", size = 10.5, face = "bold"
    ),
    axis.line = element_line(
      linewidth = 0.65, color = "#252525"
    ),
    axis.ticks = element_line(
      linewidth = 0.65, color = "#252525"
    ),
    legend.position = "top",
    legend.justification = "right",
    legend.text = element_text(
      size = 10.5, face = "bold", color = "#252525"
    ),
    legend.key.size = unit(0.45, "cm"),
    legend.margin = margin(b = -15),
    plot.margin = margin(t = 12, r = 16, b = 10, l = 12)
  )

print(p5a)

BASE_PATH <- file.path(OUT_DIR, "Supplementary_Figure_S5a_LibrarySize")

ggsave(
  paste0(BASE_PATH, ".pdf"),
  p5a, width = 7.08, height = 3.60, units = "in",
  device = cairo_pdf
)

ggsave(
  paste0(BASE_PATH, ".svg"),
  p5a, width = 7.08, height = 3.60, units = "in"
)

ggsave(
  paste0(BASE_PATH, ".tiff"),
  p5a, width = 6, height = 4, units = "in",
  dpi = DPI, compression = "lzw"
)

cat(paste0("\n", paste(rep("=", 65), collapse = ""), "\n"))
cat("✔ Supplementary Figure S5a exported successfully:\n")
cat(sprintf("  - PDF : %s.pdf\n", BASE_PATH))
cat(sprintf("  - SVG : %s.svg\n", BASE_PATH))
cat(sprintf("  - TIFF: %s.tiff (350 dpi, LZW)\n", BASE_PATH))
cat(paste0(paste(rep("=", 65), collapse = ""), "\n"))

###############################################################################
# ============================================================
# Supplementary Figure S5b | Excluded sample library complexity
# GSE254844 | Control.AS6 and NWB.DZ9
#
# Plotting and export settings:
# 1. Zero data occlusion: moved sample labels to safe open quadrant
# 2. Enlarged typography (Titles 13pt, Axes 12.5pt, Ticks 10.5pt, bold)
# 3. Standard bold 'b' panel identifier
# 4. Strict scientific LaTeX log10 scale labels
# 5. Output to unified directory: PDF / SVG / TIFF (350 dpi, LZW)
# ============================================================

library(ggplot2)

OUT_DIR <- get_output_dir()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DPI <- 350

qc <- read.csv(
  file.path(OUT_DIR, "S5a_GSE254844_LibrarySize_QC_Audit.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

qc$Status <- ifelse(
  grepl("^Retained", trimws(qc$Status)),
  "Retained",
  ifelse(grepl("^Excluded", trimws(qc$Status)), "Excluded", NA_character_)
)
stopifnot(nrow(qc) == 124, sum(qc$Status == "Retained") == 122, sum(qc$Status == "Excluded") == 2)

targets <- c("Control.AS6", "NWB.DZ9")
retained <- qc[qc$Status == "Retained", ]

reference <- do.call(
  rbind,
  lapply(targets, function(s) {
    data.frame(
      Sample = retained$Sample,
      LibrarySize = retained$LibrarySize,
      DetectedGenes = retained$DetectedGenes,
      Panel = s
    )
  })
)

highlight <- qc[match(targets, qc$Sample), c("Sample", "LibrarySize", "DetectedGenes")]
highlight$Panel <- factor(highlight$Sample, levels = targets)
reference$Panel <- factor(reference$Panel, levels = targets)

label_df <- data.frame(
  Panel = factor(targets, levels = targets),
  x = rep(max(qc$LibrarySize) * 0.95, 2),
  y = rep(35, 2),
  label = c(
    "Control.AS6\n181 counts · 102 genes",
    "NWB.DZ9\n70 counts · 35 genes"
  ),
  color = c("#D65A8A", "#6A3D9A")
)

TEAL   <- "#20A39E"
ROSE   <- "#D65A8A"
PURPLE <- "#6A3D9A"
DARK   <- "#252525"

med_lib  <- median(retained$LibrarySize)
med_gene <- median(retained$DetectedGenes)

p5b <- ggplot() +
  geom_point(
    data = reference,
    aes(x = LibrarySize, y = DetectedGenes),
    size = 2.0,
    color = "#94A3B8",
    alpha = 0.55
  ) +
  geom_vline(
    xintercept = med_lib,
    color = TEAL,
    linetype = "dashed",
    linewidth = 0.65,
    alpha = 0.75
  ) +
  geom_hline(
    yintercept = med_gene,
    color = TEAL,
    linetype = "dashed",
    linewidth = 0.65,
    alpha = 0.75
  ) +
  geom_point(
    data = highlight,
    aes(x = LibrarySize, y = DetectedGenes, color = Panel),
    size = 7.5,
    shape = 21,
    fill = "#FFFFFF",
    stroke = 2.2
  ) +
  geom_point(
    data = highlight,
    aes(x = LibrarySize, y = DetectedGenes, color = Panel),
    size = 11.5,
    shape = 1,
    stroke = 1.0,
    alpha = 0.8
  ) +
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = label, color = Panel),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 0,
    size = 4.2,
    fontface = "bold",
    lineheight = 1.05,
    fill = alpha("#FFFFFF", 0.90),
    label.size = 0,
    label.padding = unit(0.35, "lines"),
    show.legend = FALSE
  ) +
  facet_wrap(~Panel, nrow = 1) +
  scale_color_manual(
    values = c("Control.AS6" = ROSE, "NWB.DZ9" = PURPLE)
  ) +
  scale_x_log10(
    breaks = c(1e2, 1e3, 1e4, 1e5, 1e6),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  coord_cartesian(
    xlim = c(35, max(qc$LibrarySize) * 1.5),
    ylim = c(18, max(qc$DetectedGenes) * 1.8)
  ) +
  labs(
    title = "Excluded samples: library complexity",
    subtitle = "124 spatial samples | 122 retained references (dashed: reference medians)",
    x = expression(bold("Library size (total counts, ")*bolditalic(log)[10]*bold(")")),
    y = expression(bold("Detected genes (")*bolditalic(log)[10]*bold(")"))
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 13.5, face = "bold", color = DARK, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.0, color = "#555555", margin = margin(b = 8)),
    strip.background = element_blank(),
    strip.text = element_text(size = 12.5, face = "bold", color = DARK, margin = margin(b = 6)),
    axis.title = element_text(size = 11.5, face = "bold", color = DARK),
    axis.text = element_text(size = 10.0, face = "bold", color = DARK),
    axis.line = element_line(linewidth = 0.65, color = DARK),
    axis.ticks = element_line(linewidth = 0.65, color = DARK),
    legend.position = "none",
    panel.spacing = grid::unit(0.8, "cm"),
    plot.margin = margin(t = 12, r = 16, b = 10, l = 12)
  )

print(p5b)

BASE_PATH <- file.path(OUT_DIR, "Supplementary_Figure_S5b_ExcludedSamples")

ggsave(
  paste0(BASE_PATH, ".pdf"),
  p5b, width = 7.08, height = 3.65, units = "in",
  device = cairo_pdf
)

ggsave(
  paste0(BASE_PATH, ".svg"),
  p5b, width = 7.08, height = 3.65, units = "in"
)

ggsave(
  paste0(BASE_PATH, ".tiff"),
  p5b, width = 7.08, height = 3.65, units = "in",
  dpi = DPI, compression = "lzw"
)

cat(paste0("\n", paste(rep("=", 65), collapse = ""), "\n"))
cat("✔ Supplementary Figure S5b exported successfully:\n")
cat(sprintf("  - PDF : %s.pdf\n", BASE_PATH))
cat(sprintf("  - SVG : %s.svg\n", BASE_PATH))
cat(sprintf("  - TIFF: %s.tiff (350 dpi, LZW)\n", BASE_PATH))
cat(paste0(paste(rep("=", 65), collapse = ""), "\n"))

################################################################################
# =============================================================
# Supplementary Figure S5c | Donor sex-marker expression audit
# GSE254844 | XIST + 5 Y-linked genes across 5 donors
#
# Plotting and export settings:
# 1. Zero text overlap: clean donor labels, removed congested multiline titles
# 2. Enlarged typography (Gene labels 11pt bold, Ticks 11pt, clean colorbar)
# 3. Explicit biological sex status annotation
# 4. Standard bold 'c' panel identifier
# 5. Output: PDF / SVG / TIFF (350 dpi, LZW)
# =============================================================

library(ggplot2)

DATA_DIR <- get_data_dir()
OUT_DIR <- get_output_dir()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DPI <- 350

# -------------------------------------------------------------
meta0 <- read.csv(
  file.path(DATA_DIR, "GSE254844_Geoseq-metadata.csv.gz"),
  row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
)

counts0 <- read.table(
  file.path(DATA_DIR, "GSE254844_processed_counts.txt.gz"),
  header = TRUE, row.names = 1, sep = "", check.names = FALSE
)

rownames(meta0) <- colnames(counts0)

excluded <- c("Control.AS6", "NWB.DZ9")
keep <- !(colnames(counts0) %in% excluded)
counts <- counts0[, keep, drop = FALSE]
meta <- meta0[keep, , drop = FALSE]

meta$Donor <- sub("-[0-9]+$", "", as.character(meta$SampleID))
meta$Sex_GEO <- tools::toTitleCase(tolower(trimws(as.character(meta$sex))))

# -------------------------------------------------------------
genes <- c("XIST", "RPS4Y1", "DDX3Y", "USP9Y", "KDM5D", "UTY")
lib <- colSums(counts)
expr <- as.matrix(counts[genes, , drop = FALSE])
cpm <- sweep(expr, 2, lib, "/", check.margin = FALSE) * 1e6
log2cpm <- log2(cpm + 1)

donors <- c("SC1", "SOA1", "SOA2", "SOA3", "SOA4")

donor_expr <- sapply(donors, function(d) {
  rowMeans(log2cpm[, meta$Donor == d, drop = FALSE])
})

# -------------------------------------------------------------
heat <- expand.grid(Gene = genes, Donor = donors, stringsAsFactors = FALSE)
heat$Expression <- mapply(function(g, d) donor_expr[g, d], heat$Gene, heat$Donor)

heat$Gene <- factor(heat$Gene, levels = rev(genes))
heat$Donor <- factor(heat$Donor, levels = donors)

donor_labels <- c(
  "SC1\n(Female)",
  "SOA1\n(Female)",
  "SOA2\n(Female)",
  "SOA3\n(Male)",
  "SOA4\n(Male)"
)
names(donor_labels) <- donors

# -------------------------------------------------------------
## -------------------------------------------------------------
# -------------------------------------------------------------
p5c <- ggplot(heat, aes(x = Donor, y = Gene, fill = Expression)) +
  geom_tile(color = "#FFFFFF", linewidth = 1.2) +
  scale_fill_gradientn(
    colours = c("#F1F5F9", "#CBD5E1", "#94D2BD", "#20A39E", "#0F4C5C"),
    name = expression(bold("Mean\n")*bolditalic(log)[2]*bold("(CPM+1)")),
    limits = c(0, 6.2),
    breaks = seq(0, 6, 2)
  ) +
  scale_x_discrete(
    labels = donor_labels,
    position = "bottom"
  ) +
  coord_fixed(ratio = 0.85) +
  labs(
    title = "Donor sex-marker expression audit",
    subtitle = "GSE254844 | 122 spatial samples | Corrected 5 donor profiles",
    x = NULL, y = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 13.5, face = "bold", color = "#252525", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.0, color = "#555555", margin = margin(b = 14)),
    axis.text.x = element_text(size = 11.0, face = "bold", color = "#252525", lineheight = 1.15),

    axis.text.y = element_text(size = 11.0, face = "bold.italic", color = "#252525"),

    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9.5, face = "bold", color = "#252525", lineheight = 1.1),
    legend.text = element_text(size = 9.0, color = "#252525", face = "bold"),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(1.1, "cm"),
    plot.margin = margin(t = 12, r = 16, b = 10, l = 12)
  )

print(p5c)

# -------------------------------------------------------------
BASE_PATH <- file.path(OUT_DIR, "Supplementary_Figure_S5c_DonorSexAudit")

ggsave(
  paste0(BASE_PATH, ".pdf"),
  p5c, width = 6.8, height = 4.6, units = "in",
  device = cairo_pdf
)

ggsave(
  paste0(BASE_PATH, ".svg"),
  p5c, width = 6.8, height = 4.6, units = "in"
)

ggsave(
  paste0(BASE_PATH, ".tiff"),
  p5c, width = 6.8, height = 4.6, units = "in",
  dpi = DPI, compression = "lzw"
)

cat(paste0("\n", paste(rep("=", 65), collapse = ""), "\n"))
cat("✔ Supplementary Figure S5c exported successfully:\n")
cat(sprintf("  - PDF : %s.pdf\n", BASE_PATH))
cat(sprintf("  - SVG : %s.svg\n", BASE_PATH))
cat(sprintf("  - TIFF: %s.tiff (350 dpi, LZW)\n", BASE_PATH))
cat(paste0(paste(rep("=", 65), collapse = ""), "\n"))
#########################################################
#
# ============================================================
# Supplementary Figure S5d | Zonal DEG distribution
# GSE254844 | DESeq2 poscounts across AS / SZ / MZ / DZ
#
# Plotting and export settings:
# 1. Pure rigorous DESeq2 (sfType = "poscounts") zonal DEG profile (no counts+1)
# 2. External bold count labels to eliminate narrow-bar text collisions in DZ
# 3. Unified color palette, Arial typography, and Source Data CSV export
# 4. Output: PDF / SVG / TIFF (350 dpi, LZW) using consistent S5 panel dimensions
# ============================================================

library(ggplot2)
library(dplyr)

DATA_DIR <- "/home/kanga/OANC_project_20260403/Spatial Multi-omics 20261024"
OUT_DIR  <- "/home/kanga/OANC_project_20260403/CP_FIG_20261107"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DPI <- 350

# 开关：TRUE = 同时展示上调与下调 (对应 Methods 的 |log2FC| >= 1.0)；
#       FALSE = 仅展示上调基因 (2,061 / 1,165 / 500 / 91) TRUE
SHOW_DOWNREGULATED <- FALSE

ZONES <- c("AS", "SZ", "MZ", "DZ")
ROSE  <- "#D65A8A"  # Upregulated in OA
TEAL  <- "#20A39E"  # Downregulated in OA
DARK  <- "#252525"

# ------------------------------------------------------------
# 1. 仅读取标准 DESeq2 poscounts 四层结果表
# ------------------------------------------------------------
read_deg <- function(path) {
  stopifnot(file.exists(path))
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  df$Gene <- trimws(as.character(df$Gene))
  df$log2FoldChange <- suppressWarnings(as.numeric(df$log2FoldChange))
  df$padj <- suppressWarnings(as.numeric(df$padj))
  df
}

count_sig <- function(df, direction = "Up") {
  keep <- is.finite(df$padj) & df$padj < 0.05 & is.finite(df$log2FoldChange)
  if (direction == "Up") {
    keep <- keep & df$log2FoldChange >= 1.0
  } else {
    keep <- keep & df$log2FoldChange <= -1.0
  }
  length(unique(df$Gene[keep]))
}

deg_summary <- do.call(rbind, lapply(ZONES, function(z) {
  f_pos <- file.path(DATA_DIR, paste0("4A_DEG_OA_vs_Ctrl_", z, "_Rigorous_poscounts.csv"))
  df <- read_deg(f_pos)
  data.frame(
    Zone = z,
    Upregulated   = count_sig(df, "Up"),
    Downregulated = count_sig(df, "Down"),
    stringsAsFactors = FALSE
  )
}))

# 导出 S5d 对应的统计汇总数据
write.csv(
  deg_summary,
  file.path(OUT_DIR, "S5d_GSE254844_Poscounts_Zonal_DEGs.csv"),
  row.names = FALSE
)
print(deg_summary, row.names = FALSE)

# ------------------------------------------------------------
# 2. 构建绘图数据表
# ------------------------------------------------------------
if (SHOW_DOWNREGULATED) {
  plot_df <- bind_rows(
    data.frame(Zone = deg_summary$Zone, Direction = "Upregulated",   Genes = deg_summary$Upregulated),
    data.frame(Zone = deg_summary$Zone, Direction = "Downregulated", Genes = deg_summary$Downregulated)
  )
  plot_df$Direction <- factor(plot_df$Direction, levels = c("Downregulated", "Upregulated"))
  sub_expr <- expression(
    bold("OA vs Control | DESeq2 poscounts | |")*bolditalic(log)[2]*bold("FC| ≥ 1.0, BH-adjusted ")*bolditalic(P)*bold(" < 0.05")
  )
} else {
  plot_df <- data.frame(
    Zone = deg_summary$Zone,
    Direction = "Upregulated",
    Genes = deg_summary$Upregulated
  )
  sub_expr <- expression(
    bold("OA vs Control | DESeq2 poscounts | ")*bolditalic(log)[2]*bold("FC ≥ 1.0, BH-adjusted ")*bolditalic(P)*bold(" < 0.05")
  )
}

plot_df$Zone <- factor(plot_df$Zone, levels = rev(ZONES))
plot_df$Label <- format(plot_df$Genes, big.mark = ",", trim = TRUE)

# ------------------------------------------------------------
# 3. 绘制条形图
# ------------------------------------------------------------
if (SHOW_DOWNREGULATED) {
  dodge_pos <- position_dodge(width = 0.72)
  p5d <- ggplot(plot_df, aes(x = Genes, y = Zone, fill = Direction)) +
    geom_col(position = dodge_pos, width = 0.64, color = "#FFFFFF", linewidth = 0.35) +
    geom_text(
      aes(x = Genes + 35, label = Label),
      hjust = -0.18,
      size = 3.5,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("Upregulated" = ROSE, "Downregulated" = TEAL),
      breaks = c("Upregulated", "Downregulated"),
      labels = c("Upregulated in OA", "Downregulated in OA")
    ) +
    scale_color_manual(
      values = c("Upregulated" = ROSE, "Downregulated" = TEAL)
    )
} else {
  p5d <- ggplot(plot_df, aes(x = Genes, y = Zone)) +
    geom_col(fill = ROSE, width = 0.58, color = NA) +
    geom_text(
      aes(label = Label),
      hjust = -0.18,
      size = 3.8,
      color = DARK,
      fontface = "bold"
    )
}

p5d <- p5d +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Zonal differentially expressed genes in OA",
    subtitle = sub_expr,
    x = "Number of significant differentially expressed genes",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 13.5, face = "bold", color = DARK, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.0, color = "#555555", margin = margin(b = 12)),
    axis.text.y = element_text(size = 11.5, face = "bold", color = DARK),
    axis.text.x = element_text(size = 10.0, face = "bold", color = DARK),
    axis.title.x = element_text(size = 11.0, face = "bold", color = DARK, margin = margin(t = 8)),
    axis.line = element_line(linewidth = 0.65, color = DARK),
    axis.ticks = element_line(linewidth = 0.65, color = DARK),
    legend.position = if (SHOW_DOWNREGULATED) "bottom" else "none",
    legend.text = element_text(size = 10.0, face = "bold", color = DARK),
    legend.key.size = unit(0.45, "cm"),
    legend.margin = margin(t = 6),
    plot.margin = margin(t = 12, r = 20, b = 10, l = 12)
  )

print(p5d)

# ------------------------------------------------------------
# 4. 导出三格式文件至统一目录 (350 DPI)
# ------------------------------------------------------------
BASE_PATH <- file.path(OUT_DIR, "Supplementary_Figure_S5d_ZonalDEGs")

ggsave(
  paste0(BASE_PATH, ".pdf"),
  p5d, width = 7.08, height = 3.80, units = "in",
  device = cairo_pdf
)

ggsave(
  paste0(BASE_PATH, ".svg"),
  p5d, width = 7.08, height = 3.80, units = "in"
)

ggsave(
  paste0(BASE_PATH, ".tiff"),
  p5d, width = 7.08, height = 3.80, units = "in",
  dpi = DPI, compression = "lzw"
)

cat(paste0("\n", paste(rep("=", 65), collapse = ""), "\n"))
cat("✔ Supplementary Figure S5d (Zonal DEGs, poscounts) exported successfully:\n")
cat(sprintf("  - PDF : %s.pdf\n", BASE_PATH))
cat(sprintf("  - SVG : %s.svg\n", BASE_PATH))
cat(sprintf("  - TIFF: %s.tiff (350 dpi, LZW)\n", BASE_PATH))
cat(paste0(paste(rep("=", 65), collapse = ""), "\n"))


