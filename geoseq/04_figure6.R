# Figure 6a–d based on supplied Step5, latest main manuscript calculation.
# Uses counts, meta and score_df4 from Figure 5d.
# Figure 6a-b | Step 1
# Focal adhesion vs candidate collagen remodeling
# ============================================================

library(dplyr)

OUT_DIR <- get_output_dir()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

collagen_genes <- c(
  "COL1A1", "COL1A2", "MMP13", "MMP1"
)

stopifnot(
  all(collagen_genes %in% rownames(counts)),
  exists("score_df4"),
  "Focal_Adhesion" %in% names(score_df4)
)

counts_aligned <- counts[
  , rownames(meta), drop = FALSE
]

lib_size <- colSums(counts_aligned)

expr <- as.matrix(
  counts_aligned[collagen_genes, , drop = FALSE]
)

log_cpm <- log2(
  sweep(expr, 2, lib_size, "/") * 1e6 + 1
)

stopifnot(all(apply(log_cpm, 1, sd) > 0))

expr_z <- t(scale(t(log_cpm)))

collagen_score <- colMeans(expr_z)

fig6_df <- data.frame(
  Sample = colnames(counts_aligned),
  Collagen_Remodeling = as.numeric(collagen_score)
) %>%
  left_join(
    score_df4 %>%
      select(
        Sample, Patient_ID,
        Condition, Zone, Focal_Adhesion
      ),
    by = "Sample"
  )

stopifnot(
  nrow(fig6_df) == 122,
  !anyNA(fig6_df)
)

oa_df <- fig6_df %>%
  filter(Condition %in% c("NWB", "WB"))

nwb_df6 <- fig6_df %>%
  filter(Condition == "NWB")

stopifnot(
  nrow(oa_df) == 102,
  nrow(nwb_df6) == 52
)

calc_cor <- function(dat, label) {
  test <- cor.test(
    dat$Focal_Adhesion,
    dat$Collagen_Remodeling,
    method = "pearson"
  )
  
  data.frame(
    Cohort = label,
    Spatial_n = nrow(dat),
    Donors = n_distinct(dat$Patient_ID),
    Pearson_r = unname(test$estimate),
    P = test$p.value
  )
}

cor_results6 <- bind_rows(
  calc_cor(oa_df, "OA_all"),
  calc_cor(nwb_df6, "NWB")
)

cat("\n=== Figure 6a-b Pearson correlations ===\n")
print(cor_results6)

adjusted_cor <- function(dat, label) {
  dat$Patient_ID <- factor(dat$Patient_ID)
  dat$Zone <- factor(dat$Zone)
  
  x <- residuals(lm(
    Focal_Adhesion ~ Patient_ID + Zone,
    data = dat
  ))
  
  y <- residuals(lm(
    Collagen_Remodeling ~ Patient_ID + Zone,
    data = dat
  ))
  
  data.frame(
    Cohort = label,
    Adjusted_r = cor(x, y)
  )
}

adjusted_results6 <- bind_rows(
  adjusted_cor(oa_df, "OA_all"),
  adjusted_cor(nwb_df6, "NWB")
)

cat("\n=== Donor/zone-adjusted residual r ===\n")
print(adjusted_results6)

write.csv(
  fig6_df,
  file.path(
    OUT_DIR,
    "Figure6ab_RealModuleScores_Audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  cor_results6,
  file.path(
    OUT_DIR,
    "Figure6ab_Correlation_Audit.csv"
  ),
  row.names = FALSE
)

cat("\n✅ Figure 6a-b Step 1 completed\n")
# 
# 




# ============================================================
# Figure 6a-b | Real module association plots
# OA n=102 / NWB n=52
# ============================================================
# =============================================================================
# Figure 6a-b | Focal adhesion vs Collagen remodeling module association
# Nature Communications Strict Publication Standard
# OA cohort (n = 102) & NWB cohort (n = 52)
# =============================================================================

library(ggplot2)
library(dplyr)
library(grid)

OUT_DIR <- get_output_dir()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DPI <- 350

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
test_oa <- cor.test(oa_df$Focal_Adhesion, oa_df$Collagen_Remodeling, method = "pearson")
r_oa <- test_oa$estimate
p_oa <- test_oa$p.value
n_oa <- nrow(oa_df)

test_nwb <- cor.test(nwb_df6$Focal_Adhesion, nwb_df6$Collagen_Remodeling, method = "pearson")
r_nwb <- test_nwb$estimate
p_nwb <- test_nwb$p.value
n_nwb <- nrow(nwb_df6)

format_p <- function(p) {
  if (p < 0.001) {
    return("italic(P) < 0.001")
  } else {
    return(sprintf("italic(P) == '%.3f'", p))
  }
}

label_oa <- sprintf("list(italic(r) == '%.3f', %s, italic(n) == %d)", r_oa, format_p(p_oa), n_oa)
label_nwb <- sprintf("list(italic(r) == '%.3f', %s, italic(n) == %d)", r_nwb, format_p(p_nwb), n_nwb)

plot6 <- bind_rows(
  oa_df %>% mutate(Cohort_Label = "OA cohort (all layers)", Panel_Tag = "a"),
  nwb_df6 %>% mutate(Cohort_Label = "Lateral condyle (NWB)", Panel_Tag = "b")
)

plot6$Cohort_Label <- factor(
  plot6$Cohort_Label,
  levels = c("OA cohort (all layers)", "Lateral condyle (NWB)")
)

stats_df <- data.frame(
  Cohort_Label = levels(plot6$Cohort_Label),
  Panel_Tag = c("a", "b"),
  Label_Expr = c(label_oa, label_nwb)
)

x_range <- range(plot6$Focal_Adhesion, na.rm = TRUE)
y_range <- range(plot6$Collagen_Remodeling, na.rm = TRUE)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
p6ab <- ggplot(plot6, aes(x = Focal_Adhesion, y = Collagen_Remodeling)) +
  geom_point(
    color = "#6A3D9A",
    fill = "#6A3D9A",
    shape = 21,
    size = 2.2,
    stroke = 0.25,
    alpha = 0.65
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#4A235A",
    fill = "#6A3D9A",
    linewidth = 0.9,
    alpha = 0.15
  ) +
  facet_wrap(~Cohort_Label, nrow = 1) +
  geom_text(
    data = stats_df,
    aes(x = -Inf, y = Inf, label = Label_Expr),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = -0.08,
    vjust = 1.45,
    size = 3.0,
    color = "#252525"
  ) +
  scale_x_continuous(
    limits = x_range,
    expand = expansion(mult = 0.08)
  ) +
  scale_y_continuous(
    limits = y_range,
    expand = expansion(mult = 0.08)
  ) +
  labs(
    x = "Focal adhesion module score",
    y = "Collagen remodeling score"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 9.5,
      face = "bold",
      color = "#252525",
      margin = margin(b = 6)
    ),
    axis.title = element_text(
      size = 8.5,
      face = "bold",
      color = "#252525"
    ),
    axis.text = element_text(
      size = 7.5,
      color = "#252525"
    ),
    axis.line = element_line(
      linewidth = 0.5,
      color = "#252525"
    ),
    axis.ticks = element_line(
      linewidth = 0.5,
      color = "#252525"
    ),
    panel.spacing = unit(1.0, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 8, l = 10)
  )

print(p6ab)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
BASE <- file.path(OUT_DIR, "Figure6ab_FocalAdhesion_CollagenAssociation_NC_Final")

ggsave(
  paste0(BASE, ".pdf"),
  p6ab,
  width = 7.08,
  height = 2.85,
  units = "in",
  device = cairo_pdf
)

ggsave(
  paste0(BASE, ".svg"),
  p6ab,
  width = 7.08,
  height = 2.85,
  units = "in"
)

ggsave(
  paste0(BASE, ".tiff"),
  p6ab,
  width = 7.08,
  height = 2.85,
  units = "in",
  dpi = DPI,
  compression = "lzw"
)

cat("\n✅ Figure 6a-b NC compliant figure successfully exported to:\n")
cat(BASE, "\n")
# 


################################################################################

# ============================================================
# Figure 6c | V5 manuscript-defined molecular indices
# PG_mol vs Col_mol | OA n=102
# PDF / SVG / TIFF | 350 dpi
# ============================================================

library(dplyr)
library(ggplot2)

OUT_DIR <- get_output_dir()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DPI <- 350

pg_synthesis <- c("ACAN", "SOX9")
pg_degradation <- c("ADAMTS4", "ADAMTS5", "CEMIP2")

col_remodeling <- c("COL1A1", "COL1A2", "MMP13", "MMP1")
col_native <- c("COL2A1", "COL9A1", "COL11A1")

gene_sets <- list(
  PG_Synthesis = pg_synthesis,
  PG_Degradation = pg_degradation,
  Collagen_Remodeling = col_remodeling,
  Native_Collagen = col_native
)

all_genes <- unique(unlist(gene_sets))

missing <- setdiff(all_genes, rownames(counts))
if (length(missing) > 0) {
  stop("Missing genes: ", paste(missing, collapse = ", "))
}

counts_aligned <- counts[, rownames(meta), drop = FALSE]
lib_size <- colSums(counts_aligned)

stopifnot(all(is.finite(lib_size)), all(lib_size > 0))

expr <- as.matrix(
  counts_aligned[all_genes, , drop = FALSE]
)

log_cpm <- log2(
  sweep(expr, 2, lib_size, "/") * 1e6 + 1
)

gene_sd <- apply(log_cpm, 1, sd)
stopifnot(all(is.finite(gene_sd)), all(gene_sd > 0))

expr_z <- t(scale(t(log_cpm)))

score <- function(genes) {
  colMeans(expr_z[genes, , drop = FALSE])
}

PG_mol <- score(pg_synthesis) - score(pg_degradation)
Col_mol <- score(col_remodeling) - score(col_native)

fig6c_df <- data.frame(
  Sample = colnames(counts_aligned),
  Patient_ID = as.character(meta$Patient_ID),
  Condition = as.character(meta$Condition),
  Zone = as.character(meta$Zone),
  PG_mol = as.numeric(PG_mol),
  Col_mol = as.numeric(Col_mol)
)

oa6c <- fig6c_df %>%
  filter(Condition %in% c("NWB", "WB"))

stopifnot(
  nrow(oa6c) == 102,
  n_distinct(oa6c$Patient_ID) == 4
)

ct <- cor.test(
  oa6c$PG_mol,
  oa6c$Col_mol,
  method = "pearson"
)

r <- unname(ct$estimate)
pval <- ct$p.value

cat("\n=== Figure 6c V5 results ===\n")
cat(sprintf("n = %d\n", nrow(oa6c)))
cat(sprintf("Pearson r = %.6f\n", r))
cat(sprintf("P = %.6g\n", pval))


oa6c$Patient_ID <- factor(oa6c$Patient_ID)
oa6c$Zone <- factor(oa6c$Zone)

rx <- residuals(lm(
  PG_mol ~ Patient_ID + Zone, data = oa6c
))
ry <- residuals(lm(
  Col_mol ~ Patient_ID + Zone, data = oa6c
))

cat(sprintf(
  "Donor/zone-adjusted residual r = %.6f\n",
  cor(rx, ry)
))

label <- sprintf(
  "r = %+.3f | P = %.3g\nn = %d spatial samples",
  r, pval, nrow(oa6c)
)

p6c <- ggplot(
  oa6c,
  aes(x = PG_mol, y = Col_mol)
) +
  geom_point(
    color = "#5E6C84",
    size = 2.2,
    alpha = 0.55
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#5E6C84",
    fill = "#5E6C84",
    alpha = 0.12,
    linewidth = 1.05
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = label,
    hjust = -0.08,
    vjust = 1.25,
    size = 4.2,
    fontface = "bold",
    color = "#252525"
  ) +
  labs(
    x = "PG retention index",
    y = "Collagen injury index"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(
      size = 14, face = "bold"
    ),
    axis.text = element_text(size = 12),
    plot.margin = margin(12, 14, 10, 12)
  )

print(p6c)

BASE <- file.path(
  OUT_DIR,
  "Figure6c_PG_Collagen_V5_Indices"
)

write.csv(
  fig6c_df,
  paste0(BASE, "_Scores.csv"),
  row.names = FALSE
)

ggsave(
  paste0(BASE, ".pdf"),
  p6c,
  width = 5.2,
  height = 4.8,
  units = "in",
  dpi = DPI
)

ggsave(
  paste0(BASE, ".svg"),
  p6c,
  width = 5.2,
  height = 4.8,
  units = "in",
  dpi = DPI
)

ggsave(
  paste0(BASE, ".tiff"),
  p6c,
  width = 5.2,
  height = 4.8,
  units = "in",
  dpi = DPI,
  compression = "lzw"
)

cat("\n✅ Figure 6c V5 calculation and export completed\n")


# ============================================================
# Figure 6d | PG-Collagen Molecular States
# Real OA Geo-seq data: 102 spatial samples, 4 donors
# PDF / SVG / TIFF | 350 dpi
# ============================================================

library(dplyr)
library(ggplot2)

OUT_DIR <- get_output_dir()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DPI <- 350

stopifnot(
  exists("oa6c"),
  all(c(
    "PG_mol", "Col_mol",
    "Condition", "Patient_ID"
  ) %in% names(oa6c)),
  nrow(oa6c) == 102
)

plot6d <- oa6c %>%
  mutate(
    Condition = factor(
      as.character(Condition),
      levels = c("NWB", "WB")
    ),
    Quadrant = case_when(
      PG_mol >= 0 & Col_mol >= 0 ~ "Q1",
      PG_mol <  0 & Col_mol >= 0 ~ "Q2",
      PG_mol <  0 & Col_mol <  0 ~ "Q3",
      PG_mol >= 0 & Col_mol <  0 ~ "Q4"
    )
  )

quadrant_summary <- plot6d %>%
  count(Quadrant, name = "n") %>%
  mutate(Percent = 100 * n / nrow(plot6d))

quadrant_summary <- data.frame(
  Quadrant = c("Q1", "Q2", "Q3", "Q4")
) %>%
  left_join(
    quadrant_summary,
    by = "Quadrant"
  ) %>%
  mutate(
    n = tidyr::replace_na(n, 0L),
    Percent = tidyr::replace_na(Percent, 0),
    Label = sprintf(
      "n = %d (%.1f%%)", n, Percent
    )
  )

cat("\n=== Figure 6d quadrant summary ===\n")
print(quadrant_summary)

condition_quadrants <- plot6d %>%
  count(Condition, Quadrant, name = "n")

cat("\n=== Quadrants by anatomical region ===\n")
print(condition_quadrants)

axis_max <- max(
  abs(c(plot6d$PG_mol, plot6d$Col_mol)),
  na.rm = TRUE
) * 1.15

quad_labels <- quadrant_summary %>%
  mutate(
    X = c(0.58, -0.58, -0.58, 0.58) * axis_max,
    Y = c(0.72, 0.72, -0.66, -0.66) * axis_max,
    Title = c(
      "High PG retention\nHigh collagen injury",
      "Low PG retention\nHigh collagen injury",
      "Low PG retention\nLow collagen injury",
      "High PG retention\nLow collagen injury"
    )
  )

cols <- c(
  NWB = "#20A39E",
  WB = "#D65A8A"
)

p6d <- ggplot(
  plot6d,
  aes(
    x = PG_mol,
    y = Col_mol,
    color = Condition
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "#777777"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "#777777"
  ) +
  
  geom_point(
    size = 2.4,
    alpha = 0.68
  ) +
  
  geom_text(
    data = quad_labels,
    aes(
      x = X,
      y = Y,
      label = Title
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold",
    lineheight = 0.95,
    color = "#444444"
  ) +
  
  geom_text(
    data = quad_labels,
    aes(
      x = X,
      y = Y - 0.12 * axis_max,
      label = Label
    ),
    inherit.aes = FALSE,
    size = 3.7,
    fontface = "bold",
    color = "#252525"
  ) +
  
  scale_color_manual(
    values = cols,
    labels = c(
      "NWB (n=52)",
      "WB (n=50)"
    )
  ) +
  
  coord_cartesian(
    xlim = c(-axis_max, axis_max),
    ylim = c(-axis_max, axis_max),
    expand = FALSE
  ) +
  
  labs(
    x = "PG retention index",
    y = "Collagen injury index",
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    text = element_text(
      family = "sans",
      color = "#252525"
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    axis.text = element_text(
      size = 12
    ),
    legend.position = "bottom",
    legend.text = element_text(
      size = 12
    ),
    plot.margin = margin(
      12, 16, 10, 12
    )
  )

print(p6d)

write.csv(
  plot6d,
  file.path(
    OUT_DIR,
    "Figure6d_Quadrant_SampleAssignment.csv"
  ),
  row.names = FALSE
)

write.csv(
  quadrant_summary,
  file.path(
    OUT_DIR,
    "Figure6d_Quadrant_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  condition_quadrants,
  file.path(
    OUT_DIR,
    "Figure6d_Quadrant_ByCondition.csv"
  ),
  row.names = FALSE
)

BASE <- file.path(
  OUT_DIR,
  "Figure6d_PG_Collagen_Quadrant"
)

ggsave(
  paste0(BASE, ".pdf"),
  p6d,
  width = 6.0,
  height = 5.6,
  units = "in",
  dpi = DPI
)

ggsave(
  paste0(BASE, ".svg"),
  p6d,
  width = 6.0,
  height = 5.6,
  units = "in",
  dpi = DPI
)

ggsave(
  paste0(BASE, ".tiff"),
  p6d,
  width = 6.0,
  height = 5.6,
  units = "in",
  dpi = DPI,
  compression = "lzw"
)

cat("\n✅ Figure 6d completed\n")
cat("Output directory:", OUT_DIR, "\n")
