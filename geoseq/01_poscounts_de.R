# Final DESeq2 primary method. No pseudocount addition.
# Based on Geo-seq Stage1 final poscounts arm and Stage2 paired contrast.
suppressPackageStartupMessages(library(DESeq2))
zones <- c("AS", "SZ", "MZ", "DZ")
result_summary <- list()
for (z in zones) {
  idx <- as.character(meta$Zone) == z
  sub_m <- meta[idx, , drop=FALSE]
  sub_c <- as.matrix(counts[, idx, drop=FALSE])
  keep <- rowSums(sub_c)>0
  sub_m$Trait <- factor(as.character(sub_m$Trait), levels=c("Control","OA"))
  dds <- DESeqDataSetFromMatrix(countData=sub_c[keep,,drop=FALSE], colData=sub_m, design=~Trait)
  dds <- DESeq(dds,sfType="poscounts",quiet=TRUE)
  res <- as.data.frame(results(dds,contrast=c("Trait","OA","Control")))
  res$Gene <- rownames(res);res$Zone <- z
  res <- res[order(res$padj),]
  write.csv(res,file.path(OUT_DIR,paste0("GSE254844_OA_vs_Control_",z,"_poscounts.csv")),row.names=FALSE)
  result_summary[[z]] <- data.frame(Zone=z,N_microregions=ncol(sub_c),
    Up=sum(!is.na(res$padj)&res$padj<.05&res$log2FoldChange>=1),
    Down=sum(!is.na(res$padj)&res$padj<.05&res$log2FoldChange<= -1))
  # OA versus control has one control donor, limiting donor-level inference.
}
write.csv(do.call(rbind,result_summary),file.path(OUT_DIR,"GSE254844_DE_poscounts_Summary.csv"),row.names=FALSE)
# Within-donor WB vs NWB design, separately for each zone.
for (z in zones) {
  idx <- as.character(meta$Zone)==z & as.character(meta$Trait)=="OA"
  sub_m <- meta[idx,,drop=FALSE]
  sub_m$Patient_ID <- factor(sub_m$Patient_ID)
  sub_m$Condition <- factor(as.character(sub_m$Condition),levels=c("NWB","WB"))
  sub_c <- as.matrix(counts[,idx,drop=FALSE])
  keep <- rowSums(sub_c)>0
  dds <- DESeqDataSetFromMatrix(countData=sub_c[keep,,drop=FALSE],colData=sub_m,design=~Patient_ID+Condition)
  dds <- DESeq(dds,sfType="poscounts",quiet=TRUE)
  res <- as.data.frame(results(dds,contrast=c("Condition","WB","NWB")))
  res$Gene <- rownames(res);res$Zone <- z
  write.csv(res,file.path(OUT_DIR,paste0("GSE254844_WB_vs_NWB_",z,"_poscounts_paired.csv")),row.names=FALSE)
}
