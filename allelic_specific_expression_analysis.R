library(dplyr)
library(edgeR)

data_filter <- read.csv("expression.csv", header = TRUE, row.names = 1)
group <- select(data_filter, contains("_AF"))
cat("\n=== AF tissue ===\n")
cat("Total samples:", ncol(group), "\n")

# Check parental samples
mat_samples <- grep("_MAT", colnames(group), value = TRUE)
pat_samples <- grep("_PAT", colnames(group), value = TRUE)
cat("Maternal:", length(mat_samples), " Paternal:", length(pat_samples), "\n")
if (length(mat_samples) == 0 || length(pat_samples) == 0) stop("Insufficient parental samples")

y <- DGEList(counts = group, genes = rownames(group))
isexpr <- rowSums(cpm(y) > 1) >= 6
y <- y[isexpr, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
analyze_cross <- function(cross_group, cross_name) {
  if (ncol(cross_group) == 0) {
    cat(cross_name, "cross: no data\n")
    return(data.frame())
  }
  y_cross <- DGEList(counts = cross_group, genes = rownames(cross_group))
  y_cross <- y_cross[isexpr, , keep.lib.sizes = FALSE]
  y_cross <- calcNormFactors(y_cross, method = "TMM")
  
  po <- ifelse(grepl("_MAT", colnames(cross_group)), 1, 0)
  sex <- ifelse(grepl("_F_", colnames(cross_group)), 1, 0)
  design <- model.matrix(~PO + Sex, data = data.frame(PO = po, Sex = sex))
  
  y_cross <- estimateDisp(y_cross, design, robust = TRUE)
  fit <- glmFit(y_cross, design)
  coef_idx <- grep("PO", colnames(design))
  if (length(coef_idx) == 0) return(data.frame())
  
  lrt <- glmLRT(fit, coef = coef_idx[1])
  res <- lrt$table
  res$FDR <- p.adjust(res$PValue, method = "fdr")
  res_sig <- filter(res, FDR < 0.05, abs(logFC) >= 1)
  
  # Write outputs
  write.csv(res, paste0("AF_AllGenes_", cross_name, ".csv"), quote = FALSE, row.names = TRUE)
  write.csv(res_sig, paste0("AF_FDR005_FC1_", cross_name, "_removed.csv"), quote = FALSE, row.names = TRUE)
  cat(cross_name, "cross significant:", nrow(res_sig), "\n")
  return(res_sig)
}

lh_group <- group %>% select(starts_with("LH_"))
hl_group <- group %>% select(starts_with("HL_"))
group_lh_final <- analyze_cross(lh_group, "LH")
group_hl_final <- analyze_cross(hl_group, "HL")

if (nrow(group_lh_final) == 0 || nrow(group_hl_final) == 0) {
  stop("No significant genes in one or both crosses")
}

common_genes <- intersect(rownames(group_lh_final), rownames(group_hl_final))
cat("\nGenes significant in both crosses:", length(common_genes), "\n")
if (length(common_genes) == 0) stop("No common significant genes")

# Extract expression data for common genes
lh_common <- lh_group %>% filter(rownames(lh_group) %in% common_genes)
hl_common <- hl_group %>% filter(rownames(hl_group) %in% common_genes)

pseudocount <- 0.5

# Parental effect: paternal/maternal in each cross
lh_pat <- lh_common %>% select(contains("_PAT"))
lh_mat <- lh_common %>% select(contains("_MAT"))
hl_pat <- hl_common %>% select(contains("_PAT"))
hl_mat <- hl_common %>% select(contains("_MAT"))

ratio_lh <- (rowMeans(lh_pat, na.rm=TRUE) + pseudocount) / (rowMeans(lh_mat, na.rm=TRUE) + pseudocount)
ratio_hl <- (rowMeans(hl_pat, na.rm=TRUE) + pseudocount) / (rowMeans(hl_mat, na.rm=TRUE) + pseudocount)
parent_effect <- common_genes[((ratio_lh > 1 & ratio_hl > 1) | (ratio_lh < 1 & ratio_hl < 1))]

# Breed effect: L/H in each cross (L is maternal in HL, paternal in LH)
ratio_breed_lh <- (rowMeans(lh_pat, na.rm=TRUE) + pseudocount) / (rowMeans(lh_mat, na.rm=TRUE) + pseudocount)  
ratio_breed_hl <- (rowMeans(hl_mat, na.rm=TRUE) + pseudocount) / (rowMeans(hl_pat, na.rm=TRUE) + pseudocount)  
breed_effect <- common_genes[((ratio_breed_lh > 1 & ratio_breed_hl > 1) | (ratio_breed_lh < 1 & ratio_breed_hl < 1))]

cat("\n=== Summary ===\n")
cat("Parental effect genes:", length(parent_effect), "\n")
cat("Breed effect genes:", length(breed_effect), "\n")
cat("Overlap:", length(intersect(parent_effect, breed_effect)), "\n")

write.csv(data.frame(Gene = parent_effect), "AF_parent_effect_PO_genes.csv", quote = FALSE, row.names = FALSE)
write.csv(data.frame(Gene = breed_effect), "AF_breed_effect_PO_genes.csv", quote = FALSE, row.names = FALSE)
