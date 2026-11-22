# _______________________________________________________________________________
# Proteomics workflow: DEP with DeqMS processing 
# Conditions in file are:
# B = biofilm 
# BC = B + C = biofilm with cellobiose
# C = FL = freeliving culture
# CC = FL + C = freeliving with cellobiose
# _______________________________________________________________________________

#Loading Dataset and libraries
#protein_file <- "AP25015 grouped new databse altered settings.xlsx"
#annotation_file <- "annotation_LS.xlsx"

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(SummarizedExperiment)
library(DEP)
library(limma)
library(DEqMS)
library(ggplot2)
library(clusterProfiler)

prot_raw <- read_excel(protein_file, sheet = "Proteins")
anno_raw <- read_excel(annotation_file, sheet = 1)

anno <- anno_raw
if (!"Accession" %in% colnames(anno)) {
  colnames(anno)[1] <- "Accession"
}

anno <- anno %>%
  distinct(Accession, .keep_all = TRUE)

anno_keep <- c( "Accession", "Final_annoation", "Description", "Prokka_annoation",
  "KEGG old rr", "KEGG_ko", "seed_ortholog")

anno_keep <- anno_keep[anno_keep %in% colnames(anno)]

anno2 <- anno %>% select(all_of(anno_keep))

dat <- prot_raw %>% left_join(anno2, by = "Accession")

abundance_cols <- grep("^Abundance: F\\d+: Sample, ", colnames(dat), value = TRUE)
abundance_cols

sample_names <- abundance_cols

sample_info <- tibble( label = sample_names, condition = case_when(
    str_detect(label, "Sample, B$")  ~ "B", str_detect(label, "Sample, BC$") ~ "BC",
    str_detect(label, "Sample, C$")  ~ "C", str_detect(label, "Sample, CC$") ~ "CC",
    TRUE ~ NA_character_ )) %>% group_by(condition) %>%
  mutate(replicate = row_number()) %>% ungroup() %>%
  mutate(condition = factor(condition, levels = c("B", "BC", "C", "CC")))

# Using peptide counts for DEqMS; if unavailable, using PSM count
count_col <- if ("# Peptides" %in% colnames(dat)) "# Peptides" else "# PSMs"

analysis_df <- dat %>% mutate(gene_name = Accession, protein_id = Accession,
    n_peptides = .data[[count_col]]) %>%
  select( protein_id, gene_name, Accession, Description,
    any_of(c("Final_annoation", "Prokka_annoation", "KEGG old rr", "KEGG_ko", "seed_ortholog")),
    n_peptides, all_of(abundance_cols))

# Removing duplicated protein IDs if present
analysis_df <- analysis_df %>% distinct(protein_id, .keep_all = TRUE)

# Creating intensity matrix

expr <- analysis_df %>% select(all_of(abundance_cols)) %>% as.data.frame()

rownames(expr) <- analysis_df$protein_id

expr[expr == 0] <- NA

expr_log2 <- log2(expr)

# Defining structural absence
# Rule:
# If all 4 replicates of a condition are NA; structural absence (non-expression)
# If at least one replicate is present; other NAs are true missing values

cond_levels <- levels(sample_info$condition)

is_structural_absence <- matrix(FALSE, nrow = nrow(expr_log2), ncol = length(cond_levels))
rownames(is_structural_absence) <- rownames(expr_log2)
colnames(is_structural_absence) <- cond_levels

for (cond in cond_levels) {
  cols <- sample_info$label[sample_info$condition == cond]
  is_structural_absence[, cond] <- rowSums(is.na(expr_log2[, cols, drop = FALSE])) == length(cols)
}

# Keeping proteins that have enough replication
# Keep proteins detected in at least 3 replicate in at least one condition
keep <- sapply(
  levels(sample_info$condition),
  function(cond) {
    cols <- sample_info$label[
      sample_info$condition == cond
    ]
    
    rowSums(
      !is.na(expr_log2[, cols, drop = FALSE])
    ) >= 3
  }
)

# Retain proteins meeting the threshold in ANY treatment
keep <- rowSums(keep) >= 1

# Apply filtering
expr_log2 <- expr_log2[keep, , drop = FALSE]
analysis_df <- analysis_df[keep, , drop = FALSE]
is_structural_absence <- is_structural_absence[keep, , drop = FALSE]


# Impute only true missing values
# structural absence should stay NA and should not be imputed as low random values

expr_imp <- expr_log2

for (i in seq_len(nrow(expr_imp))) {
  for (cond in cond_levels) {
    cols <- sample_info$label[sample_info$condition == cond]
    if (is_structural_absence[i, cond]) next
    
    vals <- as.numeric(expr_imp[i, cols])
    
    # If some values exist and some are missing, impute as within-condition missing
    if (any(is.na(vals)) && any(!is.na(vals))) {
      observed <- vals[!is.na(vals)]
      imp_mean <- mean(observed) - 1.8 * sd(observed)
      imp_sd   <- 0.3 * sd(observed)
      if (is.na(imp_sd) || imp_sd == 0) imp_sd <- 0.2
      if (is.na(imp_mean)) imp_mean <- min(observed, na.rm = TRUE) - 1
      
      vals[is.na(vals)] <- rnorm(sum(is.na(vals)), mean = imp_mean, sd = imp_sd)
      expr_imp[i, cols] <- vals
    }
  }
}

# Global intensity normalization between arrays on log2 scale before limma/DEqMS

expr_norm <- normalizeBetweenArrays(as.matrix(expr_imp), method = "quantile")

# Restore structural absence as NA after normalization
for (cond in cond_levels) {
  cols <- sample_info$label[sample_info$condition == cond]
  idx <- which(is_structural_absence[, cond])
  if (length(idx) > 0) {
    expr_norm[idx, cols] <- NA
  }
}

expr_model <- expr_norm

sample_mins <- apply(expr_model, 2, function(x) min(x, na.rm = TRUE))
for (j in seq_len(ncol(expr_model))) {
  expr_model[is.na(expr_model[, j]), j] <- sample_mins[j] - 1
}

design <- model.matrix(~ 0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)
design

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LIMMA with DEqMS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Limma fit

fit <- lmFit(expr_model, design)

contrast.matrix <- makeContrasts(B_vs_C = B - C, BC_vs_CC = BC - CC,
  BC_vs_B = BC - B, CC_vs_C = CC - C, BC_vs_C = BC - C, CC_vs_B = CC - B,
  levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

fit2$count <- pep_count[rownames(fit2$coefficients)]

fit3 <- spectraCounteBayes(fit2)

# DEqMS using peptide/PSM counts

pep_count <- analysis_df$n_peptides
names(pep_count) <- analysis_df$protein_id
fit2$count <- pep_count[rownames(fit2$coefficients)]
fit3 <- spectraCounteBayes(fit2)

get_deqms_results <- function(fit_obj, coef_name, analysis_df) {
  res <- outputResult(fit_obj, coef_col = which(colnames(fit_obj$coefficients) == coef_name))
  res <- as.data.frame(res)
  res$protein_id <- rownames(res)
  
  out <- analysis_df %>%
    select(
      protein_id, Accession, gene_name, Description,
      any_of(c("Final_annoation", "Prokka_annoation", "KEGG old rr", "KEGG_ko", "seed_ortholog")),
      n_peptides
    ) %>%
    left_join(res, by = "protein_id") %>%
    arrange(adj.P.Val)
  
  return(out)
}

res_B_vs_C   <- get_deqms_results(fit3, "B_vs_C", analysis_df)
res_BC_vs_CC <- get_deqms_results(fit3, "BC_vs_CC", analysis_df)
res_BC_vs_B  <- get_deqms_results(fit3, "BC_vs_B", analysis_df)
res_CC_vs_C  <- get_deqms_results(fit3, "CC_vs_C", analysis_df)
res_BC_vs_C  <- get_deqms_results(fit3, "BC_vs_C", analysis_df)
res_CC_vs_B  <- get_deqms_results(fit3, "CC_vs_B", analysis_df)

# Significant proteins

sig_B_vs_C <- res_B_vs_C %>% filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) >= 0.5)
sig_BC_vs_CC <- res_BC_vs_CC %>% filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) >= 0.5)
sig_BC_vs_B <- res_BC_vs_B %>% filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) >= 0.5)
sig_CC_vs_C <- res_CC_vs_C %>% filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) >= 0.5)

# QC plots

pdf("QC_boxplot_normalized.pdf", width = 10, height = 6)
boxplot(expr_norm, las = 2, main = "Normalized log2 intensities", outline = FALSE)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

var_genes <- apply(expr_model, 1, var, na.rm = TRUE)
expr_pca <- expr_model[var_genes > 1e-6, ]
pca <- prcomp(t(expr_pca), scale. = TRUE)
pca_df <- data.frame( PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  condition = sample_info$condition, sample = sample_info$label)

condition_cols <- c("B"  = "#2d6c25", "BC" = "#72426A", "CC" = "#6fbbe6", "C"  = "grey20")

expr_pca <- as.matrix(expr_norm)
expr_pca[!is.finite(expr_pca)] <- NA
expr_pca <- expr_pca[rowSums(!is.na(expr_pca)) >= 2, , drop = FALSE]
for (i in seq_len(nrow(expr_pca))) {
  vals <- expr_pca[i, ]
  if (anyNA(vals)) {
    min_val <- min(vals, na.rm = TRUE)
    vals[is.na(vals)] <- min_val
    expr_pca[i, ] <- vals
  }
}

# Remove zero-variance proteins
rv <- apply(expr_pca, 1, var)
expr_pca <- expr_pca[rv > 0 & !is.na(rv), , drop = FALSE]

pca_res <- plot_pca_conditions(expr_pca, sample_info, condition_cols)
pca_res$plot

sum(is.na(expr_norm))
sum(is.nan(as.matrix(expr_norm)))
sum(is.infinite(as.matrix(expr_norm)))
range(expr_norm, na.rm = TRUE)

plot_pca_conditions <- function(expr_matrix, sample_info, condition_cols) {
  library(ggplot2)
  library(dplyr)
  library(vegan)
  
  sample_info <- sample_info %>%
    mutate(condition = factor(condition, levels = names(condition_cols)))
  
  expr_matrix <- as.matrix(expr_matrix)
  expr_matrix <- expr_matrix[, sample_info$label, drop = FALSE]
  
  # Remove rows with any non-finite values
  finite_rows <- apply(expr_matrix, 1, function(x) all(is.finite(x)))
  expr_matrix <- expr_matrix[finite_rows, , drop = FALSE]
  
  # Remove zero-variance rows
  row_var <- apply(expr_matrix, 1, var)
  expr_matrix <- expr_matrix[row_var > 0 & !is.na(row_var), , drop = FALSE]
  
  # Safety check
  if (nrow(expr_matrix) == 0) {
    stop("No proteins left for PCA after removing non-finite and zero-variance rows.")
  }
  
  pca <- prcomp(t(expr_matrix), scale. = TRUE)
  
  percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  
  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  ) %>%
    left_join(sample_info, by = c("Sample" = "label"))
  
  dist_mat <- dist(t(expr_matrix), method = "euclidean")
  permanova <- adonis2(dist_mat ~ condition, data = sample_info, permutations = 999)
  pval <- permanova$`Pr(>F)`[1]
  pval_txt <- ifelse(pval < 0.001, "PERMANOVA p < 0.001", paste0("PERMANOVA p = ", signif(pval, 3)))
  
 g <- ggplot(pca_df, aes(PC1, PC2, colour = condition, fill = condition)) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, show.legend = FALSE) +
    geom_point(size = 2, alpha = 0.9) +
    scale_colour_manual(values = condition_cols) +
    scale_fill_manual(values = condition_cols) +
    theme_bw(base_size = 13) +
    labs(title = "PCA", subtitle = pval_txt,
      x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
      y = paste0("PC2 (", round(percent_var[2], 1), "%)") ) +
    theme( panel.grid = element_blank(), legend.title = element_blank())
  
  list(plot = g, pca = pca, data = pca_df, permanova = permanova, matrix_used = expr_matrix)
}

pca_res <- plot_pca_conditions(expr_norm, sample_info, condition_cols)
pca_res$plot

prep_pca_matrix <- function(expr_matrix, sample_info_subset) {
  expr_pca <- as.matrix(expr_matrix)
  expr_pca <- expr_pca[, sample_info_subset$label, drop = FALSE]
  expr_pca[!is.finite(expr_pca)] <- NA
  expr_pca <- expr_pca[rowSums(!is.na(expr_pca)) >= 3, , drop = FALSE]
  for (i in seq_len(nrow(expr_pca))) {
    vals <- expr_pca[i, ]
    if (anyNA(vals)) {
      min_val <- min(vals, na.rm = TRUE)
      vals[is.na(vals)] <- min_val
      expr_pca[i, ] <- vals
    }
  }

  rv <- apply(expr_pca, 1, var, na.rm = TRUE)
  expr_pca <- expr_pca[rv > 0 & !is.na(rv), , drop = FALSE]
  
  expr_pca
}

dev.off()

plot_pca_subset <- function(expr_matrix, sample_info_subset, condition_cols, title = "PCA") {
  expr_matrix <- expr_matrix[, sample_info_subset$label, drop = FALSE]
  
  pca <- prcomp(t(expr_matrix), scale. = TRUE)
  percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  pca_df <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]) %>% left_join(sample_info_subset, by = c("Sample" = "label"))
  
  dist_mat <- dist(t(expr_matrix), method = "euclidean")
  permanova <- adonis2(dist_mat ~ condition, data = sample_info_subset, permutations = 999)
  pval <- permanova$`Pr(>F)`[1]
  pval_txt <- ifelse(pval < 0.001, "PERMANOVA p < 0.001", paste0("PERMANOVA p = ", signif(pval, 3)))
  
  ggplot(pca_df, aes(PC1, PC2, colour = condition, fill = condition)) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, show.legend = FALSE) +
    geom_point(size = 2, alpha = 0.95) +
    scale_colour_manual(values = condition_cols) +
    scale_fill_manual(values = condition_cols) +
    theme_bw(base_size = 13) +
    labs(
      title = title,
      subtitle = pval_txt,
      x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
      y = paste0("PC2 (", round(percent_var[2], 1), "%)")
    ) +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank()
    )
}

ggplot(pca_df, aes(PC1, PC2, colour = condition, fill = condition)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, show.legend = FALSE) +
  geom_point(size = 2, alpha = 0.95) +
  scale_colour_manual(values = condition_cols) +
  scale_fill_manual(values = condition_cols) +
  theme_bw(base_size = 13) +
  theme( panel.grid = element_blank(), legend.title = element_blank() )

# CC vs C only
sample_info_CC_C <- sample_info %>%
  filter(condition %in% c("CC", "C")) %>%
  mutate(condition = factor(condition, levels = c("C", "CC")))
expr_pca_CC_C <- prep_pca_matrix(expr_norm, sample_info_CC_C)

 plot_pca_subset(
  expr_matrix = expr_pca_CC_C,
  sample_info_subset = sample_info_CC_C,
  condition_cols = condition_cols[c("C", "CC")],
  title = "PCA: C vs CC")

# B vs BC only
sample_info_B_BC <- sample_info %>%
  filter(condition %in% c("B", "BC")) %>%
  mutate(condition = factor(condition, levels = c("B", "BC")))

expr_pca_B_BC <- prep_pca_matrix(expr_norm, sample_info_B_BC)

 plot_pca_subset(
  expr_matrix = expr_pca_B_BC,
  sample_info_subset = sample_info_B_BC,
  condition_cols = condition_cols[c("B", "BC")],
  title = "PCA: B vs BC")

#PERMANOVA RESULTS

run_permanova <- function(expr_matrix, sample_info_subset, method = "euclidean") {
 expr_matrix <- expr_matrix[, sample_info_subset$label, drop = FALSE]
 dist_mat <- dist(t(expr_matrix), method = method)
 perm <- adonis2(dist_mat ~ condition, data = sample_info_subset, permutations = 999)
 res <- data.frame( term = rownames(perm), Df = perm$Df, SumOfSqs = perm$SumOfSqs,
    R2 = perm$R2, F = perm$F, p_value = perm$`Pr(>F)`)
  res
}


perm_CC_C <- run_permanova(expr_pca_CC_C, sample_info_CC_C)
perm_CC_C

perm_B_BC <- run_permanova(expr_pca_B_BC, sample_info_B_BC)
perm_B_BC

perm_summary <- bind_rows(perm_CC_C %>% mutate(comparison = "C_vs_CC"),
  perm_B_BC %>% mutate(comparison = "B_vs_BC"))
perm_summary 

perm_all <- run_permanova(expr_pca, sample_info)
perm_all


plot_volcano <- function(df, title) {
  df$y <- -log10(df$adj.P.Val)
  df$col <- ifelse(df$adj.P.Val < 0.05 & abs(df$logFC) >= 0.5, "green3", "grey70")
  
  ggplot(df, aes(logFC, y, colour = col)) +
    geom_point() +
    scale_colour_identity() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_y_continuous(limits = c(0, NA)) +
    theme_bw() +
    labs(title = title, x = "log2FC", y = "-log10 adj.P")
}

plot_volcano(res_B_vs_C,   "B vs C")
plot_volcano(res_BC_vs_CC,"BC vs CC")
plot_volcano(res_BC_vs_B, "BC vs B")
plot_volcano(res_CC_vs_C, "CC vs C")
plot_volcano(res_BC_vs_C, "BC vs C")
plot_volcano(res_CC_vs_B, "CC vs B")

#Adding significance code

add_significance_column <- function(res_df, expr_log2, sample_info, cond1, cond2, p_cutoff = 0.05) {
  
  cols1 <- sample_info$label[sample_info$condition == cond1]
  cols2 <- sample_info$label[sample_info$condition == cond2]
  
  present1 <- rowSums(!is.na(expr_log2[res_df$protein_id, cols1, drop = FALSE])) > 0
  present2 <- rowSums(!is.na(expr_log2[res_df$protein_id, cols2, drop = FALSE])) > 0
  
  res_df$significance <- ifelse(
    !present2 & present1, paste("Exclusively expressed in", cond1),
    ifelse(
      !present1 & present2, paste("Exclusively expressed in", cond2),
      ifelse(
        is.na(res_df$adj.P.Val) | res_df$adj.P.Val >= p_cutoff, "Non-significant",
        ifelse(res_df$logFC > 0, paste("More in", cond1), paste("More in", cond2))
      )
    )
  )
  
  res_df
}

#Apply to all comparisons

res_B_vs_C   <- add_significance_column(res_B_vs_C,   expr_log2, sample_info, "B",  "C")
res_BC_vs_CC <- add_significance_column(res_BC_vs_CC, expr_log2, sample_info, "BC", "CC")
res_BC_vs_B  <- add_significance_column(res_BC_vs_B,  expr_log2, sample_info, "BC", "B")
res_CC_vs_C  <- add_significance_column(res_CC_vs_C,  expr_log2, sample_info, "CC", "C")
res_BC_vs_C  <- add_significance_column(res_BC_vs_C,  expr_log2, sample_info, "BC", "C")
res_CC_vs_B  <- add_significance_column(res_CC_vs_B,  expr_log2, sample_info, "CC", "B")

sig_B_vs_C   <- subset(res_B_vs_C,   significance != "Non-significant")
sig_BC_vs_CC <- subset(res_BC_vs_CC, significance != "Non-significant")
sig_BC_vs_B  <- subset(res_BC_vs_B,  significance != "Non-significant")
sig_CC_vs_C  <- subset(res_CC_vs_C,  significance != "Non-significant")
sig_BC_vs_C  <- subset(res_BC_vs_C,  significance != "Non-significant")
sig_CC_vs_B  <- subset(res_CC_vs_B,  significance != "Non-significant")


# _______________________________________________________________________________
# Proteomics workflow: DEP with DeqMS processing 
# Conditions in file are:
# pH 5 pH 7 and pH 6
# _______________________________________________________________________________

dat_57 <- read_excel(file_57, sheet = "Proteins")
dat_6  <- read_excel(file_6,  sheet = "Proteins")
anno_raw <- read_excel(annotation_file, sheet = 1)

anno <- anno_raw
if (!"Accession" %in% colnames(anno)) colnames(anno)[1] <- "Accession"
abund_57 <- grep("^Abundance: F\\d+: Sample, ", colnames(dat_57), value = TRUE)
abund_6  <- grep("^Abundance: F\\d+: Sample, ", colnames(dat_6),  value = TRUE)

count_col_57 <- if ("# Peptides" %in% colnames(dat_57)) "# Peptides" else "# PSMs"
count_col_6  <- if ("# Peptides" %in% colnames(dat_6))  "# Peptides" else "# PSMs"

count_col_57 <- grep("Peptide|PSM", colnames(dat_57), value = TRUE)[1]
count_col_6  <- grep("Peptide|PSM", colnames(dat_6),  value = TRUE)[1]

print(count_col_57)
print(count_col_6)

if (is.na(count_col_57)) stop("No peptide/PSM column found in dat_57")
if (is.na(count_col_6))  stop("No peptide/PSM column found in dat_6")

dat_57_small <- dat_57 %>%
  select(Accession, Description, all_of(count_col_57), all_of(abund_57))

colnames(dat_57_small)[colnames(dat_57_small) == count_col_57] <- "n_peptides_57"

dat_6_small <- dat_6 %>%
  select(Accession, Description, all_of(count_col_6), all_of(abund_6))

colnames(dat_6_small)[colnames(dat_6_small) == count_col_6] <- "n_peptides_6"

dat <- full_join(dat_57_small, dat_6_small, by = "Accession", suffix = c("_57", "_6")) %>%
  mutate(Description = coalesce(Description_57, Description_6),
    n_peptides = coalesce(n_peptides_57, n_peptides_6)) %>%
  select(Accession, Description, n_peptides, all_of(abund_57), all_of(abund_6))

anno_keep <- c("Accession", "Final_annoation", "Description", "Prokka_annoation",
  "KEGG old rr", "KEGG_ko", "KEGG_ko.1", "KEGG_Pathway", "COG_category", "seed_ortholog")
anno_keep <- anno_keep[anno_keep %in% colnames(anno)]

dat <- dat %>%
  left_join(anno %>% select(all_of(anno_keep)), by = "Accession")

sample_names <- c(abund_57, abund_6)

sample_info <- tibble(label = sample_names,condition = case_when(
    str_detect(label, "Sample, 5$") ~ "pH5", str_detect(label, "Sample, 7$") ~ "pH7",
    str_detect(label, "Sample, C$") ~ "pH6", TRUE ~ NA_character_ )) %>%
  group_by(condition) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  mutate(condition = factor(condition, levels = c("pH5", "pH6", "pH7")))

analysis_df <- dat %>%
  mutate(protein_id = Accession, gene_name = Accession) %>%
  distinct(protein_id, .keep_all = TRUE)


expr <- analysis_df %>%
  select(all_of(sample_names)) %>%
  as.data.frame()

rownames(expr) <- analysis_df$protein_id
expr[expr == 0] <- NA

expr_log2 <- log2(expr)

cond_levels <- levels(sample_info$condition)

is_structural_absence <- matrix(FALSE, nrow = nrow(expr_log2), ncol = length(cond_levels))
rownames(is_structural_absence) <- rownames(expr_log2)
colnames(is_structural_absence) <- cond_levels

for (cond in cond_levels) {
  cols <- sample_info$label[sample_info$condition == cond]
  is_structural_absence[, cond] <- rowSums(is.na(expr_log2[, cols, drop = FALSE])) == length(cols)
}

keep <- sapply(
  levels(sample_info$condition),
  function(cond) {
    cols <- sample_info$label[
      sample_info$condition == cond
    ]
    
    rowSums(
      !is.na(expr_log2[, cols, drop = FALSE])
    ) >= 3
  }
)

# Retain proteins meeting the threshold in ANY treatment
keep <- rowSums(keep) >= 1

# Apply filtering
expr_log2 <- expr_log2[keep, , drop = FALSE]
analysis_df <- analysis_df[keep, , drop = FALSE]
is_structural_absence <- is_structural_absence[keep, , drop = FALSE]


expr_imp <- expr_log2

for (i in seq_len(nrow(expr_imp))) {
  for (cond in cond_levels) {
    cols <- sample_info$label[sample_info$condition == cond]
    
    if (is_structural_absence[i, cond]) next
    
    vals <- as.numeric(expr_imp[i, cols])
    
    if (any(is.na(vals)) && any(!is.na(vals))) {
      observed <- vals[!is.na(vals)]
      imp_mean <- mean(observed) - 1.8 * sd(observed)
      imp_sd   <- 0.3 * sd(observed)
      
      if (is.na(imp_sd) || imp_sd == 0) imp_sd <- 0.2
      if (is.na(imp_mean)) imp_mean <- min(observed, na.rm = TRUE) - 1
      
      vals[is.na(vals)] <- rnorm(sum(is.na(vals)), mean = imp_mean, sd = imp_sd)
      expr_imp[i, cols] <- vals
    }
  }
}

expr_norm <- normalizeBetweenArrays(as.matrix(expr_imp), method = "quantile")

for (cond in cond_levels) {
  cols <- sample_info$label[sample_info$condition == cond]
  idx <- which(is_structural_absence[, cond])
  if (length(idx) > 0) expr_norm[idx, cols] <- NA
}

expr_model <- expr_norm

sample_mins <- apply(expr_model, 2, function(x) min(x, na.rm = TRUE))
for (j in seq_len(ncol(expr_model))) {
  expr_model[is.na(expr_model[, j]), j] <- sample_mins[j] - 1
}

design <- model.matrix(~ 0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)

fit <- lmFit(expr_model, design)

contrast.matrix <- makeContrasts(
  pH5_vs_pH7 = pH5 - pH7,
  pH5_vs_pH6 = pH5 - pH6,
  pH7_vs_pH6 = pH7 - pH6,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

pep_count <- as.numeric(analysis_df$n_peptides)
names(pep_count) <- analysis_df$protein_id
fit2$count <- pep_count[rownames(fit2$coefficients)]

fit3 <- spectraCounteBayes(fit2)

get_deqms_results <- function(fit_obj, coef_name, analysis_df) {
  coef_idx <- which(colnames(fit_obj$coefficients) == coef_name)
  
  res <- outputResult(fit_obj, coef_col = coef_idx)
  res <- as.data.frame(res)
  res$protein_id <- rownames(res)
  
  analysis_df %>%
    select(
      protein_id, Accession, gene_name, Description, n_peptides,
      any_of(c("Final_annoation", "Prokka_annoation", "KEGG old rr",
               "KEGG_ko", "KEGG_ko.1", "KEGG_Pathway", "COG_category", "seed_ortholog"))
    ) %>%
    left_join(res, by = "protein_id") %>%
    arrange(adj.P.Val)
}

res_pH5_vs_pH7 <- get_deqms_results(fit3, "pH5_vs_pH7", analysis_df)
res_pH5_vs_pH6 <- get_deqms_results(fit3, "pH5_vs_pH6", analysis_df)
res_pH7_vs_pH6 <- get_deqms_results(fit3, "pH7_vs_pH6", analysis_df)

add_significance_column <- function(res_df, expr_log2, sample_info, cond1, cond2,
                                    p_cutoff = 0.05, fc_cutoff = 0.5) {
  
  cols1 <- sample_info$label[sample_info$condition == cond1]
  cols2 <- sample_info$label[sample_info$condition == cond2]
  
  present1 <- rowSums(!is.na(expr_log2[res_df$protein_id, cols1, drop = FALSE])) > 0
  present2 <- rowSums(!is.na(expr_log2[res_df$protein_id, cols2, drop = FALSE])) > 0
  
  res_df$significance <- ifelse(
    !present2 & present1, paste("Exclusively expressed in", cond1),
    ifelse(
      !present1 & present2, paste("Exclusively expressed in", cond2),
      ifelse(
        is.na(res_df$adj.P.Val) | res_df$adj.P.Val >= p_cutoff | abs(res_df$logFC) < fc_cutoff,
        "Non-significant",
        ifelse(res_df$logFC > 0, paste("More in", cond1), paste("More in", cond2))
      )
    )
  )
  
  res_df
}

res_pH5_vs_pH7 <- add_significance_column(res_pH5_vs_pH7, expr_log2, sample_info, "pH5", "pH7")
res_pH5_vs_pH6 <- add_significance_column(res_pH5_vs_pH6, expr_log2, sample_info, "pH5", "pH6")
res_pH7_vs_pH6 <- add_significance_column(res_pH7_vs_pH6, expr_log2, sample_info, "pH7", "pH6")

sig_pH5_vs_pH7 <- subset(res_pH5_vs_pH7, significance != "Non-significant")
sig_pH5_vs_pH6 <- subset(res_pH5_vs_pH6, significance != "Non-significant")
sig_pH7_vs_pH6 <- subset(res_pH7_vs_pH6, significance != "Non-significant")

# PCA

condition_cols <- c("pH5" = "#f69dad","pH6" = "grey40","pH7" = "#a9abe0")

expr_pca <- as.matrix(expr_norm)
expr_pca[!is.finite(expr_pca)] <- NA
expr_pca <- expr_pca[rowSums(!is.na(expr_pca)) >= 2, , drop = FALSE]

for (i in seq_len(nrow(expr_pca))) {
  vals <- expr_pca[i, ]
  if (anyNA(vals)) vals[is.na(vals)] <- min(vals, na.rm = TRUE)
  expr_pca[i, ] <- vals
}

rv <- apply(expr_pca, 1, var, na.rm = TRUE)
expr_pca <- expr_pca[rv > 0 & !is.na(rv), , drop = FALSE]

pca <- prcomp(t(expr_pca), scale. = TRUE)
percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>%
  left_join(sample_info, by = c("Sample" = "label"))

dist_mat <- dist(t(expr_pca), method = "euclidean")
perm_all <- adonis2(dist_mat ~ condition, data = sample_info, permutations = 999)
pval <- perm_all$`Pr(>F)`[1]
pval_txt <- ifelse(pval < 0.001, "PERMANOVA p < 0.001", paste0("PERMANOVA p = ", signif(pval, 3)))

ggplot(pca_df, aes(PC1, PC2, colour = condition, fill = condition)) +
  geom_point(size = 3.5, alpha = 0.95) +
  scale_colour_manual(values = condition_cols) +
  scale_fill_manual(values = condition_cols) +
  theme_bw(base_size = 13) +
  labs(
    title = "PCA: pH proteome comparison",
    subtitle = pval_txt,
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )

bd <- betadisper(dist_mat, sample_info$condition)
anova(bd)

dist_mat <- dist(t(expr_pca), method = "euclidean")

# PERMANOVA

perm_all <- adonis2(dist_mat ~ condition, data = sample_info, permutations = 999)
perm_all

perm_all_df <- data.frame(
  term = rownames(perm_all),
  Df = perm_all$Df,
  SumOfSqs = perm_all$SumOfSqs,
  R2 = perm_all$R2,
  F = perm_all$F,
  p_value = perm_all$`Pr(>F)`
)

perm_all_df

#Volcano plots 

plot_volcano <- function(df, title) {
  df$y <- -log10(df$adj.P.Val)
  df$col <- ifelse(df$adj.P.Val < 0.05 & abs(df$logFC) >= 0.5, "green3", "grey60")
  
  ggplot(df, aes(logFC, y, colour = col)) +
    geom_point() +
    scale_colour_identity() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_y_continuous(limits = c(0, NA)) +
    theme_bw() +
    labs(title = title, x = "log2FC", y = "-log10 adj.P")
}

plot_volcano(res_pH5_vs_pH7, "pH5 vs pH7")
plot_volcano(res_pH5_vs_pH6, "pH5 vs pH6")
plot_volcano(res_pH7_vs_pH6, "pH7 vs pH6")

add_significance_column <- function(res_df, expr_log2, sample_info, cond1, cond2,
                                    p_cutoff = 0.05, fc_cutoff = 0.5) {
  
  cols1 <- sample_info$label[sample_info$condition == cond1]
  cols2 <- sample_info$label[sample_info$condition == cond2]
  
  present1 <- rowSums(!is.na(expr_log2[res_df$protein_id, cols1, drop = FALSE])) > 0
  present2 <- rowSums(!is.na(expr_log2[res_df$protein_id, cols2, drop = FALSE])) > 0
  
  res_df$significance <- ifelse(
    !present2 & present1, paste("Exclusively expressed in", cond1),
    ifelse(
      !present1 & present2, paste("Exclusively expressed in", cond2),
      ifelse(
        is.na(res_df$adj.P.Val) | res_df$adj.P.Val >= p_cutoff | abs(res_df$logFC) < fc_cutoff,
        "Non-significant",
        ifelse(res_df$logFC > 0, paste("More in", cond1), paste("More in", cond2))
      )
    )
  )
  
  res_df
}

res_pH5_vs_pH7 <- add_significance_column(res_pH5_vs_pH7, expr_log2, sample_info, "pH5", "pH7")
res_pH5_vs_pH6 <- add_significance_column(res_pH5_vs_pH6, expr_log2, sample_info, "pH5", "pH6")
res_pH7_vs_pH6 <- add_significance_column(res_pH7_vs_pH6, expr_log2, sample_info, "pH7", "pH6")

sig_pH5_vs_pH7 <- subset(res_pH5_vs_pH7, significance != "Non-significant")
sig_pH5_vs_pH6 <- subset(res_pH5_vs_pH6, significance != "Non-significant")
sig_pH7_vs_pH6 <- subset(res_pH7_vs_pH6, significance != "Non-significant")


