# Necessary Packages --------------------------------------------------------
library(ballgown)
library(NormalyzerDE)
library(FactoMineR)
library(factoextra)
library(limma)
library(clusterProfiler)
library(org.At.tair.db)
library(VennDiagram)
library(WGCNA)
library(cluster)
library(ggplot2)
library(dplyr)

# Set Paths -----------------------------------------------------------------
input_path <- "~/write/your/path"
output_path <- "~/write/your/output_path"

# Create necessary directories if they don't exist
dirs_to_create <- c("pca", "resultados/ecofun1/t_vs_m", "resultados/ecofun1/ex_vs_ant_ex", 
                    "resultados/ecofun2/t_vs_m", "resultados/ecofun2/ex_vs_ant_ex")
for (dir in dirs_to_create) {
  dir.create(file.path(output_path, dir), recursive = TRUE, showWarnings = FALSE)
}

# Set Working Directory
setwd(input_path)

# First Year Analysis -------------------------------------------------------

# Load experimental design
experimental_design <- read.csv(file.path(input_path, "ecofun_experimental_design.csv"))

# Load ballgown data
ecofun_data <- ballgown(dataDir = ".", samplePattern = "sample", pData = experimental_design)

# Extract gene expression data
gene_expression <- gexpr(ecofun_data)
colnames(gene_expression) <- c(
  "ex1velm_1","ex1velm_2","ex1velm_3","ex1velt_1","ex1velt_2","ex1velt_3",
  "ex2velm_1","ex2velm_2","ex2velm_3","ex2velm_4","ex2velt_1","ex2velt_2",
  "ex2velt_3","ex3velm_1","ex3velm_2","ex3velt_1","ex3velt_2","ex3velt_3",
  "ex4velm_1","ex4velm_2","ex4velm_3","ex4velt_1","ex4velt_2","ex4velt_3",
  "ex1bonm_1","ex1bonm_2","ex1bonm_3","ex1bont_1","ex1bont_2","ex1bont_3",
  "ex2bonm_1","ex2bonm_2","ex2bonm_3","ex2bont_1","ex2bont_2","ex2bont_3",
  "ex3bonm_1","ex3bonm_2","ex3bonm_3","ex3bont_1","ex3bont_2","ex3bont_3",
  "ex1mojm_1","ex1mojm_2","ex1mojm_3","ex1mojt_1","ex1mojt_2","ex1mojt_3",
  "ex2mojm_1","ex2mojm_2","ex2mojm_3","ex2mojt_1","ex2mojt_2","ex2mojt_3",
  "ex3mojm_1","ex3mojm_2","ex3mojm_3","ex3mojt_1","ex3mojt_2","ex3mojt_3",
  "ex1col0m_1","ex1col0m_2","ex1col0m_3","ex1col0t_1","ex1col0t_2","ex1col0t_3",
  "ex2col0m_1","ex2col0m_2","ex2col0m_3","ex2col0t_1","ex2col0t_2","ex2col0t_3",
  "ex3col0m_1","ex3col0m_2","ex3col0m_3","ex3col0t_1","ex3col0t_2","ex3col0t_3",
  "ex4col0m_1","ex4col0m_2","ex4col0m_3","ex4col0t_1","ex4col0t_2","ex4col0t_3")

# Write gene expression data to file
write.table(gene_expression, file.path(output_path, "ecofun_gene_expression.tsv"), quote = FALSE, sep = "\t")

# Read the gene expression data
ecofun_gene_expression <- read.table(file.path(output_path, "ecofun_gene_expression.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Boxplot of raw data
png(filename = file.path(output_path, "raw_data_boxplot.png"), width = 2000)
boxplot(log2(ecofun_gene_expression + 1), col = rainbow(ncol(ecofun_gene_expression)), ylab = "log2(FPKM + 1)", cex.lab = 1.5, las = 2)
dev.off()

# Remove degraded sample (column 10)
ecofun_gene_expression <- ecofun_gene_expression[, -10]

# Normalization
ecofun_gene_expression_df <- data.frame(geneID = rownames(ecofun_gene_expression), ecofun_gene_expression + 1)
write.table(ecofun_gene_expression_df, file.path(output_path, "ecofun_gene_expression_df.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

# Update experimental design
design <- data.frame(sample = colnames(ecofun_gene_expression), group = experimental_design$ecotype[-10])
write.table(design, file.path(output_path, "ecofun_design.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

# Perform normalization with NormalyzerDE
normalyzer(jobName = "ecofun", designPath = file.path(output_path, "ecofun_design.tsv"), dataPath = file.path(output_path, "ecofun_gene_expression_df.tsv"), outputDir = output_path)

# Read normalized and log-transformed gene expression data
log_normalized_ecofun_gene_expression <- read.table(file.path(output_path, "ecofun/Quantile-normalized.txt"), header = TRUE, stringsAsFactors = FALSE)
geneID <- log_normalized_ecofun_gene_expression$geneID
log_normalized_ecofun_gene_expression <- as.matrix(log_normalized_ecofun_gene_expression[, -1])
rownames(log_normalized_ecofun_gene_expression) <- geneID

# Differential Gene Expression Analysis -------------------------------------
DEGs <- function(year, contrast, output_name) {
  # Set file paths
  root <- output_path
  file_norm <- paste0(root, "/datos/Quantile-normalized_ecofun", year, ".txt")
  file_exp_design <- paste0(root, "/datos/exp_design_e", year, ".csv")
  
  if (year == 1) {
    root2 <- "resultados/ecofun1/"
  } else {
    root2 <- "resultados/ecofun2/"
  }
  
  # Load and clean normalized gene expression data
  log_n_exp <- read.table(file_norm, header = TRUE, stringsAsFactors = FALSE)
  geneID <- log_n_exp$geneID
  log_n_exp <- as.matrix(log_n_exp[, -1])
  rownames(log_n_exp) <- geneID
  
  # Load experimental design
  exp_design <- as.matrix(read.csv(file_exp_design))
  
  # Fit linear model
  linear_fit <- lmFit(log_n_exp, exp_design)
  
  # Define contrast matrices
  contrast_matrix1 <- makeContrasts(
    ex1velt - ex1velm, ex2velt - ex2velm, ex3velt - ex3velm, ex4velt - ex4velm,
    ex1bont - ex1bonm, ex2bont - ex2bonm, ex3bont - ex3bonm, ex1mojt - ex1mojm,
    ex2mojt - ex2mojm, ex3mojt - ex3mojm, ex1col0t - ex1col0m, ex2col0t - ex2col0m,
    ex3col0t - ex3col0m, ex4col0t - ex4col0m,
    levels = c("ex1velm", "ex1velt", "ex2velm", "ex2velt", "ex3velm", "ex3velt", "ex4velm", "ex4velt",
               "ex1bonm", "ex1bont", "ex2bonm", "ex2bont", "ex3bonm", "ex3bont", "ex1mojm", "ex1mojt",
               "ex2mojm", "ex2mojt", "ex3mojm", "ex3mojt", "ex1col0m", "ex1col0t", "ex2col0m", "ex2col0t",
               "ex3col0m", "ex3col0t", "ex4col0m", "ex4col0t")
  )
  
  contrast_matrix2 <- makeContrasts(
    ex2velm - ex1velm, ex2velt - ex1velt, ex3velm - ex1velm, ex3velt - ex1velt, ex3velm - ex2velm, ex3velt - ex2velt,
    ex4velm - ex1velm, ex4velt - ex1velt, ex4velm - ex2velm, ex4velt - ex2velt, ex4velm - ex3velm, ex4velt - ex3velt,
    ex2bonm - ex1bonm, ex2bont - ex1bont, ex3bonm - ex1bonm, ex3bont - ex1bont, ex3bonm - ex2bonm, ex3bont - ex2bont,
    ex2mojm - ex1mojm, ex2mojt - ex1mojt, ex3mojm - ex1mojm, ex3mojt - ex1mojt, ex3mojm - ex2mojm, ex3mojt - ex2mojt,
    ex2col0m - ex1col0m, ex2col0t - ex1col0t, ex3col0m - ex1col0m, ex3col0t - ex1col0t, ex3col0m - ex2col0m, ex3col0t - ex2col0t,
    ex4col0m - ex1col0m, ex4col0t - ex1col0t, ex4col0m - ex2col0m, ex4col0t - ex2col0t, ex4col0m - ex3col0m, ex4col0t - ex3col0t,
    levels = c("ex1velm", "ex1velt", "ex2velm", "ex2velt", "ex3velm", "ex3velt", "ex4velm", "ex4velt",
               "ex1bonm", "ex1bont", "ex2bonm", "ex2bont", "ex3bonm", "ex3bont", "ex1mojm", "ex1mojt",
               "ex2mojm", "ex2mojt", "ex3mojm", "ex3mojt", "ex1col0m", "ex1col0t", "ex2col0m", "ex2col0t",
               "ex3col0m", "ex3col0t", "ex4col0m", "ex4col0t")
  )
  
  # Select the contrast matrix
  if (contrast == 1) {
    contrasts_linear_fit <- contrasts.fit(linear_fit, contrasts = contrast_matrix1)
    root3 <- "t_vs_m/"
  } else {
    contrasts_linear_fit <- contrasts.fit(linear_fit, contrasts = contrast_matrix2)
    root3 <- "ex_vs_ant_ex/"
  }
  
  # Calculate fold changes and p-values
  contrasts_results <- eBayes(contrasts_linear_fit)
  file_out <- paste0(root, "resultados/", root2, root3, output_name, "compl_degs.rds")
  saveRDS(contrasts_results, file = file_out)
  
  # Extract and save results for each contrast
  for (e in 1:length(colnames(contrasts_results))) {
    number_comparison <- e
    file_name <- paste0(root, "resultados/", root2, root3, gsub(" - ", "_", colnames(contrasts_results)[e]), "_degs.csv")
    
    # Extract comparison results
    comparison_results <- topTable(contrasts_results, number = 33602, coef = number_comparison, sort.by = "logFC")
    
    # Select genes with logFC > 1 or logFC < -1 and adjusted p-value < 0.05
    comparison_results <- comparison_results[which(comparison_results$adj.P.Val <= 0.05 & (comparison_results$logFC > 1 | comparison_results$logFC < -1)), ]
    write.csv(comparison_results, file = file_name)
  }
}

# First Year: Afternoon vs Morning
DEGs(year = 1, contrast = 1, output_name = "ecof1_tm_")
# First Year: Extraction vs Previous
DEGs(year = 1, contrast = 2, output_name = "ecof1_exs_")
# Second Year: Afternoon vs Morning
DEGs(year = 2, contrast = 1, output_name = "ecof2_tm_")
# Second Year: Extraction vs Previous
DEGs(year = 2, contrast = 2, output_name = "ecof2_exs_")

# Gene Set Enrichment Analysis (GSEA) ---------------------------------------
degs <- function(year, con, nam) {
  # Load flowering gene info
  flow_genes <- read.table(file.path(input_path, "Flowering_time_genes_306_22dec23.csv"), header = TRUE, sep = ";")
  
  # Set path construction
  root <- output_path
  n <- c("ex_vs_ant_ex/", "t_vs_m/")
  
  if (year == 1) {
    root2 <- "ecofun1/"
  } else {
    root2 <- "ecofun2/"
  }
  
  root3 <- n[con]
  con_names <- list.files(paste0(root, "resultados/", root2, root3), pattern = ".csv")
  
  # Load the DEGs data
  contrast <- read.csv(paste0(root, "resultados/", root2, root3, con_names[nam]), header = TRUE)
  
  # Obtain the most differentially expressed (up and down) genes
  con_sig_up <- contrast[which(contrast$logFC > 1.5 & contrast$adj.P.Val < 0.01), ]
  con_sig_dw <- contrast[which(contrast$logFC < -1.5 & contrast$adj.P.Val < 0.01), ]
  
  # DEGs list
  genes_up <- con_sig_up$X
  genes_dw <- con_sig_dw$X
  
  # Manage flowering gene sets data
  pos_fgenes <- flow_genes[which(flow_genes$Effect.on.FT == "Positive"), c(6, 3)]
  neg_fgenes <- flow_genes[which(flow_genes$Effect.on.FT == "Negative"), c(6, 3)]
}

# PCA Analysis --------------------------------------------------------------
root <- output_path
year <- 1
if (year == 1) {
  root2 <- "resultados/ecofun1/"
} else {
  root2 <- "resultados/ecofun2/"
}
root3 <- paste0("ex_vs_ant_ex/ecof", year, "_exs_compl_degs.rds")
root4 <- paste0("t_vs_m/ecof", year, "_tm_compl_degs.rds")

degs1 <- readRDS(file = file.path(root, root2, root3))
degs2 <- readRDS(file = file.path(root, root2, root4))

# Collect significant genes
genes_tot <- vector()
for (e in 1:length(colnames(degs1))) {
  number_comparison <- e
  comparison_results <- topTable(degs1, number = 33602, coef = number_comparison, sort.by = "logFC")
  comparison_results <- comparison_results[which(comparison_results$adj.P.Val <= 0.01 & (comparison_results$logFC > 1.5 | comparison_results$logFC < -1.5)), ]
  genes <- rownames(comparison_results)
  genes_tot <- union(genes, genes_tot)
}

genes_tot2 <- vector()
for (e in 1:length(colnames(degs2))) {
  number_comparison <- e
  comparison_results <- topTable(degs2, number = 33602, coef = number_comparison, sort.by = "logFC")
  comparison_results <- comparison_results[which(comparison_results$adj.P.Val <= 0.01 & (comparison_results$logFC > 1.5 | comparison_results$logFC < -1.5)), ]
  genes2 <- rownames(comparison_results)
  genes_tot2 <- union(genes2, genes_tot2)
}

# Final list of significant genes
genes_degs1 <- union(genes_tot, genes_tot2)
genes_degs2 <- union(genes_tot, genes_tot2)

write.csv(genes_degs1, file = file.path(root, root2, "genes_degs_ecofun1.csv"))
write.csv(genes_degs2, file = file.path(root, root2, "genes_degs_ecofun2.csv"))

# Load and prepare data for PCA
file_norm_eco <- paste0(root, "/datos/Quantile-normalized_ecofun", year, ".txt")
file_exp_design <- paste0(root, "/datos/exp_design_e", year, ".csv")

log_n_exp <- read.table(file_norm_eco, header = TRUE, stringsAsFactors = FALSE)
geneID <- log_n_exp$geneID
log_n_exp <- as.matrix(log_n_exp[, -1])
rownames(log_n_exp) <- geneID

exp_design <- as.matrix(read.csv(file_exp_design))

linear_fit <- lmFit(log_n_exp, exp_design)
data <- as.data.frame(t(linear_fit[["coefficients"]]))
data <- data[, colnames(data) %in% genes_degs2]

# Remove genes with zero variance
var0 <- apply(data, 2, sd)
data <- data[, var0 != 0]

# Extract metadata
extraer_info <- function(df) {
  info <- str_match(rownames(df), "(ex[1-4])(vel|moj|bon|col0)(m|t)")
  data.frame(row.names = rownames(df), Muestreo = info[, 2], Ecotipo = info[, 3], Momento = info[, 4])
}

metadata <- extraer_info(data)

# Perform PCA
data_standardized <- scale(data)
pca_results <- prcomp(data_standardized, scale. = TRUE)
summary(pca_results)

# Extract PCA results and add metadata
pca_data <- as.data.frame(pca_results$x[, 1:4])
pca_data$muestreo <- metadata$Muestreo
pca_data$ecotipo <- metadata$Ecotipo
pca_data$momento <- metadata$Momento
pca_data$exp <- rownames(metadata)
pca_data$mm <- paste0(pca_data$muestreo, pca_data$momento)
pca_data$me <- paste0(pca_data$muestreo, pca_data$ecotipo)

# Plot PCA
ggplot(pca_data, aes(PC1, PC2, label = exp, color = ecotipo)) +
  geom_point() +
  geom_text(hjust = 0, vjust = 0) +  # Add labels
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  ggtitle("PCA Plot")

# Save PCA results
saveRDS(pca_results, file = file.path(root, root2, "pca/pca_res_ecof2.rds"))
write.csv(pca_data, file = file.path(root, root2, "pca/pca_res_ecof2.csv"), row.names = TRUE)
