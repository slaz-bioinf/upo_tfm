# TFM: Endometrial transcriptomic signatures of female infertility uncovered 
# through RNA-seq and covariate-adjusted differential expression analysis
# 
# Author: Sofia Lazaro
# 
# RNA-seq pipeline using limma-voom and selected covariates

# Load packages
library(edgeR)
library(limma)
library(openxlsx)

# Set working directory
setwd("D:/Escritorio/AnalisisR")
source("D:/Escritorio/AnalisisR/RNAseq_utils.R")

# Define files
rna_seq_counts <- "D:/Escritorio/AnalisisR/str-ReadCount.tab"
samples <- "D:/Escritorio/AnalisisR/targets.txt"
samples_cov <- "D:/Escritorio/AnalisisR/targets_cov.csv"

##################### Preliminary analysis of RNA-seq counts ##################### 

# Folder to store QC results
label <- "QC_outliers_final"
path <- project()

# Read RNA-seq counts
rawdata<-read.delim(file=rna_seq_counts,header=TRUE,row.names=1,check.names = F)
colnames(rawdata)<-gsub("_nat","",colnames(rawdata))
rawdata<-rawdata[,sort(colnames(rawdata))]

# Read sample names and experimental group classifications
targets<-read.table(file = samples,header = T,stringsAsFactors=F)
targets
colnames(rawdata)
targets$Filename  
targets$Name

table(targets$Type)
targets$Filename<-as.vector(targets$Filename)

# Join sample name and group for plotting 
targets$Name_Group <- paste0(targets$Filename,"_",targets$Type)

# Order sample names in the raw data table to match those in the targets table
rawdata<-rawdata[,targets$Filename]

# Create DGElist object
dge <- DGEList(counts = rawdata, group = targets$Type)
dim(dge)

# Filter low expression genes with default parameters
keep <- filterByExpr(dge, group = targets$Type)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dim(dge)

# Apply TMM normalization method to dge object, calcNormFactors only calculates factors to be used downstream
dge <- calcNormFactors(dge)

# Quality analyisis
Lau_QC(label=label,dge=dge,data=rawdata)



##################### Limma-voom analysis #####################

# Label of the project
cov_adj <- TRUE
label <- "RNA-seq_results_final"
if (cov_adj){
  label <- paste0(label, "-cov")
} else {
  label <- paste0(label, "-no-cov-adj")
}

# Create fold to store the results
path <- project()

# Read RNA seq counts
rawdata<-read.delim(file=rna_seq_counts,header=TRUE,row.names=1,check.names = F)
colnames(rawdata)<-gsub("_nat","",colnames(rawdata))
rawdata<-rawdata[,sort(colnames(rawdata))]

# Read metadata table with the final selection of samples (36) and covariates
targets<-read.table(file = samples_cov, header = T, stringsAsFactors=F, sep = ",")
targets
colnames(rawdata)
targets$Filename
targets$Name


# Transform character columns to factor
for (col in names(targets)) {
  if (is.character(targets[[col]])) {
    targets[[col]] <- as.factor(targets[[col]])
  }
}

# Convert Filename column to character again
targets$Filename<-as.vector(targets$Filename)

# Order samples name in rawdata table
rawdata<-rawdata[,targets$Filename]

# Create new columns with the translated names of selected covariates for biplot representation
targets["Childhood passive smoking exposure"] <- targets["Ambiente.fumador.infancia"] 
targets["Daily sitting time"] <- targets["Horas.sentada.dia.cod"]
targets["Smoke status"] <- targets["Fumadora.si.no"]


#################### GROUP 4 #################### 
# Define the name to be used for the results files
if (cov_adj){
  analysis <- "Cov-4-class"
} else {
  analysis <- "4-class"
}
name <- file.path(label,analysis)

# Define the selected covariates
covariates <- c("BMI", "Childhood passive smoking exposure", "Daily sitting time")

# Create DGElist object
dge <- DGEList(counts = rawdata, group = targets$Type)
dim(dge)

# Create the experimental design
group <- factor(targets$Type)

# Create a design matrix with samples as rows and experimental groups as columns,
# indicating membership with 0/1 for each sample
if (cov_adj){
  design <- model.matrix(~0 + group + BMI + Ambiente.fumador.infancia + Horas.sentada.dia.cod, data = targets)
} else {
  design <- model.matrix(~0 + group)
}
colnames(design)

# Filter low expression genes with default parameters
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dim(dge)

# Apply TMM normalization method to dge object, calcNormFactors only calculates factors to be used downstream
dge <- calcNormFactors(dge)

# Voom transforms the data to be used with limma, plot TRUE to show the graph
v <- voom(dge, design, plot = TRUE)

#  Create a linear model for each gene
fit <- lmFit(v, design)


# BIPLOT

# Obtain numeric columns to be projected in the PCA plot
biplot_cols <- covariates
biplot_bool_full <- colnames(targets) %in% biplot_cols

pdf(paste0(name,"-biplot.pdf"), paper = "A4r")
fun_biplot(v$E,targets,biplot_bool_full,analysis)
dev.off()

# Define Cutoffs for DEGs
p.value_cutoff <- 0.01
logFC_cutoff <- 0.58 #1 for 2x expression and 0.58 for 1.5X expression

# 1. Comparison Unexplained_infertility vs. Male_Factor
# Create contrast matrix
contr1 <- makeContrasts(UnI.vs.MF = groupUnexplained_infertility - groupMale_factor, levels = design)
contr1

# Apply contrasts to the linear model
fit1 <- contrasts.fit(fit, contr1)
fit1 <- eBayes(fit1)

# Extract all genes sorted by p-value
top_genes_1 <- topTable(fit1, coef = "UnI.vs.MF", number=Inf, sort.by="P")
top_genes_1

# Filter DEGs with logFC and pvalue cutoffs
deg_sig1 <- top_genes_1[abs(top_genes_1$logFC) > logFC_cutoff & top_genes_1$P.Value < p.value_cutoff, ]

# Define files name
name <- file.path(label, paste0(analysis,"-", colnames(contr1)))

# If there are significant DEGs
if (nrow(deg_sig1)>0){
  # Add gene ids to deg table
  deg_sig1 <- fun_mapIds(deg_sig1)

  # Export whole DEG table (optional) and significant DEGs to txt
  # write.xlsx(top_genes_3, file = paste0(name,"-DEGs.xlsx"), row.names = TRUE)
  write.xlsx(deg_sig1, file = paste0(name,"-DEGs-sig.xlsx"), row.names = TRUE)

  #Volcano
  pdf(paste0(name,"-volcano.pdf"), paper = "A4r")
  fun_volcano(top_genes_1, p.value_cutoff,logFC_cutoff,title=colnames(contr1))
  dev.off()

  #Enrichment
  pdf(paste0(name,"-Enrichment.pdf"), paper = "A4r")
  enr1 <- fun_enrichment(top_genes_1, FCcutoff= 1.5, annotations = 10, output_name = name)
  dev.off()

}else{
  print(paste("No significant DEGs with p-value <", p.value_cutoff, 
               "and |log2FC| >", logFC_cutoff, 
               "were found in comparison", colnames(contr1)))
}

# 2. Comparison Endometriosis vs. Male_Factor
# Create contrast matrix
contr2 <- makeContrasts(End.vs.MF = groupEndometriosis - groupMale_factor, levels = design)
contr2

# Apply contrasts to the linear model
fit2 <- contrasts.fit(fit, contr2)
fit2 <- eBayes(fit2)

# Extract all genes sorted by p-value
top_genes_2 <- topTable(fit2, coef = "End.vs.MF", number=Inf, sort.by="P")
top_genes_2

# Filter DEGs with logFC and pvalue cutoffs
deg_sig2 <- top_genes_2[abs(top_genes_2$logFC) > logFC_cutoff & top_genes_2$P.Value < p.value_cutoff, ]

# Define files name
name <- file.path(label, paste0(analysis,"-", colnames(contr2)))

# If there are significant DEGs
if (nrow(deg_sig2)>0){

  # Add gene ids to deg table
  deg_sig2 <- fun_mapIds(deg_sig2)

  # Export whole DEG table (optional) and significant DEGs to txt
  # write.xlsx(top_genes_2, file = paste0(name,"-DEGs.xlsx"), row.names = TRUE)
  write.xlsx(deg_sig2, file = paste0(name,"-DEGs-sig.xlsx"), row.names = TRUE)

  #Volcano
  pdf(paste0(name,"-volcano.pdf"), paper = "A4r")
  fun_volcano(top_genes_2, p.value_cutoff,logFC_cutoff,title=colnames(contr2))
  dev.off()
  
  #Enrichment
  pdf(paste0(name,"-Enrichment.pdf"), paper = "A4r")
  enr2 <- fun_enrichment(top_genes_2, FCcutoff= 1.5, annotations = 10, output_name = name)
  dev.off()
  
} else{
  print(paste("No significant DEGs with p-value <", p.value_cutoff, 
              "and |log2FC| >", logFC_cutoff, 
              "were found in comparison", colnames(contr2)))
}

# 3. Comparison RIF vs. Male_Factor
# Create contrast matrix
contr3 <- makeContrasts(RIF.vs.MF = groupRIF - groupMale_factor, levels = design)
contr3

# Apply contrasts to the linear model
fit3 <- contrasts.fit(fit, contr3)
fit3 <- eBayes(fit3)

# Extract all genes sorted by p-value
top_genes_3 <- topTable(fit3, coef = "RIF.vs.MF", number=Inf, sort.by="P")
top_genes_3

# Filter DEGs with logFC and pvalue cutoffs
deg_sig3 <- top_genes_3[abs(top_genes_3$logFC) > logFC_cutoff & top_genes_3$P.Value < p.value_cutoff, ]

# Define files name
name <- file.path(label, paste0(analysis,"-", colnames(contr3)))

# If there are significant DEGs
if (nrow(deg_sig3)>0){
  # Add gene ids to deg table
  deg_sig3 <- fun_mapIds(deg_sig3)
  
  # Export whole DEG table (optional) and significant DEGs to txt
  # write.xlsx(top_genes_3, file = paste0(name,"-DEGs.xlsx"), row.names = TRUE)
  write.xlsx(deg_sig3, file = paste0(name,"-DEGs-sig.xlsx"), row.names = TRUE)
  
  #Volcano
  pdf(paste0(name,"-volcano.pdf"), paper = "A4r")
  fun_volcano(top_genes_3, p.value_cutoff,logFC_cutoff,title=colnames(contr3))
  dev.off()
  
  #Enrichment
  
  pdf(paste0(name,"-Enrichment.pdf"), paper = "A4r")
  enr3 <- fun_enrichment(top_genes_3, FCcutoff= 1.5, annotations = 10, output_name = name)
  dev.off()
} else{
  print(paste("No significant DEGs with p-value <", p.value_cutoff, 
              "and |log2FC| >", logFC_cutoff, 
              "were found in comparison", colnames(contr3)))
}


#################### GROUP 2 #################### 
# Define the name to be used for the results files
if (cov_adj){
  analysis <- "Cov-2-class"
} else {
  analysis <- "2-class"
}
name <- file.path(label,analysis)

#Reorder levels to match color group
targets$DiseaseClass <- factor(
  targets$DiseaseClass,
  levels = c("Pathology", "NoPathology")
)

# Define the selected covariates
covariates <- c("Smoke status", "Childhood passive smoking exposure")

# Create DGElist object
dge <- DGEList(counts = rawdata, group = targets$DiseaseClass)
dim(dge)

# Create the experimental design
group <- factor(targets$DiseaseClass)

# Create a design matrix with samples as rows and experimental groups as columns,
# indicating membership with 0/1 for each sample
if (cov_adj){
design <- model.matrix(~0 + group + Fumadora.si.no + Ambiente.fumador.infancia, data = targets)
} else{
  design <- model.matrix(~0 + group)
}
colnames(design)

# Filter low expression genes with default parameters
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dim(dge)

# Apply TMM normalization method to dge object, calcNormFactors only calculates factors to be used downstream
dge <- calcNormFactors(dge)

# Voom transforms the data to be used with limma, plot TRUE to show the graph
v <- voom(dge, design, plot = TRUE)

#  Create a linear model for each gene
fit <- lmFit(v, design)


# Biplot
# Obtain numeric columns to be projected in the PCA plot
biplot_cols <- covariates
biplot_bool_full <- colnames(targets) %in% biplot_cols

pdf(paste0(name,"-biplot.pdf"), paper = "A4r")
fun_biplot(v$E,targets,biplot_bool_full,analysis)
dev.off()

# Define Cutoffs for DEGs
p.value_cutoff <- 0.01
logFC_cutoff <- 0.58 #1 for 2x expression and 0.58 for 1.5X expression

# 7. Comparison Pathology vs. NoPathology
# Create contrast matrix
contr7 <- makeContrasts(Pat.vs.NoP = groupPathology - groupNoPathology, levels = design)
contr7

# Apply contrasts to the linear model
fit7 <- contrasts.fit(fit, contr7)
fit7 <- eBayes(fit7)

# Extract all genes sorted by p-value
top_genes_7 <- topTable(fit7, coef = "Pat.vs.NoP", number=Inf, sort.by="P")
top_genes_7

# Define files name
name <- file.path(label, paste0(analysis,"-", colnames(contr7)))

# Filter DEGs with logFC and pvalue cutoffs
deg_sig7 <- top_genes_7[abs(top_genes_7$logFC) > logFC_cutoff & top_genes_7$P.Value < p.value_cutoff, ]

# If there are significant DEGs
if (nrow(deg_sig7)>0){
  # Add gene ids to deg table
  deg_sig7 <- fun_mapIds(deg_sig7)
  
  # Export whole DEG table (optional) and significant DEGs to txt
  # write.xlsx(top_genes_7, file = paste0(name,"-DEGs.xlsx"), row.names = TRUE)
  write.xlsx(deg_sig7, file = paste0(name,"-DEGs-sig.xlsx"), row.names = TRUE)
  
  #Volcano
  pdf(paste0(name,"-volcano.pdf"), paper = "A4r")
  fun_volcano(top_genes_7, p.value_cutoff,logFC_cutoff,title=colnames(contr7))
  dev.off()
  
  #Enrichment p-value 0.01
  pdf(paste0(name,"-Enrichment.pdf"), paper = "A4r")
  enr7 <- fun_enrichment(top_genes_7, FCcutoff= 1.5, annotations = 10, output_name = name, padjusted = FALSE)
  dev.off()
  
  
} else{
  print(paste("No significant DEGs with p-value <", p.value_cutoff, 
              "and |log2FC| >", logFC_cutoff, 
              "were found in comparison", colnames(contr7)))
}


################## VENN Diagrams ##################
# Create a list of significant DEGs (symbols) for each comparison
# - End: Endometriosis vs Male factor
# - RIF: Recurrent implantation failure vs Male factor
# - UnI: Unexplained infertility vs Male factor
list_4 <- list(End = deg_sig2$symbol, RIF = deg_sig3$symbol, UnI = deg_sig1$symbol)

# Print the Venn diagram
fun_venn(list_4)

# For each group, combine the enriched GO IDs from:
# - Biological Process (BP)
# - Molecular Function (MF)
# - Cellular Component (CC)
list_func_3_GO <- list(End_GO = c(enr2$BP$ID, enr2$MF$ID, enr2$CC$ID), 
                       RIF_GO=c(enr3$BP$ID, enr3$MF$ID, enr3$CC$ID), 
                       UnI_GO=c(enr1$BP$ID, enr1$MF$ID, enr1$CC$ID))

# For each group, extract the enriched KEGG pathway IDs
list_func_3_KEGG <- list(End_KEGG = c(enr2$KEGG$ID, enr2$KEGG$ID, enr2$KEGG$ID), 
                         RIF_KEGG=c(enr3$KEGG$ID, enr3$KEGG$ID, enr3$KEGG$ID), 
                         UnI_KEGG=c(enr1$KEGG$ID, enr1$KEGG$ID, enr1$KEGG$ID))

# Export Venn diagram of DEGs to PDF
pdf(file.path(label, "Venn.pdf"), paper = "A4r")
fun_venn(list_4)
dev.off()

# Export Venn diagram of GO terms to PDF
pdf(file.path(label, "Venn-GO.pdf"), paper = "A4r")
fun_venn(list_func_3_GO)
dev.off()

