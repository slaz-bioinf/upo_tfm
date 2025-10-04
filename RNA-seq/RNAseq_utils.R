# TFM: Endometrial transcriptomic signatures of female infertility uncovered 
# through RNA-seq and covariate-adjusted differential expression analysis
# 
# Author: Sofia Lazaro
# 
# RNA-seq functions for RNAseq-pipeline.

# Import libraries
library("edgeR")
library("RColorBrewer")
library("genefilter")
library("ggrepel")
library("ggfortify")
library("ggplot2")
library("openxlsx")
library("clusterProfiler")
library("corrplot")
library("org.Hs.eg.db")
library("vegan")
library("ggVennDiagram")
library("enrichplot")


# Create color palettes
palette4 <- c(
  "#C97BA0", # soft magenta
  "#4BA3A3", # turquoise dark
  "#F4D7E8", # very light pink
  "#A2DCDC" # turquoise light
  
)

palette2 <- c(
  "#C97BA0", # soft magenta
  "#4BA3A3" # soft magenta
)




Lau_QC<-function(label,dge,data){
  # Save DGEList object to x
  x<-dge
  # Save experimental group column from DGEList object samples
  group<-x$samples$group
  
  # Transform ReadCount data into log counts
  log_Edu<-log(rawdata,2)
  
  # Calculate counts per million
  cpm <- cpm(x)
  # Calculate log counts per million
  lcpm <- cpm(x, log=TRUE)
  
  # Average library size in millions
  L <- mean(x$samples$lib.size) * 1e-6
  # Median library size in millions
  M <- median(x$samples$lib.size) * 1e-6
  # Print mean and median of library size
  c(L, M)
  
  summary(lcpm)
  table(rowSums(x$counts==0)==6)
  
  # Generate PDF QC report #
  
  pdf(paste(path,label,"_QC.pdf",sep=""),paper="A4")
  # Define a logcpm cut off for graphics
  lcpm.cutoff <- log2(10/M + 2/L)
  # Obtain number of samples
  nsamples <- ncol(x)
  # Assign colours to samples
  col <- brewer.pal(nsamples, "Paired")
  # Divide the plot in 1x2 plots
  par(mfrow=c(1,2))
  
  # Plot density of raw data
  # Plot of the first sample
  plot(density(log_Edu[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  # Add lines of log cpm cut off calculated before
  abline(v=lcpm.cutoff, lty=3)
  # Add the density plot for the remaining samples
  for (i in 2:nsamples){
    den <- density(log_Edu[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  # Legend
  legend("topright",legend=targets$Name, text.col=col, bty = "n", cex = 0.5)
  
  # Plot density of filtered log-CPM data
  lcpm <- cpm(x, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", legend=targets$Name, text.col=col, bty="n", cex = 0.5)
  
  # Apply TMM normalization of edgeR
  x <- calcNormFactors(x)
  # Print normalized sample info
  x$samples
  # Estimate common dispersion
  x <- estimateCommonDisp(x, robust=TRUE)
  # Estimate tagwise dispersion
  x <- estimateTagwiseDisp(x)
  
  # Create a copy with normalization factors reset to 1 (for comparison)
  x2 <- x
  x2$samples$norm.factors <- 1
  
  # Variable with the exp groups as factor
  col.group <- as.factor(group)
  # Assign colours per group
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  # Convert factor to character for plotting
  col.group <- as.character(col.group)
  
  # Divide the plot in 1x2 plots
  par(mfrow=c(1,2))
  # Recalculate lcpm (unnormalized)
  lcpm <- cpm(dge, log=TRUE)
  # Boxplot of unnormalized data
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Type, cex.axis=0.4)
  title(main="Unnormalised data",ylab="Log-cpm")
  
  # Recalculate lcpm (normalized)
  lcpm <- cpm(x, log=TRUE)
  # Boxplot of normalized data
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Type, cex.axis=0.4)
  title(main="Normalised data",ylab="Log-cpm")
  
  # Define only 1 plot layout
  par(mfrow=c(1,1))
  
  #Libray size as barplots
  barplot(x$samples$lib.size,names.arg = targets$Name,las=2, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))))
  
  #Corrplot using spearman correlation
  lcpm <- cpm(x, log= FALSE)
  corrplot(cor(lcpm,method="spearman"), method='number',type = 'upper')
  corrplot(cor(lcpm,method="spearman"), order='AOE')
  
  # MDS plot (multi-dimensional scaling) using sample groups for color
  col.group <- as.factor(group)
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  levels(col.group)
  z<-plotMDS(lcpm, labels=targets$Name, col=col.group, gene.selection = "pairwise", plot=F) # Compute MDS but don’t plot yet
  edge<-sd(z$x) # Spread of first dimension (for plot limits)
  plotMDS(lcpm, labels=targets$Name, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge)) # Plot MDS with adjusted x-axis
  title(main="A. MDS-PCoA Sample Names")
  
  # MDS on log2 CPM values
  lcpm <- cpm(x, log=TRUE)
  par(mfrow=c(1,1))
  col.group <- as.factor(group)
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  levels(col.group)
  z<-plotMDS(lcpm, labels=targets$Name, col=col.group, gene.selection = "pairwise", plot=F)
  edge<-sd(z$x)
  #cat(edge)
  plotMDS(lcpm, labels=targets$Name, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  
  title(main="A. MDS-PCoA log2 Sample Names")
  
  # PCA by sample
  data_pca<-as.matrix(x) # Convert counts to a matrix
  data_pca<-as.data.frame(t(data_pca)) # Transpose: samples to rows
  rownames(data_pca)<-targets$Filename # Label rows with sample filenames
  data_pca.PC = prcomp(data_pca) # Principal component analysis
  
  # Add metadata to the PCA data frame
  data_pca$Type<-targets$Type
  data_pca$Filename<-targets$Filename
  data_pca$Name<-targets$Name
  data_pca$Sex<-targets$Sex
  data_pca$Age<-targets$Age
  data_pca$VAS_Group<-targets$VAS_Group
  data_pca$TypeII<-targets$TypeII
  
  
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Type',xlim = c(-0.8,0.8)))
  
  # Heatmap of top 250 most variable genes
  rsd <- rowSds(as.matrix(x)) # Row standard deviations
  sel <- order(rsd, decreasing=TRUE)[1:250] # Select top 250 most variable genes
  samplenames<-as.character(targets$Name) # Sample names for labeling
  heatmap(na.omit(as.matrix(x[sel,])),margins=c(10,8),main="Heatmap 250 most diff entities",cexRow=0.01,cexCol=0.5,labCol=samplenames)
  
  # Hierarchical clustering
  par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  pr.hc.c<- hclust(na.omit(dist(t(cpm(x$counts,log=T)),method = "euclidean")))
  plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of Normalized samples of ", label, sep=""), labels=targets$Name_Group, cex=0.5)
  
  dev.off()
}

project<-function(){ # Create results folder
  # Path to the folder where results will be stored
  path=paste(getwd(),"/",label,"/",sep="")
  # Execute a system command to create the folder
  system(paste("mkdir -p ", path,sep=""))
  # Print on screen the message indicating where all files will be saved
  cat(paste("All files will be saved under folder:", path),"\n")
  # Return the path of the results folder
  return(path)
}


fun_biplot <- function(logcpm,targets,num_cols,title=""){
  
  # PCA
  pca <- prcomp(t(logcpm), scale. = TRUE)
  
  # Obtain coordinates for each sample/PC and add to the targets data.frame
  pca_scores <- as.data.frame(pca$x)
  pca_scores <- cbind(pca_scores, targets)
  
  # Proyecting vectors
  # Select numeric covariates
  env_vars <- targets[, num_cols, drop=FALSE]
  
  # Associate numeric covariates with PC
  envfit_model <- envfit(pca, env_vars, permutations = 0)
  
  # Obtain the vectors
  arrows <- as.data.frame(scores(envfit_model, "vectors"))
  arrows$Variable <- rownames(arrows)
  
  # PCA plot + arrows for numeric covariates
  
  # Calculate scale factor for arrows and labels in proportion to PCA coordinates
  scale_factor <- min(
    (max(pca_scores$PC1) - min(pca_scores$PC1)) / max(abs(arrows$PC1)),
    (max(pca_scores$PC2) - min(pca_scores$PC2)) / max(abs(arrows$PC2))
  ) * 0.8
  
  # Scale vectors
  arrows_scaled <- arrows
  arrows_scaled$PC1 <- arrows$PC1 * scale_factor
  arrows_scaled$PC2 <- arrows$PC2 * scale_factor
  
  # PCA plot + covariates vectors
  ggplot(pca_scores, aes(x = PC1, y = PC2, colour = group)) +
    geom_point(size = 3) +
    geom_text(aes(label = Filename), vjust = -1, size = 3) + 
    geom_segment(data = arrows_scaled,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black",
                 inherit.aes = FALSE) +
    geom_text(data = arrows_scaled,
              aes(x = PC1, y = PC2, label = Variable),
              vjust = 0, hjust = 0, size = 3.5,
              inherit.aes = FALSE) +
    scale_colour_manual(values = palette4) +
    theme_minimal() +
    labs(title = paste("Biplot", title))
}

fun_volcano <- function(topgenes,pvaluecutoff,logFCcutoff,title=""){
  
  p.value_cutoff <- pvaluecutoff
  logFC_cutoff <- logFCcutoff
  
  #  Create colunm negLogP with -log10 p-value (transform p-values in order to more significant values appear at the top)
  topgenes$negLogP <- -log10(topgenes$P.Value)
  
  # Define significant genes
  topgenes$Significant <- with(topgenes, P.Value < p.value_cutoff & abs(logFC) > logFC_cutoff)
  
  # Add gene names to a new column
  topgenes$Gene <- rownames(topgenes)
  
  # Filter top 10 genes to add the label in the plot (exclude ENSBL genes)
  sig_genes <- topgenes[topgenes$Significant==TRUE & !grepl("^ENS",topgenes$Gene, ignore.case = TRUE), ][1:15, ]
  
  # Volcano plot
  vplot <- ggplot(topgenes, aes(x = logFC, y = negLogP, color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey", "#4BA3A3")) +
    geom_text_repel(data = sig_genes, aes(label = Gene), size = 3)+
    theme_minimal() +
    xlab("log2 Fold Change") +
    ylab("-log10(P.Value)") +
    ggtitle(paste("Volcano Plot: ",title)) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "#333333") +
    geom_hline(yintercept = -log10(p.value_cutoff), linetype = "dashed", color = "#333333")
  
  print(vplot)
  
}
# topgenes <- top_genes_1
# FCcutoff <- 1.5
# pvaluecutoff <- 0.01
# annotations <- 10

fun_mapIds <- function(topgenes){
  # Add a new column "symbol" using the row names of the table
  topgenes$symbol <- rownames(topgenes)
  
  # Using the mapIds() function, we will add different annotations to the table 
  # The key for mapping will be the gene symbol (row names)
  
  # Add the ENSEMBL identifier for each gene
  topgenes$ENSEMBL <- mapIds(org.Hs.eg.db, keys=rownames(topgenes), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
  
  # Add the full gene name ("GENENAME")
  topgenes$genename <- mapIds(org.Hs.eg.db, keys=rownames(topgenes), column="GENENAME", keytype="SYMBOL", multiVals="first")
  
  # Add the Entrez Gene identifier ("ENTREZID")
  topgenes$entrezid <- mapIds(org.Hs.eg.db, keys=rownames(topgenes), column="ENTREZID", keytype="SYMBOL", multiVals="first")
  
  # Return the annotated data frame
  return(topgenes)
  
}



fun_enrichment <- function(
    topgenes,
    pvaluecutoff=0.01, # P-value (or adjusted p-value) threshold for DEG selection
    FCcutoff=1.5, # Fold-change threshold (absolute); will be converted to log2
    annotations = 15, # Number of categories to display in plots
    output_name="enrichment", # Base name for exported files
    padjusted = FALSE# If TRUE, use adj.P.Val; if FALSE, use raw P.Value
    ){ 
  
  # Map identifiers and gene metadata (adds symbol, ENSEMBL, genename, entrezid)
  topgenes <- fun_mapIds(topgenes)
  
  # Thresholds for defining DEGs (used for enrichment input)
  # Convert FC cutoff to log2  and set p-value threshold
  logfc_threshold <- log2(FCcutoff)
  pvalue_threshold <- pvaluecutoff
  
  # Build the DEG vector (Entrez IDs) considering the fixed thresholds
  if (padjusted){
    # Use adjusted p-values
    deg <- topgenes$entrezid[which(abs(topgenes$logFC) > logfc_threshold & topgenes$adj.P.Val < pvalue_threshold)]
              
  } else{
    # Use raw p-values
    deg <- topgenes$entrezid[which(abs(topgenes$logFC) > logfc_threshold & topgenes$P.Value < pvalue_threshold)]
  }
  
  # GO enrichment (three ontologies: BP, MF, CC)
  
  fea_GO_BP <- enrichGO(gene = as.character(deg), 
                        ont = "BP", # Biological Process
                        keyType = "ENTREZID",
                        OrgDb = org.Hs.eg.db, # Human annotation database
                        pAdjustMethod = "BH") # Benjamini–Hochberg correction
  
  fea_GO_MF <- enrichGO(gene = as.character(deg), 
                        ont = "MF", # Molecular Function
                        keyType = "ENTREZID",
                        OrgDb = org.Hs.eg.db, 
                        pAdjustMethod = "BH")
  
  fea_GO_CC <- enrichGO(gene = as.character(deg), 
                        ont = "CC", # Cellular Component
                        keyType = "ENTREZID",
                        OrgDb = org.Hs.eg.db, 
                        pAdjustMethod = "BH")
  
  # GO Plots 
  # If no terms are enriched, print a message
  if (nrow(as.data.frame(fea_GO_BP))==0) {
    print("No significant enriched GO BP terms were found")
  } else{
    # If there are enriched terms, print the plots
    print(dotplot(fea_GO_BP, showCategory = annotations, title="BP: Biological Process") +
            scale_fill_gradientn(colors = palette2, name = "p.adjust") +
            theme(
            plot.title   = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x  = element_text(size = 12, hjust = 1),
            axis.text.y  = element_text(size = 10, angle = 30),
            plot.margin = margin(5, 50, 5, 5, unit = "pt")))
  }
  
  if (nrow(as.data.frame(fea_GO_MF))==0) {
    print("No significant enriched GO MF terms were found")
  } else{
  print(dotplot(fea_GO_MF, showCategory = annotations, title="MF: Molecular Function") +
          scale_fill_gradientn(colors = palette2, name = "p.adjust") +
          theme(
            plot.title   = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x  = element_text(size = 12, hjust = 1),
            axis.text.y  = element_text(size = 10, angle = 30),
            plot.margin = margin(5, 50, 5, 5, unit = "pt")))
    
  }
  
  if (nrow(as.data.frame(fea_GO_CC))==0) {
    print("No significant enriched GO CC terms were found")
  } else{
  print(dotplot(fea_GO_CC, showCategory = annotations, title="CC: Cellular Component") +
          scale_fill_gradientn(colors = palette2, name = "p.adjust") +
          theme(
            plot.title   = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x  = element_text(size = 12, hjust = 1),
            axis.text.y  = element_text(size = 10, angle = 30),
            plot.margin = margin(5, 50, 5, 5, unit = "pt")))
  }
  
  # KEGG enrichment
  fea_KEGG <- enrichKEGG(gene = as.character(deg), 
                         organism = "hsa" # Homo sapiens
                         )
  
  
  # Plot KEGG (dotplot)
  if (nrow(as.data.frame(fea_KEGG))==0) {
    print("No significant enriched KEGG terms were found")
  } else{
  print(dotplot(fea_KEGG, showCategory = annotations) + 
          ggtitle("KEGG pathway enrichment") +
          scale_fill_gradientn(colors = palette2, name = "p.adjust") +
          theme(
            plot.title   = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x  = element_text(size = 12, hjust = 1),
            axis.text.y  = element_text(size = 10, angle = 30),
            plot.margin = margin(5, 50, 5, 5, unit = "pt")))

  }
  
  # Plot KEGG (barplot)
  if (nrow(as.data.frame(fea_KEGG))==0) {
    print("No significant enriched KEGG terms were found")
  } else{
    print(barplot(fea_KEGG, showCategory = annotations, x = "GeneRatio") + 
            ggtitle("KEGG pathway enrichment") +
            scale_fill_gradientn(colors = palette2, name = "p.adjust") +
            theme(
              plot.title   = element_text(size = 16, face = "bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.x  = element_text(size = 12, hjust = 1),
              axis.text.y  = element_text(size = 10, angle = 30),
              plot.margin = margin(5, 50, 5, 5, unit = "pt")))
    
  }
  
  # Export results in a single Excel file with 4 sheets
  df_BP   <- as.data.frame(fea_GO_BP)
  df_MF   <- as.data.frame(fea_GO_MF)
  df_CC   <- as.data.frame(fea_GO_CC)
  df_KEGG <- as.data.frame(fea_KEGG)
  
  # Create workbook
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "GO_BP")
  openxlsx::writeData(wb, "GO_BP", df_BP, withFilter = TRUE)

  openxlsx::addWorksheet(wb, "GO_MF")
  openxlsx::writeData(wb, "GO_MF", df_MF, withFilter = TRUE)

  openxlsx::addWorksheet(wb, "GO_CC")
  openxlsx::writeData(wb, "GO_CC", df_CC, withFilter = TRUE)

  openxlsx::addWorksheet(wb, "KEGG")
  openxlsx::writeData(wb, "KEGG", df_KEGG, withFilter = TRUE)

  openxlsx::saveWorkbook(wb, file = paste0(output_name,"-enrichment_table.xlsx"), overwrite = TRUE)

  out <- list(
    BP   = df_BP,
    MF   = df_MF,
    CC   = df_CC,
    KEGG = df_KEGG
  )
  
  # Return results as a named list of data frames
  return(out)
  
  }
  


fun_venn <- function(genelist) {
  # Create a Venn diagram from a list of gene sets.
  # genelist must be a named list where each element is a vector of genes (End, RIF, UnI)
  
  # Create the Venn diagram object using ggVennDiagram
  venndiag <- ggVennDiagram(genelist)+
    # Customize the fill colors with the custom palette
    scale_fill_gradientn(colors = palette4)
  
  # Print the Venn diagram so it is displayed in the output
  print (venndiag)
}














