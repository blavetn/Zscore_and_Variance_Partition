# load libraries
library(data.table)
library(ggplot2)
library(RColorBrewer)
library("variancePartition")
library(DESeq2)
library(edgeR)
library(robustbase) # for covMcd() needed for robust Mahalanobis
library(rrcov) # for PcaHubert()
library(DGEobj.utils)
library(ggrepel)

# Loading raw count data and create matrix
rawCount <- fread("complete_featureCount_exon_table.tsv", sep="\t", header=T) 
rawCount[, `:=`(Chr=NULL, Start=NULL, End=NULL, Strand=NULL, Length=NULL)]
names(rawCount) <- gsub("^X","R",names(rawCount))
setcolorder(rawCount, sort(names(rawCount)))

# create a matrix
matCount <- as.matrix(rawCount[,-1])
rownames(matCount)<-rawCount$Geneid

# log2 transform the matrix
log2Count <- log2(matCount + 1)

# CPM transform
cpmCount <- cpm(matCount,log = T)
logcpmCount <- log2(cpmCount + 1) # we use this one !

# load sample info
sample_info <- fread("sample_info.tsv", sep="\t",header=T)
sample_info[!Sample %like% "^R", Sample:=paste0("R",Sample)]
setorder(sample_info, "Sample")

detect_outlier <- function(matOfCount, sample_info, 
                          gender = c("M","F","MF"), 
                          timepoint = c("6h","24h","3d","7d","3m"),
                          scaling = TRUE,
                          z_cutoff = 2, 
                          maha_cutoff = 5){
  gender <- match.arg(gender)
  timepoint <- match.arg(timepoint)
  
  # select only info of interest
  si <- sample_info[sample_info$Timepoint == timepoint & sample_info$Gender == gender,]
  # select only count of interest
  mat <- matOfCount[,si$Sample]
  
  # calculate the variance for each gene
  rv <- rowVars(mat, useNames=FALSE)
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  # perform a PCA on the data in vstcounts for the selected genes
  prcomp_data <- prcomp(t((mat)[select,]),scale=scaling)
  var_expl <- round(summary(prcomp_data)$importance[2,] *100,1)
  pc1_var <- var_expl[1]
  pc2_var <- var_expl[2]

  prdt <- as.data.table(prcomp_data$x, keep.rownames = T)
  prdt <- merge(prdt, si, by.x="rn", by.y="Sample")

  first_pca <- ggplot(prdt, aes(PC1, PC2, shape=Treatment, label=rn)) + 
    geom_point(size=2) +
    theme_bw() + 
    geom_text_repel(max.overlaps = Inf, box.padding = 0.5) + 
    xlab(paste0("PC1 - ",pc1_var,"% variance explained")) +
    ylab(paste0("PC2 - ",pc2_var,"% variance explained")) +
    theme(plot.title = element_text(face="bold")) +
    theme(legend.position="bottom") +
    ggtitle(paste0("PCA - ",gender,timepoint," - top 500 genes (highest variance)",ifelse(scaling," - scaled","")))

  # compute Z-score on pca pc1 and pc2
  prdt[, z_PC1:=scale(PC1)]
  prdt[, z_PC2:=scale(PC2)]
  prdt[, z_out:=ifelse(abs(z_PC1) > z_cutoff | abs(z_PC2) > z_cutoff, TRUE, FALSE)]

  # 2. Classic Mahalanobis
  prdt[, mahalanobis_dist:=mahalanobis(prdt[,.(PC1,PC2)], center=colMeans(prdt[,.(PC1,PC2)]), cov = cov(prdt[,.(PC1,PC2)]))]
  prdt[, m_out:=ifelse(abs(mahalanobis_dist) > maha_cutoff, TRUE, FALSE)]

  # 3. Robust Mahalanobis (Minimum Covariance Determinant)
  mcd <- covMcd(prdt[,.(PC1,PC2)])
  prdt[, mahalanobis_robust:=mahalanobis(prdt[,.(PC1,PC2)], center = mcd$center, cov = mcd$cov)]
  prdt[, mr_out:=ifelse(abs(mahalanobis_robust) > maha_cutoff, TRUE, FALSE)]

  # TODO PcaHubert

  # PCA Z-score
  z_pca <- ggplot(prdt, aes(PC1, PC2, shape=Treatment, label=rn, color=z_out)) + 
    geom_point(size=2) +
    theme_bw() + 
    scale_color_manual(values = c("black","red")) +
    geom_text_repel(max.overlaps = Inf, box.padding = 0.5) + 
    xlab(paste0("PC1 - ",pc1_var,"% variance explained")) +
    ylab(paste0("PC2 - ",pc2_var,"% variance explained")) +
    theme(plot.title = element_text(face="bold")) +
    theme(legend.position="bottom") +
    ggtitle(paste0("Z-score PCA Outlier Detection - cutoff: ",z_cutoff))

  # PCA Mahalanobis
  m_pca <- ggplot(prdt, aes(PC1, PC2, shape=Treatment, label=rn, color=m_out)) + 
    geom_point(size=2) +
    theme_bw() + 
    scale_color_manual(values = c("black","red")) +
    geom_text_repel(max.overlaps = Inf, box.padding = 0.5) + 
    xlab(paste0("PC1 - ",pc1_var,"% variance explained")) +
    ylab(paste0("PC2 - ",pc2_var,"% variance explained")) +
    theme(plot.title = element_text(face="bold")) +
    theme(legend.position="bottom") +
    ggtitle(paste0("Mahalanobis PCA Outlier Detection - cutoff: ",maha_cutoff))

  # PCA Robust Mahalanobis
  mr_pca <- ggplot(prdt, aes(PC1, PC2, shape=Treatment, label=rn, color=mr_out)) + 
    geom_point(size=2) +
    theme_bw() + 
    scale_color_manual(values = c("black","red")) +
    geom_text_repel(max.overlaps = Inf, box.padding = 0.5) + 
    xlab(paste0("PC1 - ",pc1_var,"% variance explained")) +
    ylab(paste0("PC2 - ",pc2_var,"% variance explained")) +
    theme(plot.title = element_text(face="bold")) +
    theme(legend.position="bottom") +
    ggtitle(paste0("Robust Mahalanobis PCA Outlier Detection - cutoff: ",maha_cutoff))

  return(list(si,mat,prcomp_data,pc1_var,pc2_var,prdt,first_pca,z_pca,m_pca,mr_pca))
}

f6h <- detect_outlier(logcpmCount, sample_info, "F", "6h",T,2,5)
m6h <- detect_outlier(logcpmCount, sample_info, "M", "6h",T,2,5)
f24h <- detect_outlier(logcpmCount, sample_info, "F", "24h",T,2,5)
m24h <- detect_outlier(logcpmCount, sample_info, "M", "24h",T,2,5)
f3d <- detect_outlier(logcpmCount, sample_info, "F", "3d",T,2,5)
m3d <- detect_outlier(logcpmCount, sample_info, "M", "3d",T,2,5)
f7d <- detect_outlier(logcpmCount, sample_info, "F", "7d",T,2,5)
m7d <- detect_outlier(logcpmCount, sample_info, "M", "7d",T,2,5)
f3m <- detect_outlier(logcpmCount, sample_info, "F", "3m",T,2,5)
m3m <- detect_outlier(logcpmCount, sample_info, "M", "3m",T,2,5)

# plot Hippocampus
m6h[[7]] + geom_point(data=m6h[[6]], aes(PC1,PC2, color=Hippocampus),size=2) + scale_color_manual(values = brewer.pal(8,"Set1"))

# plot Z-score
f6h[[8]] 

# plot Mahalanobis
f6h[[9]]

# plot robust Mahalanobis
f6h[[10]]



