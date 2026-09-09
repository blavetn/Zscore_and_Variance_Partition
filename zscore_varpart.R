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

Rnor6 <- "complete_featureCount_exon_table_Rnor.tsv"
mRatBN7 <- "complete_featureCount_exon_table_mRatBN.tsv"
GRCr8 <- "complete_featureCount_exon_table_GRCr.tsv"

normalize.count <- function(rawcounttab = "complete_featureCount_exon_table_Rnor.tsv",
                            transformation=c("log","cpm","logcpm","voom")){
  
  # Loading raw count data and create matrix
  rawCount <- fread(rawcounttab, sep="\t", header=T) 
  rawCount[, `:=`(Chr=NULL, Start=NULL, End=NULL, Strand=NULL, Length=NULL)]
  names(rawCount) <- gsub("^X","R",names(rawCount))
  setcolorder(rawCount, sort(names(rawCount)))

  # create a matrix
  matCount <- as.matrix(rawCount[,-1])
  rownames(matCount)<-rawCount$Geneid

  # transformation
  if(transformation == "log"){
    normCount <- log2(matCount + 1)
  }else if(transformation == "log"){
    normCount<-cpm(matCount)
  }else if(transformation == "logcpm"){
    cpmCount<-cpm(matCount)
    normCount<-log2(cpmCount + 1)
  }else{
    normCount<-limma::voom(matCount)
  }

  return(normCount)
}

rnor<-normalize.count(paste0("data/",Rnor6),"logcpm")
mrat<-normalize.count(paste0("data/",mRatBN7),"logcpm")
grcr<-normalize.count(paste0("data/",GRCr8),"logcpm")

# load sample info
sample_info <- fread("data/sample_info.tsv", sep="\t",header=T)
sample_info[!Sample %like% "^R", Sample:=paste0("R",Sample)]
setorder(sample_info, "Sample")
sample_info[, Timepoint_continuous:=as.integer(factor(Timepoint, levels=c("6h","24h","3d","7d","3m")))]
sample_info[, Novaseq_Run_continuous:=as.integer(factor(Novaseq_Run, levels=c("A","B","C","D","E","F","G","I","J","K","L","M","N","O","P","Q")))]


detect_outlier <- function(matOfCount, sample_info, 
                          gender = c("M","F","MF"), 
                          timepoint = c("6h","24h","3d","7d","3m"),
                          scaling = TRUE,
                          z_cutoff = 2, 
                          maha_cutoff = 5){
  gender <- match.arg(gender)
  timepoint <- match.arg(timepoint)
  
  # select only info of interest
  if(gender != "MF"){
    si <- sample_info[sample_info$Timepoint == timepoint & sample_info$Gender == gender,]
  }else{
    si <- sample_info[sample_info$Timepoint == timepoint,]
  }
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

  si<-as.data.frame(si)
  rownames(si)<-si$Sample
  return(list(si,mat,prcomp_data,pc1_var,pc2_var,prdt,first_pca,z_pca,m_pca,mr_pca))
}

# Rnor
rnor_f6h <- detect_outlier(rnor, sample_info, "F", "6h",T,2,5)
rnor_m6h <- detect_outlier(rnor, sample_info, "M", "6h",T,2,5)
rnor_f24h <- detect_outlier(rnor, sample_info, "F", "24h",T,2,5)
rnor_m24h <- detect_outlier(rnor, sample_info, "M", "24h",T,2,5)
rnor_f3d <- detect_outlier(rnor, sample_info, "F", "3d",T,2,5)
rnor_m3d <- detect_outlier(rnor, sample_info, "M", "3d",T,2,5)
rnor_f7d <- detect_outlier(rnor, sample_info, "F", "7d",T,2,5)
rnor_m7d <- detect_outlier(rnor, sample_info, "M", "7d",T,2,5)
rnor_f3m <- detect_outlier(rnor, sample_info, "F", "3m",T,2,5)
rnor_m3m <- detect_outlier(rnor, sample_info, "M", "3m",T,2,5)

rnor_mf6h <- detect_outlier(rnor, sample_info, "MF", "6h",T,2,5)

# mRat
mrat_f6h <- detect_outlier(mrat, sample_info, "F", "6h",T,2,5)
mrat_m6h <- detect_outlier(mrat, sample_info, "M", "6h",T,2,5)
mrat_f24h <- detect_outlier(mrat, sample_info, "F", "24h",T,2,5)
mrat_m24h <- detect_outlier(mrat, sample_info, "M", "24h",T,2,5)
mrat_f3d <- detect_outlier(mrat, sample_info, "F", "3d",T,2,5)
mrat_m3d <- detect_outlier(mrat, sample_info, "M", "3d",T,2,5)
mrat_f7d <- detect_outlier(mrat, sample_info, "F", "7d",T,2,5)
mrat_m7d <- detect_outlier(mrat, sample_info, "M", "7d",T,2,5)
mrat_f3m <- detect_outlier(mrat, sample_info, "F", "3m",T,2,5)
mrat_m3m <- detect_outlier(mrat, sample_info, "M", "3m",T,2,5)

# GRCr8
grcr_f6h <- detect_outlier(grcr, sample_info, "F", "6h",T,2,5)
grcr_m6h <- detect_outlier(grcr, sample_info, "M", "6h",T,2,5)
grcr_f24h <- detect_outlier(grcr, sample_info, "F", "24h",T,2,5)
grcr_m24h <- detect_outlier(grcr, sample_info, "M", "24h",T,2,5)
grcr_f3d <- detect_outlier(grcr, sample_info, "F", "3d",T,2,5)
grcr_m3d <- detect_outlier(grcr, sample_info, "M", "3d",T,2,5)
grcr_f7d <- detect_outlier(grcr, sample_info, "F", "7d",T,2,5)
grcr_m7d <- detect_outlier(grcr, sample_info, "M", "7d",T,2,5)
grcr_f3m <- detect_outlier(grcr, sample_info, "F", "3m",T,2,5)
grcr_m3m <- detect_outlier(grcr, sample_info, "M", "3m",T,2,5)

grcr_mf6h <- detect_outlier(rnor, sample_info, "MF", "6h",T,2,8)
grcr_mf24h <- detect_outlier(rnor, sample_info, "MF", "24h",T,2,8)
grcr_mf3d <- detect_outlier(rnor, sample_info, "MF", "3d",T,2,8)
grcr_mf7d <- detect_outlier(rnor, sample_info, "MF", "7d",T,2,8)
grcr_mf3m <- detect_outlier(rnor, sample_info, "MF", "3m",T,2,8)

# plot Hippocampus
rnor_m6h[[7]] + geom_point(data=rnor_m6h[[6]], aes(PC1,PC2, color=Hippocampus),size=2) + scale_color_manual(values = brewer.pal(8,"Set1"))

# plot Z-score
rnor_f6h[[8]] 

# plot Mahalanobis
rnor_f6h[[9]]

# plot robust Mahalanobis
rnor_f6h[[10]]

## Variance Partinioning
df.si <- as.data.frame(sample_info)
rownames(df.si)<-df.si$Sample

# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical,
# so model them as random effects
# Note the syntax used to specify random effects

# Use the data as they are: 
formul <- ~ (1 | Litter) + (1 | Hippocampus) +
            (1 | Timepoint) + (1 | Gender) +
            (1 | Treatment) + Library_Prep_batch +
            (1 | Novaseq_Run)

# alternate 1 - change Library_Prep_batch as factor
formulF <- ~ (1 | Litter) + (1 | Hippocampus) +
            (1 | Timepoint) + (1 | Gender) +
            (1 | Treatment) + (1 | Library_Prep_batch) +
            (1 | Novaseq_Run)

# alternate 2 - change timepoint and novaseq run as continuous variables
formulC <- ~ (1 | Litter) + (1 | Hippocampus) +
            Timepoint_continuous + (1 | Gender) +
            (1 | Treatment) + Library_Prep_batch +
            Novaseq_Run_continuous

# reduce to be used on comparison level ex: TF_6h vs CF_6h
formul_comp <- ~ (1 | Hippocampus) +
            (1 | Treatment) + Library_Prep_batch 

# reduce to be used on all sample from same timepoint ex: TF_6h, CF_6h, TM_6h, CM_6h
formul_comp_sex <- ~ (1 | Gender) + (1 | Hippocampus) +
            (1 | Treatment) + Library_Prep_batch +
            (1 | Litter)

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified,
#     a linear mixed model is used
# If all variables are modeled as fixed effects,
#       a linear model is used
# each entry in results is a regression model fit on a single gene
# 2) extract variance fractions from each model fit
# for each gene, returns fraction of variation attributable
#       to each variable
# Interpretation: the variance explained by each variables
# after correcting for all other variables
# Note that geneExpr can either be a matrix,
# and EList output by voom() in the limma package,
# or an ExpressionSet
# log2 transformed
varPart_rnor <- fitExtractVarPartModel(rnor, formul, df.si)
varPart_mrat <- fitExtractVarPartModel(mrat, formul, df.si)
varPart_grcr <- fitExtractVarPartModel(grcr, formul, df.si)
varPartC_rnor <- fitExtractVarPartModel(rnor, formulC, df.si)
varPartC_mrat <- fitExtractVarPartModel(mrat, formulC, df.si)
varPartC_grcr <- fitExtractVarPartModel(grcr, formulC, df.si)
varPartF_rnor <- fitExtractVarPartModel(rnor, formulF, df.si2)
varPartF_mrat <- fitExtractVarPartModel(mrat, formulF, df.si2)
varPartF_grcr <- fitExtractVarPartModel(grcr, formulF, df.si2)


varPart_rnor_mf6h <- fitExtractVarPartModel(rnor_mf6h[[2]], formul_comp_sex, as.data.frame(rnor_mf6h[[1]]))


# violin plot of contribution of each variable to total variance
plotVarPart(varPart_rnor)
plotVarPart(varPart_mrat)
plotVarPart(varPart_grcr)
plotVarPart(varPartC_rnor)
plotVarPart(varPartC_mrat)
plotVarPart(varPartC_grcr)
plotVarPart(varPartF_rnor)
plotVarPart(varPartF_mrat)
plotVarPart(varPartF_grcr)

# save result in RDS
saveRDS(varPart_rnor, "results/varPart_rnor.RDS")
saveRDS(varPart_mrat, "results/varPart_mrat.RDS")
saveRDS(varPart_grcr, "results/varPart_grcr.RDS")
saveRDS(varPartC_rnor, "results/varPartC_rnor.RDS")
saveRDS(varPartC_mrat, "results/varPartC_mrat.RDS")
saveRDS(varPartC_grcr, "results/varPartC_grcr.RDS")
saveRDS(varPartF_rnor, "results/varPartF_rnor.RDS")
saveRDS(varPartF_mrat, "results/varPartF_mrat.RDS")
saveRDS(varPartF_grcr, "results/varPartF_grcr.RDS")

save.image("results/results.RData")
