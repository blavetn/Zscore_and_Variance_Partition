# Variance partitioning analysis
library(data.table)
library(ggplot2)
library("variancePartition")
library(DESeq2)
library(edgeR)

# Loading raw count data and create matrix
rawCount <- fread("complete_featureCount_exon_table.tsv", sep="\t", header=T) 
rawCount[, `:=`(Chr=NULL, Start=NULL, End=NULL, Strand=NULL, Length=NULL)]
names(rawCount) <- gsub("^X","R",names(rawCount))
setcolorder(rawCount, sort(names(rawCount)))

# create a matrix
matCount <- as.matrix(rawCount[,-1])
rownames(matCount)<-rawCount$Geneid

# transformation
log2Count <- log2(matCount + 1)
cpmCount<-cpm(matCount)
logcpmCount<-log2(cpmCount + 1)
voomCount<-limma::voom(matCount)


# load sample info
sample_info <- fread("sample_info.tsv", sep="\t",header=T)
sample_info[!Sample %like% "^R", Sample:=paste0("R",Sample)]
# sample_info<-sample_info[Sample %in% names(rawCount),]
# Batch labels are read as numeric but are categorical: convert to factors
# sample_info$Library_Prep_batch <- factor(sample_info$Library_Prep_batch)
sample_info[, continousTimepoint:=as.integer(factor(Timepoint, levels=c("6h","24h","3d","7d","3m")))]
sample_info[, NS_run:=as.integer(factor(Novaseq_Run, levels=c("A","B","C","D","E","F","G","I","J","K","L","M","N","O","P","Q")))]

sample_info <- as.data.frame(sample_info)
rownames(sample_info)<-sample_info$Sample
sample_info <- sample_info[order(sample_info$Sample),]

# load DESeq2 table to get normalized count
ds.files <- list.files(path = "DESeq2/",
                      pattern = "DESeq2.tsv",
                      recursive = T,
                      full.names = T) 

dt <- fread(ds.files[1])
dtt <- dt[, c("Ensembl_Id", grep("normCounts", names(dt), value = TRUE)), with = FALSE] # keep gene and normCounts
names(dtt) <- gsub("_normCounts","",names(dtt)) # remove normCounts from sample name
names(dtt) <- gsub("Ensembl_Id","Geneid",names(dtt))

for(f in 2:length(ds.files)){
  dt <- fread(ds.files[f])
  dt_sub <- dt[, c("Ensembl_Id", grep("normCounts", names(dt), value = TRUE)), with = FALSE] # keep gene and normCounts
  names(dt_sub) <- gsub("_normCounts","",names(dt_sub)) # remove normCounts from sample name
  names(dt_sub) <- gsub("Ensembl_Id","Geneid",names(dt_sub))
  dtt<-merge(dtt,dt_sub,by="Geneid")
}

matt <- as.matrix(dtt[,-1])
rownames(matt) <- dtt$Geneid

# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical,
# so model them as random effects
# Note the syntax used to specify random effects
form <- ~ (1 | Treatment) + (1 | Gender) + (1 | Hippocampus) + (1 | Timepoint) + (1 | Litter) + Library_Prep_batch + (1 | Novaseq_Run)
formR <- ~ (1 | Treatment) + (1 | Gender) + (1 | Hippocampus) + (1 | Timepoint) + (1 | Litter) + Run + Library_Prep_batch

formC <- ~ NS_run + continousTimepoint + (1 | Treatment) + (1 | Gender) + (1 | Hippocampus) + (1 | Litter) + Library_Prep_batch 
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
varPart_cpm <- fitExtractVarPartModel(cpmCount, form, sample_info)
# varPartDS <- fitExtractVarPartModel(matt, form, sample_info)

# log2 transformed
varPart_log <- fitExtractVarPartModel(log2Count, form, sample_info)
varPartC_log <- fitExtractVarPartModel(log2Count, formC, sample_info)
varPartR_log <- fitExtractVarPartModel(log2Count, formR, sample_info)
# violin plot of contribution of each variable to total variance
plotVarPart(varPart_log)
plotVarPart(varPartC_log)
plotVarPart(varPartR_log)
# save result in RDS
saveRDS(varPart_log, "results/varPart_log.RDS")
saveRDS(varPartC_log, "results/varPartC_log.RDS")
saveRDS(varPartR_log, "results/varPartR_log.RDS")

varPart_logcpm <- fitExtractVarPartModel(logcpmCount, form, sample_info)
# varPartVoom <- fitExtractVarPartModel(voomCount, form, sample_info)

varPart_C <- fitExtractVarPartModel(cpmCount, formC, sample_info)
# varPartDS_C <- fitExtractVarPartModel(matt, formC, sample_info)

# violin plot of contribution of each variable to total variance
plotVarPart(varPart)
# plotVarPart(varPartDS)
plotVarPart(varPart_C)
# plotVarPart(varPartDS_C)
compute_varpart<-function(exprObj, formula, data, gender, timepoint){
  if(gender == "MF"){
    df <- data[data$Timepoint == timepoint,]
  }else{
    df <- data[data$Timepoint == timepoint & data$Gender == gender,]
  }
  mat <- exprObj[,df$Sample]
  varPart <- fitExtractVarPartModel(mat, formula, df)
  return(varPart)
}
