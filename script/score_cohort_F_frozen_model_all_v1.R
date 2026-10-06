#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(glmnet)})
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "project_paths.R"))
root <- data_root
model <- readRDS(file.path(root,"frozen_model","C_frozen_model_v1.rds")); feats <- model$feature_order
man <- fread(file.path(root,"cohort_F/derived/GSE243474_sample_manifest_v1.tsv"))
man <- man[assay_inferred=="MeDIP", .(sample=geo_accession, characteristics)]
man[, group := sub("^.*cell type: ([^;]+).*$", "\\1", characteristics)]
files <- list.files(file.path(root,"cohort_F/derived/C_model_representation_v1"), pattern="counts.tsv$", full.names=TRUE)
ids <- sub("\\.counts.tsv$", "", basename(files)); X <- matrix(0,nrow=length(feats),ncol=length(files),dimnames=list(feats,ids))
for (f in files) { d<-fread(f,header=FALSE); z<-match(feats,d[[4]]); X[,sub("\\.counts.tsv$","",basename(f))] <- ifelse(is.na(z),0,as.numeric(d[[5]][z])) }
pred <- as.numeric(predict(model$model,newx=t(log2(X+1)),type="response",s="lambda.1se"))
out <- data.table(sample=colnames(X), group=man$group[match(colnames(X),man$sample)], model_score=pred, nonzero_fraction=colMeans(X>0), total_counts=colSums(X))
fwrite(out,file.path(root,"cohort_F/qc/Cohort_F_frozen_model_scores_all_v1.tsv"),sep="\t")
print(out[, .(n=.N, mean_score=mean(model_score), median_score=median(model_score)), by=group])
