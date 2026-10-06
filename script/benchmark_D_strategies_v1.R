#!/usr/bin/env Rscript
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(glmnet) })
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "project_paths.R"))
root <- file.path(data_root, "cohort_D", "derived_v1")
auc_rank <- function(prob, y) { n1<-sum(y==1); n0<-sum(y==0); r<-rank(prob); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
run <- function(file, material) {
  material_name <- material
  x <- readRDS(file); ids <- colnames(x); y <- as.integer(grepl("patient", ids)); X <- t(log2(x+1))
  conventional <- rowMeans(X)
  set.seed(11); foldid <- sample(rep(1:5,length.out=nrow(X)))
  adaptive <- rep(NA_real_, nrow(X))
  for(k in 1:5) { tr<-foldid!=k; te<-foldid==k; fit<-cv.glmnet(X[tr,,drop=FALSE],y[tr],family="binomial",alpha=1,nfolds=5,type.measure="deviance",grouped=FALSE); adaptive[te]<-as.numeric(predict(fit,newx=X[te,,drop=FALSE],s="lambda.1se",type="response")) }
  locked <- fread(file.path(data_root, "deployment", "D_frozen_scores_v1.tsv"))[material == material_name]
  out <- rbind(data.table(material=material,strategy="conventional",n=nrow(X),auc=auc_rank(conventional,y),uses_D_labels=FALSE),
               data.table(material=material,strategy="adaptive_5fold",n=nrow(X),auc=auc_rank(adaptive,y),uses_D_labels=TRUE),
               data.table(material=material,strategy="locked_A_model",n=nrow(X),auc=auc_rank(locked$frozen_score,y),uses_D_labels=FALSE))
  out
}
res <- rbind(run(file.path(root,"D_cfDNA_model20k_v1.rds"),"cfDNA"),run(file.path(root,"D_PBL_model20k_v1.rds"),"PBL"))
fwrite(res,file.path(data_root, "deployment", "D_strategy_benchmark_v1.tsv"),sep="\t")
flags <- res[, .(locked_auc=auc[strategy=="locked_A_model"], adaptive_auc=auc[strategy=="adaptive_5fold"],
                 conventional_auc=auc[strategy=="conventional"]), by=material]
flags[, `:=`(locked_failure = locked_auc < 0.60, adaptation_dependence = adaptive_auc - locked_auc > 0.10,
             reportable_locked_signal = locked_auc >= 0.60)]
fwrite(flags,file.path(data_root, "deployment", "D_deployment_failure_flags_v1.tsv"),sep="\t")
message("[DONE] benchmark and failure flags written")
