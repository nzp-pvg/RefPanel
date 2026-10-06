#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "project_paths.R"))
d<-read.delim(file.path(data_root,'cohort_F','qc','Cohort_F_frozen_model_scores_all_v1.tsv'))
m<-read.delim(file.path(data_root,'cohort_F','derived','GSE243474_locked_manifest_v1.tsv'))
d<-merge(d,m[,c('geo_accession','primary_eligible')],by.x='sample',by.y='geo_accession'); d<-d[d$primary_eligible=='YES',]
auc<-function(y,s){r<-rank(s);(sum(r[y==1])-sum(seq_len(sum(y))))/(sum(y)*sum(!y))}
set.seed(1); B<-2000; out<-list(); k<-1
for(g in setdiff(unique(d$group),'Healthy')){x<-d[d$group%in%c(g,'Healthy'),]; y<-x$group==g; est<-auc(y,x$model_score); z<-rep(NA_real_,B); for(b in seq_len(B)){ii<-c(sample(which(y),sum(y),TRUE),sample(which(!y),sum(!y),TRUE)); z[b]<-auc(y[ii],x$model_score[ii])}; q<-quantile(z,c(.025,.975),na.rm=TRUE); out[[k]]<-data.frame(group=g,n_total=nrow(x),n_case=sum(y),n_healthy=sum(!y),auc=est,ci_low=q[1],ci_high=q[2]); k<-k+1}
out<-do.call(rbind,out); write.table(out,file.path(data_root,'cohort_F','qc','Cohort_F_frozen_deployment_bootstrap_v1.tsv'),sep='\t',row.names=FALSE,quote=FALSE); print(out)
