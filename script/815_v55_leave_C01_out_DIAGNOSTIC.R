#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(vegan);library(permute);library(jsonlite)})
options(scipen=999);setDTthreads(2)
R <- Sys.getenv("CHILE_PROJECT_ROOT", unset=".")
D <- file.path(R,'analysis/C03_C18_correction_2026-10-28/exclude_both')
O <- file.path(D,'DIAGNOSTIC_leave_C01_out_2026-10-29')
dir.create(O,recursive=TRUE,showWarnings=FALSE)
stopifnot(!file.exists(file.path(O,'SUCCESS.json')))
EV <- c('pH','Conductivity','moisture','C','N','Feo','Alo','Mno','Sio','NH4','NO3','Pt')
CV <- setdiff(EV,c('C','N','Feo','Alo','Sio','Pt'))
ids15 <- c('C01',sprintf('C%02d',4:17));ids9 <- sprintf('C%02d',4:12)
# Fixed before any model is evaluated: no result-dependent branch selection.
write_json(list(status='fixed_before_run',primary='PCA15_fixed',secondary='PCA14_refit_report_only',
    modes=c('original_depth_matched','complete_variables'),profiles=c('gene','KO'),
    samples=ids9,seed=42,permutations=999,BH_family='gene and KO within each PCA variant and missingness mode',
    rule='PCA15_fixed/original_depth_matched/gene q < 0.05 selects A; otherwise B',
    design='Nine independent surface soil-pit libraries, SG/LC/NA n=3 each; no depth tests',
    response='Unchanged Hellinger-transformed relative gene/KO profiles, Bray-Curtis distances',
    predictors='SoilPC1 + SoilPC2; full dbRDA model; no site adjustment or term selection'),
    file.path(O,'PRESPECIFIED.json'),pretty=TRUE,auto_unbox=TRUE)
md <- as.data.frame(fread(file.path(D,'inputs/metadata.tsv'),na.strings=c('','NA')))
m15 <- md[match(ids15,md$sample_id),,drop=FALSE]
stopifnot(identical(m15$sample_id,ids15),all(m15$depth_label[match(ids9,ids15)]=='5_10cm'))
ds <- readRDS(file.path(D,'derived/Hellinger_Bray_distances.rds'))
subsetd <- function(x,ids)as.dist(as.matrix(x)[ids,ids,drop=FALSE])
base <- fread(file.path(D,'adopted/environment_dbRDA.tsv'))
retained_scores <- fread(file.path(D,'adopted/environment_PCA_scores.tsv'))
rows<-list();scores_out<-list();load_out<-list();regression<-list()
set.seed(42);perms9<-shuffleSet(9,nset=999,control=how(nperm=999))
stopifnot(nrow(perms9)==999,ncol(perms9)==9)
fwrite(as.data.table(perms9),file.path(O,'DIAGNOSTIC_permutations_n9.tsv'),sep='\t')
for(mode_i in c('original_depth_matched','complete_variables')) {
  vars <- if(mode_i=='complete_variables') CV else EV
  m <- m15
  if(mode_i=='original_depth_matched')for(v in EV)for(i in which(is.na(m[[v]]))){
    donors<-which(m$site=='NB' & m$depth_label==m$depth_label[i] & !is.na(m[[v]]))
    stopifnot(length(donors)>0,m$site[i]=='NB')
    m[[v]][i]<-mean(m[[v]][donors])
  }
  stopifnot(!anyNA(m[,vars]))
  pc15<-prcomp(m[,vars],center=TRUE,scale.=TRUE)
  rownames(pc15$x)<-ids15
  if(mode_i=='original_depth_matched'){
    expected<-as.matrix(retained_scores[match(ids15,sample_id),.(PC1,PC2)])
    stopifnot(max(abs(pc15$x[,1:2]-expected))<1e-10)
  }
  # Reproduce the existing ten-library model before adding the diagnostic.
  for(pr in c('gene','KO')){
    s<-m$sample_id[m$depth_label=='5_10cm']
    pred<-data.frame(SoilPC1=pc15$x[s,1],SoilPC2=pc15$x[s,2]);dd<-subsetd(ds[[pr]],s)
    fit<-dbrda(dd~SoilPC1+SoilPC2,data=pred);rr<-RsquareAdj(fit);set.seed(42)
    p<-anova(fit,permutations=999)$`Pr(>F)`[1]
    ref<-base[mode==mode_i & scope=='surface' & profile==pr]
    stopifnot(nrow(ref)==1,abs(rr$r.squared-ref$R2)<1e-12,abs(p-ref$p)<1e-12)
    regression[[length(regression)+1]]<-data.table(mode=mode_i,profile=pr,R2=rr$r.squared,p=p,max_delta=max(abs(c(rr$r.squared-ref$R2,p-ref$p))))
  }
  for(variant in c('PCA15_fixed','PCA14_refit')){
    pc<-pc15
    if(variant=='PCA14_refit'){
      pc<-prcomp(m[m$sample_id!='C01',vars],center=TRUE,scale.=TRUE)
      rownames(pc$x)<-m$sample_id[m$sample_id!='C01']
    }
    scores_out[[length(scores_out)+1]]<-data.table(mode=mode_i,variant=variant,sample_id=rownames(pc$x),PC1=pc$x[,1],PC2=pc$x[,2])
    load_out[[length(load_out)+1]]<-data.table(mode=mode_i,variant=variant,variable=rownames(pc$rotation),PC1=pc$rotation[,1],PC2=pc$rotation[,2])
    for(pr in c('gene','KO')){
      pred<-data.frame(SoilPC1=pc$x[ids9,1],SoilPC2=pc$x[ids9,2]);dd<-subsetd(ds[[pr]],ids9)
      stopifnot(qr(model.matrix(~SoilPC1+SoilPC2,pred))$rank==3)
      fit<-dbrda(dd~SoilPC1+SoilPC2,data=pred);rr<-RsquareAdj(fit)
      at<-anova(fit,permutations=perms9)
      rows[[length(rows)+1]]<-data.table(status='DIAGNOSTIC',variant=variant,mode=mode_i,profile=pr,
        n=9L,PCA_n=if(variant=='PCA15_fixed')15L else 14L,n_variables=length(vars),
        R2=rr$r.squared,adjusted_R2=rr$adj.r.squared,F=at$F[1],p=at$`Pr(>F)`[1],
        permutation_rows=999L,seed=42L,samples=paste(ids9,collapse=';'),
        role=if(variant=='PCA15_fixed')'primary influence diagnostic' else 'secondary report-only diagnostic')
    }
  }
}
ans<-rbindlist(rows);ans[,q:=p.adjust(p,'BH'),by=.(variant,mode)]
setorder(ans,variant,mode,profile)
ans[,source_tsv_line:=.I+1L]
fwrite(ans,file.path(O,'DIAGNOSTIC_environment_dbRDA_leave_C01_out.tsv'),sep='\t')
fwrite(rbindlist(scores_out),file.path(O,'DIAGNOSTIC_soil_PCA_scores.tsv'),sep='\t')
fwrite(rbindlist(load_out),file.path(O,'DIAGNOSTIC_soil_PCA_loadings.tsv'),sep='\t')
fwrite(rbindlist(regression),file.path(O,'QA_V54_surface_model_reproduction.tsv'),sep='\t')
decision<-ans[variant=='PCA15_fixed' & mode=='original_depth_matched' & profile=='gene']
stopifnot(nrow(decision)==1)
write_json(list(status='complete',branch=if(decision$q<.05)'A' else 'B',decision=as.list(decision),
               new_tests=nrow(ans),V54_regression_max_delta=max(rbindlist(regression)$max_delta),
               native_Word='not attempted'),file.path(O,'SUCCESS.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(O,'sessionInfo.txt'))
print(ans)
