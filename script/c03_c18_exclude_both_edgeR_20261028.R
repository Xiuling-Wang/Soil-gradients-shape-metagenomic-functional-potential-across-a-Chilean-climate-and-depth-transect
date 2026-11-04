#!/usr/bin/env Rscript
# Preserve the existing full-gene-catalog offset / fractional-KO QL method.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(jsonlite)})
setDTthreads(2);set.seed(42);options(scipen=999)
OUT<-file.path(Sys.getenv("CHILE_PROJECT_ROOT", unset="."), "analysis/C03_C18_correction_2026-10-28/exclude_both")
IN<-file.path(OUT,'inputs');DER<-file.path(OUT,'derived');RES<-file.path(OUT,'results')
dir.create(RES,showWarnings=FALSE);stopifnot(file.exists(file.path(DER,'MATRICES_SUCCESS.json')))
write<-function(x,n)fwrite(x,file.path(RES,n),sep='\t',na='NA')
sets<-fromJSON(file.path(DER,'analysis_sets.json'));norms<-fread(file.path(DER,'TMM_offsets.tsv'))
meta<-as.data.frame(fread(file.path(IN,'metadata.tsv'),na.strings=c('','NA')))
kc<-as.data.table(readRDS(file.path(IN,'KO_fractional_raw_counts.rds')))
ann<-as.data.table(readRDS(file.path(IN,'KO_annotation.rds')))
audits<-list();tables<-list();summaries<-list();cpm_tables<-list();filtered<-list()
for(name in names(sets)){
 message('edgeR QL model: ',name)
 ss<-sets[[name]];m<-meta[match(ss,meta$sample_id),]
 is_surface<-grepl('surface',name)
 if(name=='B_HYPOTHETICAL_surface'){m$site[m$sample_id=='C18']<-'AZ';m$depth_label[m$sample_id=='C18']<-'5_10cm';m$site_pit[m$sample_id=='C18']<-'UNKNOWN_HYPOTHETICAL'}
 if(is_surface){
   m$group<-droplevels(factor(m$site,levels=c('AZ','SG','LC','NB')));formula<-~group
 }else{
   m$group<-factor(m$depth_label,levels=c('5_10cm','20_40cm','60_80cm'));m$pit<-droplevels(factor(m$site_pit));formula<-~pit+group
 }
 design<-model.matrix(formula,m);stopifnot(qr(design)$rank==ncol(design),nrow(design)-ncol(design)>=2)
 coef<-grep('^group',colnames(design));stopifnot(length(coef)==nlevels(m$group)-1)
 nt<-norms[analysis==name];nt<-nt[match(ss,nt$sample_id)];stopifnot(identical(nt$sample_id,ss))
 mat<-as.matrix(kc[,..ss]);rownames(mat)<-kc$KO
 y<-DGEList(counts=mat,lib.size=nt$full_catalog_library_size);y$samples$norm.factors<-nt$norm_factor
 keep<-filterByExpr(y,design=design,min.count=10)
 filtered[[name]]<-data.table(analysis=name,KO=rownames(mat),retained=keep)
 y<-y[keep,,keep.lib.sizes=TRUE]
 stopifnot(max(abs(as.numeric(getOffset(y))-log(nt$effective_library_size)))<1e-10)
 y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
 qlf<-glmQLFTest(fit,coef=coef)
 gettab<-function(q,test,contrast){
   t<-as.data.table(q$table,keep.rownames='KO');setnames(t,'logFC','log2FC',skip_absent=TRUE)
   if(!'log2FC'%in%names(t))t[,log2FC:=NA_real_]
   if(test=='omnibus')t[,log2FC:=NA_real_]
   t[,`:=`(analysis=name,test=test,contrast=contrast,FDR=p.adjust(PValue,'BH'))];t
 }
 t<-gettab(qlf,'omnibus',if(is_surface)'site' else 'depth_adjusted_for_pit')
 groupcpm<-cpm(y,normalized.lib.sizes=TRUE,log=FALSE)
 gt<-data.table(KO=rownames(groupcpm))
 for(gr in levels(m$group))gt[,(paste0('mean_CPM_',gr)):=rowMeans(groupcpm[,m$sample_id[m$group==gr],drop=FALSE])]
 gc<-grep('^mean_CPM_',names(gt),value=TRUE);gmat<-as.matrix(gt[,..gc])
 gt[,highest_group:=sub('mean_CPM_','',gc[max.col(gmat,ties.method='first')])]
 tests<-list(t)
 pairs<-combn(levels(m$group),2,simplify=FALSE)
 for(pair in pairs){
   vec<-rep(0,ncol(design));names(vec)<-colnames(design)
   if(pair[1]!=levels(m$group)[1])vec[paste0('group',pair[1])]<- -1
   if(pair[2]!=levels(m$group)[1])vec[paste0('group',pair[2])]<-1
   tests[[length(tests)+1]]<-gettab(glmQLFTest(fit,contrast=vec),'pairwise',paste(pair[2],'vs',pair[1],sep='_'))
 }
 tab<-rbindlist(tests,fill=TRUE);tab<-merge(tab,gt,by='KO',all.x=TRUE)
 tab<-merge(tab,ann,by='KO',all.x=TRUE)
 tables[[name]]<-tab[test=='omnibus']
 write(tab,paste0('KO_DA_',name,'_all_tests.tsv.gz'))
 write(data.table(sample_id=ss,design),paste0('KO_design_',name,'.tsv'))
 cpt<-as.data.table(groupcpm,keep.rownames='KO');write(cpt,paste0('KO_TMM_CPM_',name,'.tsv.gz'))
 summaries[[name]]<-tab[,.(tested=.N,FDR_lt_0_05=sum(FDR<.05),FDR_lt_0_10=sum(FDR<.10),min_FDR=min(FDR),raw_p_lt_0_05=sum(PValue<.05)),by=.(analysis,test,contrast)]
 audits[[name]]<-data.table(analysis=name,n=nrow(m),group_sizes=paste(names(table(m$group)),as.integer(table(m$group)),sep=':',collapse=';'),
  model=paste(deparse(formula),collapse=' '),rank=ncol(design),residual_df=nrow(design)-ncol(design),tested_KOs=nrow(y),common_dispersion=y$common.dispersion,
  min_norm_factor=min(nt$norm_factor),max_norm_factor=max(nt$norm_factor),normalization='subset-specific full-gene-catalog TMM offsets',
  status=if(name=='B_HYPOTHETICAL_surface')'HYPOTHETICAL_identity_and_independence_not_verified' else 'sensitivity')
}
all<-rbindlist(tables,fill=TRUE);write(rbindlist(audits),'KO_model_audit.tsv');write(rbindlist(summaries),'KO_test_summary.tsv')
write(rbindlist(filtered),'KO_filter_membership.tsv.gz')
write(all,'KO_omnibus_all_models.tsv.gz')

ref<-fread(file.path(IN,'reference/KO_DA_all_tests.tsv'))
valid<-list()
for(pair in list(c('reference_surface','surface_site'),c('reference_depth','NA_depth_pit_adjusted_sensitivity'))){
 x<-merge(all[analysis==pair[1]],ref[analysis==pair[2]&test_type=='omnibus'],by='KO',suffixes=c('_new','_old'))
 nn<-if(pair[1]=='reference_surface')7841 else 7801
 stopifnot(nrow(x)==nn,nrow(all[analysis==pair[1]])==nn)
 dif<-max(abs(x$FDR_new-x$FDR_old));df<-max(abs(x$F_new-x$F_old))
 stopifnot(dif<1e-7,df<1e-6)
 valid[[length(valid)+1]]<-data.table(analysis=pair[1],n=nn,max_abs_FDR_difference=dif,max_abs_F_difference=df,new_significant=sum(x$FDR_new<.05),reference_significant=sum(x$FDR_old<.05))
}
write(rbindlist(valid),'baseline_edgeR_validation.tsv')

changes<-list()
for(pair in list(c('A_surface','reference_surface'),c('AB_depth','reference_depth'),c('AB_complete_depth','reference_depth'))){
 old<-all[analysis==pair[2],.(KO,old_FDR=FDR,old_P=PValue,old_F=F,old_highest=highest_group)]
 new<-all[analysis==pair[1],.(KO,new_FDR=FDR,new_P=PValue,new_F=F,new_highest=highest_group)]
 z<-merge(old,new,by='KO',all=TRUE)
 z[,`:=`(analysis=pair[1],reference=pair[2],same_question=pair[1]!='A_surface')]
 z[,change:=fcase(is.na(old_FDR),'newly_tested',is.na(new_FDR),'filtered_after_subset',old_FDR<.05&new_FDR<.05,'retained_FDR_support',old_FDR<.05&new_FDR>=.05,'lost_FDR_support',old_FDR>=.05&new_FDR<.05,'gained_FDR_support',default='neither_FDR_supported')]
 changes[[length(changes)+1]]<-z
}
ch<-rbindlist(changes);write(ch,'KO_reference_vs_scenarios.tsv.gz');write(ch[,.N,by=.(analysis,same_question,change)],'KO_change_counts.tsv')
focal<-fread(file.path(IN,'reference/significant_curated_KOs.tsv'))$KO
ft<-all[KO%in%focal];ft[,paper_label:=fcase(KO=='K02040','pstS',KO=='K05922','hydB (original paper label)',KO=='K00436','hoxH',KO=='K17218','sqr',KO=='K06137','pqqC',default=KO)]
write(ft,'five_original_focal_markers.tsv')
writeLines(capture.output(sessionInfo()),file.path(RES,'edgeR_sessionInfo.txt'))
write_json(list(status='complete',models=names(sets),baseline_verified=TRUE,seed=42,completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(RES,'EDGER_SUCCESS.json'),pretty=TRUE,auto_unbox=TRUE)
message('All KO models and baseline verification COMPLETE')
