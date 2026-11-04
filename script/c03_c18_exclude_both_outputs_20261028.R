#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(vegan);library(edgeR);library(jsonlite)})
setDTthreads(2);set.seed(42);options(scipen=999)
ROOT<-normalizePath('.')
OLD<-file.path(ROOT,'analysis/C03_C18_correction_2026-10-28')
OUT<-file.path(OLD,'exclude_both');IN<-file.path(OUT,'inputs');DER<-file.path(OUT,'derived')
RES<-file.path(OUT,'results');DEST<-file.path(OUT,'adopted');QA<-file.path(OUT,'qa')
write<-function(x,n)fwrite(x,file.path(DEST,n),sep='\t',na='NA')
ids<-c('C01',sprintf('C%02d',4:17));excluded<-c('C03','C18')
analyses<-c('A_surface','AB_depth','AB_complete_depth')
for(f in c('STATISTICS_SUCCESS.json','EDGER_SUCCESS.json'))stopifnot(fromJSON(file.path(RES,f))$status=='complete')
stopifnot(fromJSON(file.path(DER,'MATRICES_SUCCESS.json'))$status=='complete')
stopifnot(fromJSON(file.path(QA,'curated_recalculation/SUCCESS.json'))$status=='complete')

message('Comparing freshly computed tables with the prior A run')
regress<-list()
files<-c('composition_Hellinger_Bray.tsv','composition_rCLR_sensitivity.tsv','KO_test_summary.tsv',
 'resource_modules_original_six_test_families.tsv','environment_dbRDA.tsv',
 'environment_imputation_donors.tsv','category_module_matched_endpoint_summary.tsv',
 'five_original_focal_markers.tsv','cross_resolution_procrustes.tsv',
 paste0('KO_DA_',analyses,'_all_tests.tsv.gz'))
for(f in files){
 a<-fread(file.path(RES,f));b<-fread(file.path(OLD,'results',f))
 if('analysis'%in%names(a))b<-b[analysis%in%unique(a$analysis)]
 if('scenario'%in%names(a))b<-b[scenario%in%unique(a$scenario)]
 stopifnot(identical(names(a),names(b)),nrow(a)==nrow(b))
 nums<-names(a)[vapply(a,is.numeric,logical(1))];other<-setdiff(names(a),nums)
 stopifnot(identical(as.data.frame(a[,..other]),as.data.frame(b[,..other])))
 stopifnot(identical(is.na(as.matrix(a[,..nums])),is.na(as.matrix(b[,..nums]))))
 delta<-max(c(0,unlist(lapply(nums,function(n)abs(a[[n]]-b[[n]])))),na.rm=TRUE)
 regress[[f]]<-data.table(file=f,rows=nrow(a),max_numeric_difference=delta,pass=delta<1e-8)
}
rv<-rbindlist(regress);fwrite(rv,file.path(QA,'historical_A_regression_validation.tsv'),sep='\t')
stopifnot(all(rv$pass))

message('Exporting the exact 15-column matrices and verifying each retained value')
meta<-fread(file.path(IN,'metadata.tsv'),na.strings=c('','NA'))[match(ids,sample_id)]
stopifnot(identical(meta$sample_id,ids))
meta[,`:=`(source_site=site,site=fifelse(site=='NB','NA',site),pit_display=sub('^NB','NA',site_pit))]
write(meta,'metadata_15_libraries.tsv')
matrix_checks<-list()
for(layer in c('gene','KO','MAG')){
 x<-switch(layer,gene=fread(file.path(IN,'gene_counts.tsv')),
  KO=as.data.table(readRDS(file.path(IN,'KO_fractional_raw_counts.rds'))),MAG=fread(file.path(IN,'mag_counts.tsv')))
 cols<-c(names(x)[1],ids);y<-x[,..cols]
 stopifnot(identical(names(y)[-1],ids),all(vapply(ids,function(s)identical(y[[s]],x[[s]]),logical(1))))
 write(y,paste0(layer,'counts_15_libraries.tsv.gz'))
 matrix_checks[[layer]]<-data.table(layer=layer,rows=nrow(y),libraries=15,max_value_change=0,all_ids_unchanged=TRUE)
 if(layer=='gene'){
   lens<-fread(file.path(IN,'gene_lengths.tsv'),header=FALSE,col.names=c('GeneName','Length'))
   ll<-lens$Length[match(x$GeneName,lens$GeneName)];stopifnot(!anyNA(ll),all(ll>0))
   den<-data.table(sample_id=ids,gene_RPK_total=vapply(ids,function(s)sum(x[[s]]/(ll/1000)),numeric(1)))
   write(den,'full_catalog_TPM_denominators.tsv')
 }
 rm(x,y);invisible(gc())
}
write(rbindlist(matrix_checks),'matrix_value_validation.tsv')
for(f in c('curated_feature_sample_TPM.tsv','curated_category_sample_TPM.tsv','curated_module_total_sample_TPM.tsv')){
 x<-fread(file.path(IN,f))[sample_id%in%ids];x[,site:=fifelse(site=='NB','NA',site)];write(x,f)
}
ft<-fread(file.path(DEST,'curated_feature_sample_TPM.tsv'))[tier=='primary'&analysis_variant=='primary_marker_fractional']
ko_rpk<-readRDS(file.path(DER,'KO_rpk.rds'))[,ids,drop=FALSE]
ko_tpm<-sweep(ko_rpk,2,den$gene_RPK_total,'/')*1e6
write(data.table(KO=rownames(ko_tpm),ko_tpm),'KO_TPM_15_libraries.tsv.gz')
kt<-ft[marker%in%rownames(ko_tpm)]
kt[,recomputed_TPM:=ko_tpm[cbind(match(marker,rownames(ko_tpm)),match(sample_id,colnames(ko_tpm)))]]
kt[,absolute_difference:=abs(TPM-recomputed_TPM)]
kt[,note:='Different estimands: curated name-filtered within-module allocation versus whole-KO catalog allocation']
fwrite(kt,file.path(QA,'whole_KO_vs_curated_marker_TPM_NOT_EQUIVALENT.tsv'),sep='\t')
if(file.exists(file.path(DEST,'curated_KO_TPM_recalculation_validation.tsv')))
 file.rename(file.path(DEST,'curated_KO_TPM_recalculation_validation.tsv'),file.path(QA,'initial_non_equivalent_TPM_check.tsv'))

message('Publishing A-only statistics; keeping original baselines in QA results')
for(f in c('composition_Hellinger_Bray.tsv','composition_rCLR_sensitivity.tsv','KO_test_summary.tsv','KO_model_audit.tsv','five_original_focal_markers.tsv','cross_resolution_procrustes.tsv','KO_omnibus_all_models.tsv.gz','KO_filter_membership.tsv.gz')){
 x<-fread(file.path(RES,f))[analysis%in%analyses]
 for(col in intersect(c('contrast','groups','group_sizes','highest_group'),names(x)))x[,(col):=gsub('NB','NA',get(col),fixed=TRUE)]
 if('mean_CPM_NB'%in%names(x))setnames(x,'mean_CPM_NB','mean_CPM_NA')
 write(x,f)
}
for(a in analyses){
 for(prefix in c('KO_DA_','KO_TMM_CPM_','KO_design_')){
   f<-paste0(prefix,a,if(prefix=='KO_DA_')'_all_tests.tsv.gz' else if(prefix=='KO_TMM_CPM_')'.tsv.gz' else '.tsv')
   x<-fread(file.path(RES,f));if('mean_CPM_NB'%in%names(x))setnames(x,'mean_CPM_NB','mean_CPM_NA')
   if('contrast'%in%names(x))x[,contrast:=gsub('NB','NA',contrast,fixed=TRUE)]
   if('highest_group'%in%names(x))x[highest_group=='NB',highest_group:='NA']
   write(x,f)
 }
}
for(f in c('environment_dbRDA.tsv','environment_imputation_donors.tsv','environment_PCA_variance.tsv'))write(fread(file.path(RES,f))[scenario=='AB_verified_DNA_soil_pairs'],f)
modules<-fread(file.path(RES,'resource_modules_original_six_test_families.tsv'))[scenario=='A_exclude_C03_C18']
modules[,groups:=gsub('NB:','NA:',groups,fixed=TRUE)]
write(modules,'resource_modules_six_test_families.tsv')
write(fread(file.path(RES,'resource_modules_complete_pits_sensitivity.tsv')),'resource_modules_complete_pits_sensitivity.tsv')
write(fread(file.path(RES,'category_module_matched_endpoint_summary.tsv'))[scenario=='AB_two_complete_pits'],'category_module_matched_endpoint_summary.tsv')
write(fread(file.path(RES,'category_module_endpoints_all_original_pits.tsv'))[pit!='NA3'],'category_module_complete_pit_endpoints.tsv')
means<-fread(file.path(RES,'category_module_group_means.tsv'))[scenario=='A_exclude_C03_C18']
means[,summary_type:=fifelse(n==1,'single_observation','group_mean')]
means[n==1,`:=`(single_value_TPM=mean_TPM,mean_TPM=NA_real_,sd_TPM=NA_real_)]
write(means,'category_module_group_summaries.tsv')
ds<-readRDS(file.path(DER,'Hellinger_Bray_distances.rds'))
for(n in names(ds))write(data.table(sample_id=ids,as.matrix(ds[[n]])[ids,ids]),paste0(n,'_Hellinger_Bray_15.tsv'))

message('Recomputing AZ descriptions, without AZ hypothesis tests')
surface<-meta[depth_label=='5_10cm']$sample_id
norms<-fread(file.path(DER,'TMM_offsets.tsv'));nt<-norms[analysis=='AZ_descriptive_surface'][match(surface,sample_id)]
kc<-as.data.table(readRDS(file.path(IN,'KO_fractional_raw_counts.rds')))
cpm<-sweep(as.matrix(kc[,..surface]),2,nt$effective_library_size,'/')*1e6;rownames(cpm)<-kc$KO
write(data.table(KO=rownames(cpm),cpm),'KO_TMM_CPM_surface10_DESCRIPTIVE.tsv.gz')
write(nt,'AZ_descriptive_full_catalog_TMM_offsets.tsv')
ranker<-function(x,keys){
 r<-x[,{
  other<-value[sample_id!='C01'];val<-value[sample_id=='C01'];mx<-max(other)
  list(C01_value=val,C01_rank_descending=min(rank(-value,ties.method='min')[sample_id=='C01']),
   other_max=mx,other_max_samples=paste(sample_id[sample_id!='C01'&value==mx],collapse=';'),
   ratio_C01_to_other_max=if(mx>0)val/mx else NA_real_,interpretation='DESCRIPTIVE_ONLY_AZ_n1')
 },by=keys];r
}
tm<-ft[sample_id%in%surface,.(module,category,marker,marker_label,sample_id,value=TPM)]
write(tm,'AZ_surface_marker_TPM_all_samples.tsv');write(ranker(tm,c('module','category','marker','marker_label')),'AZ_surface_marker_TPM_ranks.tsv')
allcpm<-melt(data.table(KO=rownames(cpm),cpm),id.vars='KO',variable.name='sample_id',value.name='value')
reg<-unique(ft[,.(module,category,marker,marker_label)])
cm<-merge(reg,allcpm,by.x='marker',by.y='KO');write(cm,'AZ_surface_marker_TMM_CPM_all_samples.tsv')
write(ranker(cm,c('module','category','marker','marker_label')),'AZ_surface_marker_TMM_CPM_ranks.tsv')
all_tests<-fread(file.path(DEST,'KO_omnibus_all_models.tsv.gz'),na.strings='')
curated_tests<-merge(unique(reg[,.(marker,marker_label,module,category)]),all_tests,by.x='marker',by.y='KO')
write(curated_tests,'all_curated_marker_KO_tests.tsv')
write(curated_tests[analysis=='A_surface'&FDR<.05],'surface_curated_markers_FDR_supported.tsv')
write(curated_tests[,.(registered_markers_tested=uniqueN(marker),FDR_supported=sum(FDR<.05)),by=analysis],
 'curated_marker_test_counts.tsv')
for(level in c('category','module')){
 f<-if(level=='category')'curated_category_sample_TPM.tsv' else 'curated_module_total_sample_TPM.tsv'
 x<-fread(file.path(DEST,f))[tier=='primary'&sample_id%in%surface]
 valcol<-if(level=='category')'category_TPM' else 'module_total_TPM';setnames(x,valcol,'value')
 write(ranker(x,if(level=='category')c('module','category') else 'module'),paste0('AZ_surface_',level,'_TPM_ranks.tsv'))
}
distrows<-list();distcmp<-list()
for(pr in names(ds)){
 d<-as.matrix(ds[[pr]])[surface,surface];pairs<-as.data.table(as.table(d));setnames(pairs,c('sample1','sample2','distance'))
 pairs<-pairs[match(sample1,surface)<match(sample2,surface)]
 pairs[,`:=`(site1=meta$site[match(sample1,meta$sample_id)],site2=meta$site[match(sample2,meta$sample_id)],profile=pr)]
 distrows[[pr]]<-pairs
 distcmp[[pr]]<-data.table(profile=pr,C01_min_distance=min(d['C01',setdiff(surface,'C01')]),
  C01_nearest_sample=setdiff(surface,'C01')[which.min(d['C01',setdiff(surface,'C01')])],
  replicated_within_site_min=min(pairs[site1==site2,distance]),replicated_within_site_max=max(pairs[site1==site2,distance]))
}
write(rbindlist(distrows),'surface_pairwise_distances_DESCRIPTIVE.tsv');write(rbindlist(distcmp),'C01_distance_context_DESCRIPTIVE.tsv')

message('Registered depth normalization diagnostic, same eight libraries and KO set')
ss<-meta[site=='NA']$sample_id;md<-as.data.frame(meta[match(ss,sample_id)])
md$pit<-factor(md$site_pit);md$group<-factor(md$depth_label,levels=c('5_10cm','20_40cm','60_80cm'))
design<-model.matrix(~pit+group,md)
keepids<-fread(file.path(DEST,'KO_filter_membership.tsv.gz'))[analysis=='AB_depth'&retained==TRUE,KO]
mat<-as.matrix(kc[match(keepids,KO),..ss]);rownames(mat)<-keepids
nt<-norms[analysis=='reference_depth'][match(ss,sample_id)]
y<-DGEList(mat,lib.size=nt$full_catalog_library_size);y$samples$norm.factors<-nt$norm_factor
stopifnot(max(abs(getOffset(y)-log(nt$effective_library_size)))<1e-10)
y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
z<-as.data.table(glmQLFTest(fit,coef=grep('^group',colnames(design)))$table,keep.rownames='KO')
z[,`:=`(FDR=p.adjust(PValue,'BH'),diagnostic='same_8_libraries_same_KOs_original_9_library_offsets')]
write(z,'depth_frozen_original_offsets_DIAGNOSTIC.tsv.gz')
primary<-fread(file.path(DEST,'KO_omnibus_all_models.tsv.gz'))[analysis=='AB_depth']
joint<-merge(primary[,.(KO,primary_FDR=FDR)],z[,.(KO,frozen_offset_FDR=FDR)],by='KO')
write(joint,'depth_normalization_DIAGNOSTIC_comparison.tsv.gz')
write(data.table(model=c('primary_subset_TMM','diagnostic_original_offsets'),tested=c(nrow(primary),nrow(z)),
 FDR_lt_0_05=c(sum(primary$FDR<.05),sum(z$FDR<.05)),both_supported=sum(joint$primary_FDR<.05&joint$frozen_offset_FDR<.05)),
 'depth_normalization_DIAGNOSTIC_summary.tsv')

message('Recomputing MAG summaries with fixed references and explicit provenance')
mag<-readRDS(file.path(DER,'MAG_rpk.rds'))[,ids,drop=FALSE];rel<-sweep(mag,2,colSums(mag),'/')
tax<-fread(file.path(IN,'mag_taxonomy.tsv'));tax<-tax[match(rownames(mag),tax[[1]])]
quality<-fread(file.path(IN,'mag_quality.tsv'));quality<-quality[match(rownames(mag),quality[['Bin Id']])]
prov<-data.table(MAG=rownames(mag),original_assembly_library=quality$Sample,
 source_status=ifelse(quality$Sample%in%excluded,'excluded_library_identity_unresolved','retained_library'),
 ecological_source=ifelse(quality$Sample%in%excluded,NA_character_,quality$Sample),reference_action='retain_fixed_reference',
 phylum=tax$phylum,completeness=quality$Completeness,contamination=quality$Contamination)
stopifnot(sum(prov$source_status=='excluded_library_identity_unresolved')==46L)
write(prov,'MAG_reference_provenance.tsv')
write(data.table(MAG=rownames(rel),rel*100),'MAG_relative_percent_15_libraries.tsv')
ph<-rowsum(rel,tax$phylum,reorder=FALSE)*100;write(data.table(phylum=rownames(ph),ph),'MAG_phylum_percent.tsv')
write(data.table(sample_id=ids,dominant_phylum=rownames(ph)[max.col(t(ph),ties.method='first')],dominant_percent=apply(ph,2,max)), 'MAG_dominant_phylum_by_sample.tsv')
arch<-which(tax$Domain=='Archaea')
write(data.table(sample_id=ids,MAG005=rel['MAG005',]*100,MAG070=rel['MAG070',]*100,archaeal_total=colSums(rel[arch,,drop=FALSE])*100),'archaeal_MAG_sample_percent.tsv')
groups<-ifelse(meta$site=='NA',paste0('NA_',meta$depth_label),meta$site)
gmeans<-sapply(unique(groups),function(g)rowMeans(rel[,ids[groups==g],drop=FALSE]))
write(data.table(MAG=rownames(rel),gmeans),'MAG_group_abundance_fractions.tsv')
write(data.table(MAG=rownames(rel),mean_fraction_15_libraries=rowMeans(rel)),'MAG_overall_abundance_fractions.tsv')
qc<-fread(file.path(ROOT,'data/processed/metagenome_output/stats_for_Xiuling/read_counts.tsv'))[Sample%in%ids]
write(qc,'read_QC_retained_libraries.tsv')
recruit<-fread(file.path(ROOT,'data/processed/metagenome_output/reports_Seimibin_for_Xiuling/genome_mapping/results_data/multiqc_samtools_stats.txt'))[match(ids,Sample)]
stopifnot(identical(recruit$Sample,ids))
rec<-recruit[,.(sample_id=Sample,QC_read_sequences=raw_total_sequences,MAG_mapped_read_sequences=reads_mapped,MAG_mapped_percent=100*reads_mapped/raw_total_sequences)]
stopifnot(max(abs(rec$MAG_mapped_percent-recruit$reads_mapped_percent))<1e-10)
write(rec,'MAG_read_recruitment_15_libraries.tsv')

message('Exporting exact soil PCA and dbRDA plotting inputs')
EV<-c('pH','Conductivity','moisture','C','N','Feo','Alo','Mno','Sio','NH4','NO3','Pt')
env<-copy(meta);donors<-fread(file.path(DEST,'environment_imputation_donors.tsv'))
for(i in seq_len(nrow(donors)))env[match(donors$sample_id[i],sample_id),(donors$variable[i]):=donors$value[i]]
stopifnot(!anyNA(env[,..EV]),!any(donors$sample_id%in%excluded))
pc<-prcomp(env[,..EV],center=TRUE,scale.=TRUE)
rownames(pc$x)<-ids
write(data.table(sample_id=ids,pc$x),'environment_PCA_scores.tsv')
write(data.table(variable=rownames(pc$rotation),pc$rotation),'environment_PCA_loadings.tsv')
write(env[,c('sample_id',EV),with=FALSE],'environment_imputed_values.tsv')
envlong<-melt(meta[,c('sample_id','site','depth_label',EV),with=FALSE],id.vars=c('sample_id','site','depth_label'),variable.name='variable')
envlong[,group:=ifelse(site=='NA',paste0('NA_',depth_label),site)]
em<-envlong[,.(observed_n=sum(!is.na(value)),value=if(all(is.na(value)))NA_real_ else mean(value,na.rm=TRUE)),by=.(variable,group)]
em[,row_z:=as.numeric(scale(value)),by=variable];write(em,'soil_observed_group_means.tsv')
for(scope in c('surface','depth')){
 ss<-if(scope=='surface')surface else meta[site=='NA']$sample_id
 d<-as.dist(as.matrix(ds$KO)[ss,ss]);pred<-data.frame(SoilPC1=pc$x[ss,1],SoilPC2=pc$x[ss,2],row.names=ss)
 mod<-dbrda(d~SoilPC1+SoilPC2,data=pred)
 coords<-scores(mod,display='sites',choices=1:2)
 write(data.table(sample_id=rownames(coords),coords),paste0('KO_',scope,'_dbRDA_coordinates.tsv'))
}
write_json(list(status='complete',libraries=ids,excluded_libraries=excluded,exclusion_is_not_identity_proof=TRUE,
 fixed_MAG_references=228,excluded_library_MAG_references=46,baseline_regression_max_difference=max(rv$max_numeric_difference),
 matrix_values_unchanged=TRUE,no_AZ_inference=TRUE,completed=as.character(Sys.time())),file.path(DEST,'ADOPTED_OUTPUTS_SUCCESS.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(QA,'outputs_sessionInfo.txt'))
message('15-library adopted outputs COMPLETE')
