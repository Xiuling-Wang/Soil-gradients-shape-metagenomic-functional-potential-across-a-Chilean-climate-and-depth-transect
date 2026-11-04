#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(vegan);library(permute);library(jsonlite)})
setDTthreads(2);set.seed(42);options(scipen=999)
OUT<-file.path(Sys.getenv("CHILE_PROJECT_ROOT", unset="."), "analysis/C03_C18_correction_2026-10-28/exclude_both")
IN<-file.path(OUT,'inputs');DER<-file.path(OUT,'derived');RES<-file.path(OUT,'results')
dir.create(RES,showWarnings=FALSE)
stopifnot(file.exists(file.path(DER,'MATRICES_SUCCESS.json')))
write<-function(x,n)fwrite(x,file.path(RES,n),sep='\t',na='NA')
S<-c('C01',sprintf('C%02d',3:18));dl<-c('5_10cm','20_40cm','60_80cm')
meta<-as.data.frame(fread(file.path(IN,'metadata.tsv'),na.strings=c('','NA')))
meta<-meta[match(S,meta$sample_id),];stopifnot(identical(meta$sample_id,S))
meta$site<-factor(meta$site,levels=c('AZ','SG','LC','NB'));meta$depth_label<-factor(meta$depth_label,levels=dl)
meta$source_site<-as.character(meta$site);meta$source_pit<-meta$site_pit;meta$source_depth<-as.character(meta$depth_label)
sets<-fromJSON(file.path(DER,'analysis_sets.json'));dists<-readRDS(file.path(DER,'Hellinger_Bray_distances.rds'))
ko<-readRDS(file.path(DER,'KO_rpk.rds'));mag<-readRDS(file.path(DER,'MAG_rpk.rds'))
configs<-list(reference_surface=list(scenario='reference_original',scope='surface',ids=sets$reference_surface),
 reference_depth=list(scenario='reference_original',scope='depth',ids=sets$reference_depth),
 A_surface=list(scenario='A_exclude_C03_C18',scope='surface',ids=sets$A_surface),
 AB_depth=list(scenario='A_and_B_shared_depth',scope='depth',ids=sets$AB_depth),
 AB_complete_depth=list(scenario='A_and_B_complete_pit_sensitivity',scope='depth',ids=sets$AB_complete_depth))
getmeta<-function(name){
 spec<-configs[[name]];m<-meta[match(spec$ids,meta$sample_id),,drop=FALSE]
 if(name=='B_HYPOTHETICAL_surface'){
   i<-which(m$sample_id=='C18');m$site[i]<-'AZ';m$depth_label[i]<-'5_10cm';m$site_pit[i]<-'UNKNOWN_HYPOTHETICAL'
   chem<-c('pH','Conductivity','moisture','C','N','Feo','Alo','Mno','Sio','NH4','NO3','Pt')
   m[i,chem]<-NA_real_
 }
 m$site<-droplevels(m$site);m$depth_label<-droplevels(m$depth_label);m$site_pit<-factor(m$site_pit)
 m$group<-if(spec$scope=='surface')m$site else m$depth_label;m
}
permutations_of<-function(v){if(length(v)==1)return(matrix(v,1));do.call(rbind,lapply(seq_along(v),function(i)cbind(v[i],permutations_of(v[-i]))))}
blockperms<-function(blocks){
 positions<-split(seq_along(blocks),factor(blocks,levels=unique(blocks)))
 p<-lapply(positions,permutations_of);g<-expand.grid(lapply(p,function(z)seq_len(nrow(z))))
 ans<-matrix(NA_integer_,nrow(g),length(blocks))
 for(i in seq_len(nrow(g))) for(j in seq_along(p)) ans[i,positions[[j]]]<-p[[j]][g[i,j],]
 ans[!apply(ans,1,function(z)all(z==seq_along(blocks))),,drop=FALSE]
}
perm_for<-function(name){if(configs[[name]]$scope=='surface')return(999L);blockperms(as.character(getmeta(name)$site_pit))}
for(name in names(configs)){
 m<-getmeta(name)
 write(data.table(analysis=name,sample_id=m$sample_id,analysis_site=as.character(m$site),analysis_depth=as.character(m$depth_label),analysis_pit=as.character(m$site_pit),source_site=m$source_site,source_depth=m$source_depth,source_pit=m$source_pit),paste0('sample_overlay_',name,'.tsv'))
 if(configs[[name]]$scope=='depth')write(as.data.table(perm_for(name)),paste0('permutations_',name,'.tsv'))
}
stopifnot(nrow(perm_for('AB_depth'))==71L,nrow(perm_for('AB_complete_depth'))==35L)
subsetd<-function(d,ids)as.dist(as.matrix(d)[ids,ids,drop=FALSE])
diagnostics<-list()
one_test<-function(d,m,perms,module=FALSE){
 warnings_seen<-character()
 capture_warning<-function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')}
 set.seed(42);a<-adonis2(d~group,data=m,permutations=perms)
 b<-withCallingHandlers(if(module)betadisper(d,m$group,type='median',add='lingoes') else betadisper(d,m$group,type='median',bias.adjust=TRUE),warning=capture_warning)
 within_ss<-sum((b$distances-ave(b$distances,b$group,FUN=mean))^2)
 distance_ss<-sum(b$distances^2)
 degenerate<-all(table(b$group)<=2L)||within_ss<=max(distance_ss,1)*1e-20
 dispersion_status<-if(degenerate)'not_estimable_degenerate_within_group_distances' else 'estimable'
 df<-dp<-NA_real_
 if(!degenerate){
   set.seed(42);bt<-withCallingHandlers(permutest(b,permutations=perms),warning=capture_warning)
   df<-bt$tab[1,'F'];dp<-bt$tab[1,'Pr(>F)']
 }
 diagnostics[[length(diagnostics)+1L]]<<-data.table(samples=paste(m$sample_id,collapse=';'),module_method=module,
   dispersion_status=dispersion_status,within_group_distance_ss=within_ss,distance_ss=distance_ss,
   warnings=paste(unique(warnings_seen),collapse=' | '))
 data.table(n=nrow(m),groups=paste(names(table(m$group)),as.integer(table(m$group)),sep=':',collapse=';'),
 R2=a$R2[1],F=a$F[1],p=a$`Pr(>F)`[1],dispersion_F=df,dispersion_p=dp,dispersion_status=dispersion_status,
 permutation_rows=if(is.matrix(perms))nrow(perms)else perms)
}
classify<-function(p,q,dp,dq){
 ans<-ifelse(p<.05,'nominal-only','unresolved/precision-limited')
 ans[!is.na(q)&q<.05&!is.na(dp)&!is.na(dq)]<-'supported'
 ans[!is.na(q)&q<.05&((!is.na(dp)&dp<.05)|(!is.na(dq)&dq<.05))]<-'dispersion-entangled'
 ans[!is.na(q)&q<.05&(is.na(dp)|is.na(dq))]<-'unresolved/precision-limited'
 ans
}

message('Direct composition and dispersion: reference and scenario subsets')
rows<-list()
for(name in names(configs))for(pr in names(dists)){
 m<-getmeta(name);x<-one_test(subsetd(dists[[pr]],m$sample_id),m,perm_for(name))
 x[,`:=`(analysis=name,scenario=configs[[name]]$scenario,scope=configs[[name]]$scope,profile=pr,method='Hellinger_Bray')]
 rows[[length(rows)+1]]<-x
}
direct<-rbindlist(rows);direct[,`:=`(q=p.adjust(p,'BH'),dispersion_q=p.adjust(dispersion_p,'BH')),by=analysis]
direct[,evidence_status:=classify(p,q,dispersion_p,dispersion_q)]
write(direct,'composition_Hellinger_Bray.tsv')
# Four-group R2 for scenario A is descriptive only: AZ has one library.
am<-meta[meta$depth_label=='5_10cm'&meta$sample_id!='C03',];am$group<-am$site
adesc<-rbindlist(lapply(names(dists),function(pr){a<-adonis2(subsetd(dists[[pr]],am$sample_id)~group,data=am,permutations=0);data.table(profile=pr,n=nrow(am),R2=a$R2[1],F=a$F[1],p=NA_real_,status='DESCRIPTIVE_ONLY_AZ_singleton')}))
write(adesc,'A_four_group_descriptive_R2_NO_inference.tsv')

message('rCLR sensitivity: re-estimating zero completion per subset')
rows<-list();rclr<-list(KO=ko,MAG=mag)
for(name in names(configs))for(pr in names(rclr)){
 m<-getmeta(name);x<-t(rclr[[pr]][,m$sample_id,drop=FALSE]);x<-x[,colSums(x)>0,drop=FALSE]
 rank<-as.integer(min(3,nrow(x)-1,ncol(x)-1));set.seed(42)
 z<-decostand(x,method='rclr',impute=TRUE,ropt=rank);stopifnot(all(is.finite(z)))
 r<-one_test(dist(z),m,perm_for(name));r[,`:=`(analysis=name,scenario=configs[[name]]$scenario,scope=configs[[name]]$scope,profile=pr,rank=rank,zero_fraction=mean(x==0),method='rCLR_Aitchison')]
 rows[[length(rows)+1]]<-r
}
rc<-rbindlist(rows);rc[,`:=`(q=p.adjust(p,'BH'),dispersion_q=p.adjust(dispersion_p,'BH')),by=analysis]
rc[,evidence_status:=classify(p,q,dispersion_p,dispersion_q)];write(rc,'composition_rCLR_sensitivity.tsv')

message('Original resource module tests and exact depth permutations')
ft<-fread(file.path(IN,'curated_feature_sample_TPM.tsv'))
ft<-ft[tier=='primary'];mods<-unique(ft$module);stopifnot(length(mods)==3)
rows<-list()
for(variant in unique(ft$analysis_variant))for(mod in mods)for(name in names(configs)){
 m<-getmeta(name);t<-ft[analysis_variant==variant&module==mod&sample_id%in%m$sample_id,.(TPM=sum(TPM)),by=.(sample_id,marker)]
 w<-dcast(t,sample_id~marker,value.var='TPM',fill=0);w<-w[match(m$sample_id,w$sample_id)]
 x<-as.matrix(w[,-1]);x<-x[,colSums(x)>0,drop=FALSE];stopifnot(all(rowSums(x)>0))
 d<-vegdist(decostand(x,'hellinger'),'bray');r<-one_test(d,m,perm_for(name),module=TRUE)
 r[,`:=`(analysis=name,variant=variant,module=mod,features=ncol(x))];rows[[length(rows)+1]]<-r
}
mr<-rbindlist(rows);mf<-list()
for(sc in c('reference_original','A_exclude_C03_C18')){
 surfcase<-c(reference_original='reference_surface',A_exclude_C03_C18='A_surface',B_HYPOTHETICAL_C18_AZ_exclude_C03='B_HYPOTHETICAL_surface')[[sc]]
 depcase<-if(sc=='reference_original')'reference_depth' else 'AB_depth'
 x<-copy(mr[analysis%in%c(surfcase,depcase)]);x[,scenario:=sc]
 x[,`:=`(q_six_tests=p.adjust(p,'BH'),dispersion_q_six_tests=p.adjust(dispersion_p,'BH')),by=variant]
 x[,evidence_status:=classify(p,q_six_tests,dispersion_p,dispersion_q_six_tests)]
 mf[[length(mf)+1]]<-x
}
write(rbindlist(mf),'resource_modules_original_six_test_families.tsv')
x<-copy(mr[analysis=='AB_complete_depth']);x[,`:=`(q_three_modules=p.adjust(p,'BH'),dispersion_q_three_modules=p.adjust(dispersion_p,'BH')),by=variant]
write(x,'resource_modules_complete_pits_sensitivity.tsv')

message('Matched endpoint summaries and sample-level category rankings')
cats<-fread(file.path(IN,'curated_category_sample_TPM.tsv'))[tier=='primary',.(module,category,sample_id,value=category_TPM)]
cats[category=='DNRA',category:='Nitrite-to-ammonium reduction']
cats<-rbind(cats,data.table(module='nitrogen_transformations',category='N fixation',sample_id=S,value=0))
stopifnot(uniqueN(cats$category)==16)
mt<-fread(file.path(IN,'curated_module_total_sample_TPM.tsv'))[tier=='primary',.(module,category='MODULE_TOTAL',sample_id,value=module_total_TPM)]
cm<-rbind(cats,mt)
ep<-list()
for(pit in 1:3){ss<-paste0('C',9+pit);dd<-paste0('C',15+pit)
 x<-merge(cm[sample_id==ss,.(module,category,surface=value)],cm[sample_id==dd,.(module,category,deep=value)],by=c('module','category'))
 x[,`:=`(pit=paste0('NA',pit),surface_id=ss,deep_id=dd,percent_change=ifelse(surface>0,100*(deep-surface)/surface,NA_real_),direction=fifelse(deep<surface,'decrease',fifelse(deep>surface,'increase','unchanged')))]
 ep[[pit]]<-x
}
ep<-rbindlist(ep);write(ep,'category_module_endpoints_all_original_pits.tsv')
egs<-rbindlist(lapply(c('reference_three_pits','AB_two_complete_pits'),function(sc){
 x<-ep[if(sc=='reference_three_pits') rep(TRUE,.N) else pit!='NA3']
 x[,.(scenario=sc,n_pairs=.N,surface_mean=mean(surface),deep_mean=mean(deep),
 percent_change=ifelse(mean(surface)>0,100*(mean(deep)-mean(surface))/mean(surface),NA_real_),
 n_decrease=sum(deep<surface),n_increase=sum(deep>surface),n_equal=sum(deep==surface)),by=.(module,category)]
}))
write(egs,'category_module_matched_endpoint_summary.tsv')
refep<-fread(file.path(IN,'reference/Category_pit_endpoints.tsv'))
cc<-merge(ep[category!='MODULE_TOTAL'],refep,by=c('module','category','pit'))
stopifnot(nrow(cc)==48,max(abs(cc$surface-cc$surface_TPM))<1e-7,max(abs(cc$deep-cc$deep_TPM))<1e-7)
means<-list()
for(sc in c('reference_original','A_exclude_C03_C18')){
 m<-meta;if(sc!='reference_original')m<-m[m$sample_id!='C03',]
 if(sc=='A_exclude_C03_C18')m<-m[m$sample_id!='C18',]
 if(grepl('HYPOTHETICAL',sc)){m$site[m$sample_id=='C18']<-'AZ';m$depth_label[m$sample_id=='C18']<-'5_10cm'}
 mm<-data.table(sample_id=m$sample_id,group=ifelse(m$site=='NB',paste('NA',m$depth_label,sep='_'),as.character(m$site)))
 x<-merge(cm,mm,by='sample_id')
 means[[sc]]<-x[,.(scenario=sc,n=.N,mean_TPM=mean(value),sd_TPM=sd(value),samples=paste(sample_id,collapse=';')),by=.(module,category,group)]
}
write(rbindlist(means),'category_module_group_means.tsv')

message('Environmental associations: C18 has no soil chemistry in hypothetical AZ')
EV<-c('pH','Conductivity','moisture','C','N','Feo','Alo','Mno','Sio','NH4','NO3','Pt')
CV<-setdiff(EV,c('C','N','Feo','Alo','Sio','Pt'))
envrows<-list();pcalogs<-list();imputelogs<-list()
for(sc in c('reference_original','AB_verified_DNA_soil_pairs')) for(mode in c('original_depth_matched','complete_variables')){
 m<-meta;if(sc!='reference_original')m<-m[!m$sample_id%in%c('C03','C18'),]
 vars<-if(mode=='complete_variables')CV else EV
 if(mode=='original_depth_matched')for(v in EV)for(i in which(is.na(m[[v]]))){
   donor<-which(m$site=='NB' & m$depth_label==m$depth_label[i] & !is.na(m[[v]]))
   stopifnot(length(donor)>0,as.character(m$site[i])=='NB')
   m[[v]][i]<-mean(m[[v]][donor])
   imputelogs[[length(imputelogs)+1]]<-data.table(scenario=sc,sample_id=m$sample_id[i],variable=v,donors=paste(m$sample_id[donor],collapse=';'),n_donors=length(donor),value=m[[v]][i])
 }
 stopifnot(!anyNA(m[,vars]));pc<-prcomp(m[,vars],center=TRUE,scale.=TRUE)
 pcalogs[[length(pcalogs)+1]]<-data.table(scenario=sc,mode=mode,PC=seq_along(pc$sdev),fraction=pc$sdev^2/sum(pc$sdev^2))
 mm<-data.frame(sample_id=m$sample_id,SoilPC1=pc$x[,1],SoilPC2=pc$x[,2]);rownames(mm)<-m$sample_id
 for(scope in c('all','surface','depth'))for(pr in c('gene','KO')){
   sub<-if(scope=='all')m else if(scope=='surface')m[m$depth_label=='5_10cm',] else m[m$site=='NB',]
   d<-subsetd(dists[[pr]],sub$sample_id);pred<-mm[sub$sample_id,]
   f<-dbrda(d~SoilPC1+SoilPC2,data=pred);r<-RsquareAdj(f);p<-NA_real_;nperm<-0L
   if(scope!='all'){
     perms<-if(scope=='surface')999L else blockperms(as.character(sub$site_pit));set.seed(42)
     p<-anova(f,permutations=perms)$`Pr(>F)`[1];nperm<-if(is.matrix(perms))nrow(perms)else perms
   }
   envrows[[length(envrows)+1]]<-data.table(scenario=sc,mode=mode,scope=scope,profile=pr,n=nrow(sub),
     n_variables=length(vars),R2=r$r.squared,adjusted_R2=r$adj.r.squared,p=p,permutation_rows=nperm)
 }
}
er<-rbindlist(envrows);er[,q:=p.adjust(p,'BH'),by=.(scenario,mode,scope)]
write(er,'environment_dbRDA.tsv');write(rbindlist(imputelogs),'environment_imputation_donors.tsv');write(rbindlist(pcalogs),'environment_PCA_variance.tsv')

message('Cross-resolution geometric agreement and MAG distribution summaries')
cr<-list()
for(name in names(configs)) for(pair in list(c('gene','KO'),c('gene','MAG'),c('KO','MAG'))){
 m<-getmeta(name);xx<-wcmdscale(subsetd(dists[[pair[1]]],m$sample_id),k=2,eig=TRUE)$points
 yy<-wcmdscale(subsetd(dists[[pair[2]]],m$sample_id),k=2,eig=TRUE)$points
 set.seed(42);t<-protest(xx,yy,permutations=perm_for(name),symmetric=TRUE)
 cr[[length(cr)+1]]<-data.table(analysis=name,pair=paste(pair,collapse='_'),n=nrow(m),r=unname(t$t0),p=t$signif)
}
cr<-rbindlist(cr);cr[,q:=p.adjust(p,'BH'),by=analysis];write(cr,'cross_resolution_procrustes.tsv')
tax<-fread(file.path(IN,'mag_taxonomy.tsv'));tax<-tax[match(rownames(mag),tax[[1]])]
rel<-sweep(mag,2,colSums(mag),'/');phy<-rowsum(rel,tax$phylum,reorder=FALSE)*100
write(data.table(phylum=rownames(phy),phy),'MAG_phylum_percent.tsv')
arch<-data.table(sample_id=S,MAG005=rel['MAG005',]*100,MAG070=rel['MAG070',]*100,
 archaeal_total=colSums(rel[tax$Domain=='Archaea',,drop=FALSE])*100)
write(arch,'archaeal_MAG_sample_percent.tsv')
write(rbindlist(diagnostics),'dispersion_diagnostics.tsv')
writeLines(capture.output(sessionInfo()),file.path(RES,'statistics_sessionInfo.txt'))
write_json(list(status='complete',seed=42,scenarios=names(configs),baseline_endpoint_max_abs_difference=max(abs(cc$deep-cc$deep_TPM)),completed=as.character(Sys.time())),file.path(RES,'STATISTICS_SUCCESS.json'),pretty=TRUE,auto_unbox=TRUE)
message('Composition, module, environmental and descriptive scenario tests COMPLETE')
