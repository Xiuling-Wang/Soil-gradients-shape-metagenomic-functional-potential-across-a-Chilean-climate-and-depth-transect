#!/usr/bin/env Rscript
# Derived matrices and TMM only, from immutable copied inputs. No source writes.
suppressPackageStartupMessages({library(data.table);library(vegan);library(edgeR);library(jsonlite)})
Sys.setenv(OMP_NUM_THREADS=2,OPENBLAS_NUM_THREADS=2)
setDTthreads(2);set.seed(42);options(scipen=999)
ROOT<-normalizePath(Sys.getenv("CHILE_PROJECT_ROOT", unset="."))
OUT<-file.path(ROOT,'analysis/C03_C18_correction_2026-10-28/exclude_both')
IN<-file.path(OUT,'inputs');RES<-file.path(OUT,'derived');dir.create(RES,showWarnings=FALSE)
S<-c('C01',sprintf('C%02d',3:18))
meta<-as.data.frame(fread(file.path(IN,'metadata.tsv'),na.strings=c('','NA')))
meta<-meta[match(S,meta$sample_id),];stopifnot(identical(meta$sample_id,S),all(meta$dna_type=='iDNA'))
depth<-meta[meta$site=='NB',];surf<-meta[meta$depth_label=='5_10cm',]
sets<-list(reference_surface=surf$sample_id,reference_depth=depth$sample_id,
 A_surface=surf$sample_id[!surf$sample_id%in%c('C01','C03')],
 AB_depth=depth$sample_id[depth$sample_id!='C18'],
 AB_complete_depth=depth$sample_id[depth$site_pit%in%c('NB1','NB2')])
write_json(sets,file.path(RES,'analysis_sets.json'),pretty=TRUE,auto_unbox=FALSE)
message('Loading immutable gene counts and lengths')
counts<-fread(file.path(IN,'gene_counts.tsv'))
stopifnot(nrow(counts)==2422720L,!anyDuplicated(counts$GeneName),identical(names(counts)[-1],S))
lengths<-fread(file.path(IN,'gene_lengths.tsv'),header=FALSE,col.names=c('GeneName','Length'))
stopifnot(!anyDuplicated(lengths$GeneName))
lens<-lengths$Length[match(counts$GeneName,lengths$GeneName)]
stopifnot(!anyNA(lens),all(lens>0));rm(lengths);invisible(gc())
gene<-as.matrix(counts[,..S]);storage.mode(gene)<-'double'
stopifnot(all(is.finite(gene)),all(gene>=0));ids<-counts$GeneName;rm(counts);invisible(gc())

message('Checking frozen KO allocation and reconstructing length-normalized KO profile')
map<-as.data.table(readRDS(file.path(IN,'KO_fractional_map.rds')))
kc<-as.data.table(readRDS(file.path(IN,'KO_fractional_raw_counts.rds')))
stopifnot(!anyDuplicated(map[,.(GeneName,KO)]),all(map$n_valid_KO>=1),!anyDuplicated(kc$KO))
mi<-match(map$GeneName,ids);stopifnot(!anyNA(mi))
ko_rpk<-matrix(0,nrow(kc),length(S),dimnames=list(kc$KO,S))
ko_count_error<-numeric(length(S))
for(j in seq_along(S)){
 mapped<-gene[mi,j]/map$n_valid_KO
 cr<-rowsum(mapped,map$KO,reorder=FALSE)
 kr<-rowsum(mapped/(lens[mi]/1000),map$KO,reorder=FALSE)
 stopifnot(setequal(rownames(cr),kc$KO))
 ko_count_error[j]<-max(abs(cr[kc$KO,1]-kc[[S[j]]]))
 ko_rpk[,j]<-kr[kc$KO,1]
}
stopifnot(max(ko_count_error)<1e-6)
fwrite(data.table(sample_id=S,max_abs_KO_cached_count_error=ko_count_error),file.path(RES,'KO_cache_validation.tsv'),sep='\t')
saveRDS(ko_rpk,file.path(RES,'KO_rpk.rds'))
rm(map,mi,kc,mapped,cr,kr);invisible(gc())

message('Re-estimating full-catalog TMM for each fitted sample subset')
norms<-list()
normsets<-c(sets,list(AZ_descriptive_surface=surf$sample_id[surf$sample_id!='C03']))
for(n in names(normsets)){
 ss<-normsets[[n]];message('TMM: ',n,' (',length(ss),')')
 g<-gene[,ss,drop=FALSE]
 y<-DGEList(counts=g);y<-calcNormFactors(y,method='TMM')
 norms[[n]]<-data.table(analysis=n,sample_id=ss,full_catalog_library_size=y$samples$lib.size,
 norm_factor=y$samples$norm.factors,effective_library_size=y$samples$lib.size*y$samples$norm.factors)
 rm(y,g);invisible(gc())
}
fwrite(rbindlist(norms),file.path(RES,'TMM_offsets.tsv'),sep='\t')

message('Building exact original Hellinger-Bray distances (all features, no reassembly)')
gene_rpk<-gene/(lens/1000)
allzero<-rowSums(gene)==0
active_counts<-rbindlist(lapply(names(sets),function(n) data.table(analysis=n,
 positive_genes=sum(rowSums(gene[,sets[[n]],drop=FALSE])>0))))
fwrite(active_counts,file.path(RES,'positive_feature_counts.tsv'),sep='\t')
stopifnot(sum(allzero)==170L)
gene_rpk<-gene_rpk[!allzero,,drop=FALSE];rm(gene,lens,ids);invisible(gc())
make_dist<-function(m) vegdist(decostand(t(m),method='hellinger'),method='bray')
dg<-make_dist(gene_rpk);rm(gene_rpk);invisible(gc())
dk<-make_dist(ko_rpk)
mc<-fread(file.path(IN,'mag_counts.tsv'))
q<-fread(file.path(IN,'mag_quality.tsv'))
stopifnot(setequal(mc$Sample,q[['Bin Id']]))
sizes<-q$Genome_Size[match(mc$Sample,q[['Bin Id']])]
mag_rpk<-as.matrix(mc[,..S])/(sizes/1000);rownames(mag_rpk)<-mc$Sample
dm<-make_dist(mag_rpk)
saveRDS(mag_rpk,file.path(RES,'MAG_rpk.rds'))
distances<-list(gene=dg,KO=dk,MAG=dm)
saveRDS(distances,file.path(RES,'Hellinger_Bray_distances.rds'))
for(n in names(distances)) fwrite(data.table(sample_id=S,as.matrix(distances[[n]])),file.path(RES,paste0(n,'_Hellinger_Bray.tsv')),sep='\t')

# Validate baseline effect sizes before allowing scenario interpretation.
ref<-fread(file.path(IN,'reference/exact_direct_location_permanova.tsv'))
labels<-c(gene='Whole gene catalog',KO='Fractional KO catalog',MAG='Genome-size-normalised MAG assemblage')
validation<-list()
for(n in names(distances)) for(z in c('surface_site','nahuelbuta_depth')){
 mm<-if(z=='surface_site')surf else depth;d<-as.dist(as.matrix(distances[[n]])[mm$sample_id,mm$sample_id])
 group<-factor(if(z=='surface_site') mm$site else mm$depth_label)
 ad<-adonis2(d~group,permutations=0)
 old<-ref[profile==labels[[n]] & subset==z]
 stopifnot(nrow(old)==1)
 validation[[length(validation)+1]]<-data.table(profile=n,subset=z,current_R2=ad$R2[1],reference_R2=old$R2,
 abs_R2_difference=abs(ad$R2[1]-old$R2),current_F=ad$F[1],reference_F=old$F)
}
v<-rbindlist(validation);fwrite(v,file.path(RES,'baseline_profile_validation.tsv'),sep='\t')
stopifnot(max(v$abs_R2_difference)<1e-10)
writeLines(capture.output(sessionInfo()),file.path(RES,'R_sessionInfo.txt'))
write_json(list(status='complete',genes=2422720,positive_original_genes=2422550,KOs=nrow(ko_rpk),MAGs=nrow(mag_rpk),
 baseline_max_abs_R2_difference=max(v$abs_R2_difference),seed=42,completed=as.character(Sys.time())),file.path(RES,'MATRICES_SUCCESS.json'),pretty=TRUE,auto_unbox=TRUE)
message('Matrices, full-catalog TMM and baseline validation COMPLETE')
