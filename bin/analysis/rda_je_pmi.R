#####
#Analysis of relatedness, AF and Fst data from PMI per-samples per-tps
#####

#load our shit
#run this in the same directory 
pkg = c('tidyverse', 'geosphere', 'data.table', 'vegan' ,'cowplot', 'parallel', 'Rcpp', 'RColorBrewer', 'UpSetR', 'sparseAHC', 'moments', 'parallel')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#remotes::install_github("khabbazian/sparseAHC")
#set num cores

#####
#GLM to see impact of 
# - site, bioclim, water cover, land cover type, land usage, habitat, treatment, on pop structure (measured with relatedness, Fis, ROH and theta and sitewise fst)
# for PMI and JE data
#####



#load bioclim data
#get bioclimatic vars from worldclim
bioclimdata <- raster::getData("worldclim",var="bio",res=10)
bioclimdata <- bioclimdata[[1:19]]
#get landcover data
landcov = raster::raster('~/Downloads/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif') #load raster file
#get water data


#for pmi and je metadata, use raster library to extract info from above data for each site

#for je
je_metadata=read.csv('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/metadata/metadata_full_sequencedonly_062022.csv')
#get coords from metadata and make into raster spatialpoints obj
je_coords <- data.frame(x=je_metadata$Long, y=je_metadata$Lat)
je_points <- sp::SpatialPoints(je_coords, proj4string = bioclimdata@crs)
#extract worldclim data with point lat long
je_bc_values <- raster::extract(bioclimdata,je_points)
je_lc_values <- raster::extract(landcov,je_points)
je_metadata = cbind(je_metadata, je_bc_values, je_lc_values)

#for pmi
pmi_meta = read.csv('~/Projects/MOVE/anopheles_pmi_cnrfp/metadata/pmi_metadata.csv')
pmi_coords <- data.frame(x=pmi_meta$long, y=pmi_meta$lat)
pmi_points <- sp::SpatialPoints(pmi_coords, proj4string = bioclimdata@crs)
#extract worldclim data with point lat long, rescale bioclim vars (to avoid discrepancy in mean and standard deviation among variables and ensure that the variable units were comparable)
pmi_bc_values <- raster::extract(bioclimdata,pmi_points)
pmi_bc_values <- scale(pmi_bc_values, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
pmi_lc_values <- raster::extract(landcov,pmi_points)
pmi_meta = cbind(pmi_meta, pmi_bc_values, pmi_lc_values)

#####
#find out which variables to include
#eliminate autocorrelated variables
#####


pmi_cormat <- cor(pmi_meta[60:81])
corrplot::corrplot(pmi_cormat, method = 'shade', order = 'AOE', diag = FALSE)
je_cormat <- cor(je_metadata[29:50])
corrplot::corrplot(je_cormat, method = 'shade', order = 'AOE', diag = FALSE)

#now, get fst thetas, roh and fis data
#read in roh stats and merge

#for je data
m_je_fis <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/f/m_form.fullset.pop.AgamP4_3L_fis')
s_je_fis <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/f/s_form.fullset.pop.AgamP4_3L_fis')

je_fismeta <- rbind(cbind(je_metadata[je_metadata$Form == 'M',], m_je_fis), cbind(je_metadata[je_metadata$Form == 'S',], s_je_fis))

je_roh <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/roh/segment_in_ROH.txt')
je_roh$seq_id <- gsub('.srt.dp.rg.bam.rohan.summary.txt', '' , je_roh$V1)
colnames(je_roh) <- c('gub', 'roh_median', 'roh_lower', 'roh_upper', 'seqid')

je_fis_roh_meta <- left_join(je_fismeta, je_roh[,2:5], by = c('seq_id'= 'seq_id'))

je_thetas <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/theta/thetas.segment_in_ROH.txt')
je_thetas$seq_id <- gsub('.srt.dp.rg.bam.rohan.summary.txt', '' , je_thetas$V1)
colnames(je_thetas) <- c('gub', 'gub2', 'theta_median', 'theta_lower', 'theta_upper', 'seqid')
je_fis_thetas_roh_meta <- left_join(je_fis_roh_meta, je_thetas[,3:6], by = c('seq_id'= 'seq_id'))




#for pmi data
pmi_roh_stats <- fread('~/Projects/MOVE/anopheles_pmi_cnrfp/data/roh/combined_roh_stats.csv')
#get inbreeding coefficients and merge with metadata
pmi_fis <- fread('/Users/tristanpwdennis/Projects/MOVE/anopheles_pmi_cnrfp/data/f/ngsf//pmi_f.txt', col.names = 'fis')
pmi_metadata_with_fis <- cbind(pmi_fis, pmi_meta)
metadata_with_fis_roh_thetas <- merge(pmi_metadata_with_fis, pmi_roh_stats, by.x = 'cgr_num', by.y = 'sampleid')
pmi_meta$year_site <- tolower(paste0(pmi_meta$Pop, pmi_meta$collection.year))
#some roh calcs failed (not yet sure why) - remove them as it gubs the median calculation
pmi_metaprun <- na.omit(metadata_with_fis_roh_thetas[,c('Pop', 'collection.year','med_rohlen', 'medtheta', 'roh_med', 'fis', 'Treatment')])
#median ROH, theta and fis stats by site and tp
agg = aggregate(metaprun,
                by = list(metaprun$Pop, metaprun$collection.year, metaprun$Treatment),
                FUN = median)

#median ROH, theta and fis stats by site and tp
count = aggregate(metaprun[,c('Pop', 'collection.year', 'Treatment')],
                  by = list(metaprun$Pop, metaprun$collection.year, metaprun$Treatment),
                  FUN = length)[,4]

agg = cbind(agg, count)

#drop shit cols
agg <- subset(agg, select = - c(Pop, collection.year, Treatment))

#make siteyear id for fst join
agg$siteyear = tolower(paste0(agg$Group.1, agg$Group.2))
#get fsts
globalstats_bysitebytp_list = list.files('~/Projects/MOVE/anopheles_pmi_cnrfp/data/fst/pmi/bytp_bysite/', pattern = 'global', full.names = TRUE)
#create names by removing path
bytp_bysite_names = data.frame(do.call(rbind, strsplit(globalstats_bysitebytp_list, '/', fixed=TRUE)))
#bind metadata to fst and have a clean
bysite_by_tp_globfst = cbind(do.call(rbind,strsplit(bytp_bysite_names$X12, "_", fixed=TRUE)), do.call(rbind,lapply(globalstats_bysitebytp_list, fread)))
bysite_by_tp_globfst$V5 = gsub('globalstat.txt','',bysite_by_tp_globfst$V5)
bysite_by_tp_globfst <- bysite_by_tp_globfst[,c(1,2, 5, 6, 7)]
colnames(bysite_by_tp_globfst) <- c('site', 'date', 'chrom', 'unweighted_fst', 'weighted_fst')
bysite_by_tp_globfst$siteyear = paste0(bysite_by_tp_globfst$site, bysite_by_tp_globfst$date)
roh_fst_thetastats<- na.omit(merge(agg, bysite_by_tp_globfst, by = 'siteyear')) #join and remove gubbins





#per study PC data
pmi_pca = eigen(read.table('/Users/tristanpwdennis/Projects/MOVE/anopheles_pmi_cnrfp/data/admix_pca/AgamP4_3L_noqfilter.pcangsd.cov'))
je_pca = eigen(read.table('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/structure/3_euch_pruned.cov'))

pmi_df = cbind(data.frame(pmi_pca$vectors), pmi_meta)
je_df = cbind(data.frame(je_pca$vectors), je_metadata)

#make polygon surrounding je samples (for eems and other analyses)
je_hull = sf::st_convex_hull(sf::st_union(sf::st_as_sf(je_points)))
plot(je_hull)
#for all pmi samples
pmi_polygon = sf::st_convex_hull(sf::st_union(sf::st_as_sf(pmi_points)))
plot(pmi_polygon)
make_pmi_tp_polygons <- function(year) {
  yearmeta <- pmi_meta[pmi_meta$collection.year == year,]
  coords <- data.frame(x=yearmeta$long, y=yearmeta$lat)
  points <- sp::SpatialPoints(coords, proj4string = r@crs)
  hull = sf::st_convex_hull(sf::st_union(sf::st_as_sf(points)))
  return(hull)
  }

years = as.list(unique(pmi_meta$collection.year))
lapply(years, make_pmi_tp_polygons)

#do pmi rda with 100k mafs
#read in mafs, subset to manageable size, remove unwanted cols, transpose and add row and colnames
pmi_mafs <- fread('/Users/tristanpwdennis/Projects/MOVE/anopheles_pmi_cnrfp/data/mafs/pmi_mafs_by_site.csv.gz')
colnames(pmi_mafs) <- c("num","chrom","pos","snpid","binduli_2008","binduli_2020","bunbuna_2014","bunbuna_2018","bunbuna_2020","dimabi_2014","gbullung_2014","gbullung_2018","kulaa_2014","kulaa_2016","kulaa_2020","nanton_2008","nanton_2014","nanton_2018","tarikpaa_2008","tarikpaa_2016","tarikpaa_2018","tugu_2020","woribugu_2014","woribugu_2016")
pmi_mafs <- pmi_mafs[complete.cases(pmi_mafs[ , 4:24]),]
pmi_mafs <- select(pmi_mafs, -num, -chrom, -pos, -snpid)

v = as.matrix(transpose(pmi_mafs))
colnames(v) <- rownames(pmi_mafs)
rownames(v) <- colnames(pmi_mafs)
#make covariate data
pmi_df$site_tp <- paste0(pmi_df$Pop, '_', pmi_df$collection.year) #make id for comparing meta to mafs
pmi_sites <-unique(select(pmi_df, site_tp, Treatment, lat, long, bio1:bio19, pmi_lc_values)) #get distinct site values from pmi meta
pmi_sites$site_tp <- tolower(pmi_sites$site_tp)
pmi_sites <- pmi_sites[pmi_sites$site_tp %in% rownames(v),] #subset to rows in maf dataset (compare to maf rownames eg sites)
#oprder both by site and tp
pmi_sites = pmi_sites[ order(pmi_sites$site_tp), ]
v = v[ order(row.names(v)), ]
#v <- v[,colSums(is.na(v))<nrow(v)] #remove all na cols (sites)

#null model
RDA0 <- rda(v ~ 1,  pmi_sites) 
#full model rda with all our variables
RDAfull <- rda(v ~  Treatment + bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19 + pmi_lc_values, pmi_sites)
## Stepwise procedure with ordiR2step function
mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
summary(mod)

je_mafs <- fread('/Users/tristanpwdennis/Projects/MOVE/bulky_larvae_paper_data/mafs/je_mafs_by_site.csv.gz')
head(je_mafs)

ibs_pmi <- fread('/Users/tristanpwdennis/Projects/MOVE/anopheles_pmi_cnrfp/data/ibs/AgamP4_3L_noqfilter.ibsMat')
ibs_pmi <- ibs_pmi[,1:219]
eigen <- eigen(ibs_pmi)
plot(eigen$vectors[,6], eigen$vectors[,7])

agc.adonis <- adonis2(ibs_pmi ~ Pop, data=pmi_meta, permutations=9999)
agc.adonis

?vegan::adonis2

#per study and TP sample polygons
#for eems and rda based analysis