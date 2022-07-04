#####
#Analysis of relatedness, AF and Fst data from PMI per-samples per-tps
#####

#load our shit
#run this in the same directory 
pkg = c('tidyverse', 'geosphere', 'data.table', 'vegan' ,'cowplot', 'parallel', 'Rcpp', 'RColorBrewer', 'UpSetR', 'lme4', 'moments', 'parallel')
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

#soil moisture (from https://earlywarning.usgs.gov/fews/product/311)
sm <- raster::raster('~/Downloads/wa_monthly_fldas_soilmoi10_40cm_tavg_1405/wa_monthly_fldas_soilmoi10_40cm_tavg_1405.tif')

#load bioclim data
#get bioclimatic vars from worldclim
bioclimdata <- raster::getData("worldclim",var="bio",res=10)
bioclimdata <- bioclimdata[[1:19]]
#get landcover data
landcov = raster::raster('~/Downloads/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif') #load raster file
#get rice data
rice <- raster::raster('~/Downloads/spam2017v2r1_ssa_harv_area.geotiff/spam2017V2r1_SSA_H_RICE_A.tif') 

#for pmi and je metadata, use raster library to extract info from above data for each site

#for je
je_metadata=read.csv('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/metadata/metadata_full_sequencedonly_062022.csv')
#get coords from metadata and make into raster spatialpoints obj
je_coords <- data.frame(x=je_metadata$Long, y=je_metadata$Lat)
je_points <- sp::SpatialPoints(je_coords, proj4string = bioclimdata@crs)
#extract worldclim data with point lat long
je_bc_values <- raster::extract(bioclimdata,je_points)
je_lc_values <- raster::extract(landcov,je_points)
je_sm_values <- raster::extract(sm,je_points)
je_rice <- raster::extract(rice, je_points)
je_metadata = cbind(je_metadata, je_bc_values, je_sm_values, je_rice, je_lc_values)

#for pmi
pmi_meta = read.csv('~/Projects/MOVE/anopheles_pmi_cnrfp/metadata/pmi_metadata.csv')
pmi_coords <- data.frame(x=pmi_meta$long, y=pmi_meta$lat)
pmi_points <- sp::SpatialPoints(pmi_coords, proj4string = bioclimdata@crs)
#extract worldclim data with point lat long, rescale bioclim vars (to avoid discrepancy in mean and standard deviation among variables and ensure that the variable units were comparable)
pmi_bc_values <- raster::extract(bioclimdata,pmi_points)
pmi_bc_values <- scale(pmi_bc_values, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
pmi_lc_values <- raster::extract(landcov,pmi_points)
pmi_sm_values <- raster::extract(sm,pmi_points)
pmi_rice <- raster::extract(rice, pmi_points)

plot(pmi_sm_values)
pmi_meta = cbind(pmi_meta, pmi_bc_values, pmi_sm_values, pmi_rice, pmi_lc_values)

#now, get individual thetas, roh and fis data
#read in roh stats and merge

#for je data
m_je_fis <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/f/m_form.fullset.pop.AgamP4_3L_fis')
s_je_fis <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/f/s_form.fullset.pop.AgamP4_3L_fis')

je_fismeta <- rbind(cbind(je_metadata[je_metadata$Form == 'M',], m_je_fis), cbind(je_metadata[je_metadata$Form == 'S',], s_je_fis))

je_roh <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/roh/segment_in_ROH.txt')
je_roh$seq_id <- gsub('.srt.dp.rg.bam.rohan.summary.txt', '' , je_roh$V1)
colnames(je_roh) <- c('gub', 'roh_median', 'roh_lower', 'roh_upper', 'seq_id')

je_fis_roh_meta <- left_join(je_fismeta, je_roh[,2:5], by = c('seq_id'= 'seq_id'))

je_thetas <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/theta/thetas.segment_in_ROH.txt')
je_thetas$seq_id <- gsub('.srt.dp.rg.bam.rohan.summary.txt', '' , je_thetas$V1)
colnames(je_thetas) <- c('gub', 'gub2', 'theta_median', 'theta_lower', 'theta_upper', 'seq_id')
je_fis_thetas_roh_meta <- left_join(je_fis_roh_meta, je_thetas[,3:6], by = c('seq_id'= 'seq_id'))

#for pmi data
pmi_roh_stats <- fread('~/Projects/MOVE/anopheles_pmi_cnrfp/data/roh/combined_roh_stats.csv')
#get inbreeding coefficients and merge with metadata
pmi_fis <- fread('/Users/tristanpwdennis/Projects/MOVE/anopheles_pmi_cnrfp/data/f/ngsf//pmi_f.txt', col.names = 'fis')
pmi_metadata_with_fis <- cbind(pmi_fis, pmi_meta)
metadata_with_fis_roh_thetas <- merge(pmi_metadata_with_fis, pmi_roh_stats, by.x = 'cgr_num', by.y = 'sampleid')
pmi_meta$year_site <- tolower(paste0(pmi_meta$Pop, pmi_meta$collection.year))
#some roh calcs failed (not yet sure why) - remove them as it gubs the median calculation




#####
#find out which variables to include
#eliminate autocorrelated variables
#####


pmi_cormat <- cor(unique(pmi_meta[60:82]))
corrplot::corrplot(pmi_cormat, method = 'shade', order = 'AOE', diag = FALSE)
#let's keep land cover, bio12 and bio10

je_cormat <- cor(na.omit(je_metadata[29:51]))
corrplot::corrplot(je_cormat, method = 'shade', order = 'AOE', diag = FALSE)
#let's keep land cover, bio12, bio1, bio5, soil moisture and


plot(je_metadata$je_rice, je_metadata$je_sm_values)

x = glmer(V1 ~  bio12 + bio1 + bio5 + je_rice + (1|Site) + (1|Form) + Long, data = je_fis_thetas_roh_meta, family = 'binomial')
x = glmer(fis ~  bio12 + bio1 + bio5 + pmi_rice + (1|Pop) + long, data = metadata_with_fis_roh_thetas, family = 'binomial')
summary(x)





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





#for eems and rda based analysis