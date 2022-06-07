####IBD segment analysis
####Tristan Dennis 10.03.22

#load env and metadata
pkg = c("tidyverse", "gridExtra", "RColorBrewer", "moments", "viridis" , "data.table", "itertools", "geosphere", 'UpSetR', 'sparseAHC', 'igraph', 'adegenet', 'vegan')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#load metadata and clean
metadata = read.csv('~/Projects//MOVE/data/anopheles_03_21_DW/metadata/location_ecozone_metadata.csv')
oldmeta = read.csv('~/Projects//MOVE/td_je_angam_2022/metadata/metadata_with_insecticide_info.csv')
newmeta = left_join(metadata, bamlist, by=c('sample' = 'sample_id'))
metastrip = metadata %>% select(Lat, Long,Site)
fullstrip = oldmeta %>% select(-Lat, -long)
metadata = left_join(fullstrip, metastrip, by = 'Site') %>% unique()
#####################
#Analysis of regular segments
#do isobd analysis, and calculate seglen over genes of interest and determine whether they are longer than average
#####################

setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_full_cov/')
regseg_list = list.files(path='/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_full_cov/', pattern = '*.ibd', recursive=TRUE)
cnam=c('ind1', 'hapindex1', 'ind2', 'hapindex2', 'chr', 'start', 'end', 'lod')
ibd_df = rbindlist(lapply(seq_along(regseg_list), function(i){fread(regseg_list[[i]], col.names = cnam) %>% mutate(fn=regseg_list[[i]])})) #fread all, adding filename as col
#calculate per seg length
ibd_df$len= ibd_df$end - ibd_df$start #find out the length
summary_ibd_table$dist[is.na(summary_ibd_table$dist)] <- 0
#summary of total and median seg length shared by individuals
summary_ibd_table = setDT(ibd_df)[,c('ind1', 'ind2', 'len')][ , .(total_seglen = sum(len), median_seglen = median(len), count_segs = .N), by = .(ind1, ind2)]

#gene coordinate table
gois = data.table(gene=c('kdr', 'rdl', 'coeae', 'ace1', 'cyp6', 'gst', 'cyp6m2.z1', 'cyp9k1'),
                  chr = c('2L', '2L', '2L', '2R', '2R', '3R', '3R', 'X'),
                  start =c(2358158, 25363652, 28548433, 3484107, 28463000, 28580000, 6900000, 15222000),
                  end=c(2431617, 25434556, 28550748, 3495790, 28568000, 28605000, 7030000, 15257000))


#join metadata to each side of the ibd dataframe
summary_ibd_table = left_join(summary_ibd_table, metadata[,c('phaseid', 'Lat', 'Long', 'Site', 'Form', 'habitat')], by=c('ind1' = 'phaseid'))%>% 
  left_join(.,metadata[,c('phaseid', 'Lat', 'Long', 'Site', 'Form', 'habitat')], by=c('ind2' = 'phaseid'))



#rbind(df, setNames(c('ind0', 'ind0', 0,0,0), names(df)))


#calculate geographic distance between pairs
summary_ibd_table$dist = distGeo(summary_ibd_table[,c('Lat.x','Long.x')], summary_ibd_table[,c('Lat.y','Long.y')])



#create within or between species comparison
#create coimparison factor levels
summary_ibd_table$formcomp = as.factor(paste0(summary_ibd_table$Form.x, summary_ibd_table$Form.y))
summary_ibd_table$ecocomp = as.factor(paste0(summary_ibd_table$habitat.x, summary_ibd_table$habitat.y))
summary_ibd_table
#get rid of samesite comps
diff_sites_ibd_table = summary_ibd_table[dist >0]

dsamp = summary_ibd_table[sample(nrow(summary_ibd_table), 1000), ]

ggplot(diff_sites_ibd_table, aes(x=))
max(diff_sites_ibd_table$total_seglen)

#subset to mform:mform comparisons
mm=diff_sites_ibd_table[formcomp == 'MM']
ss=diff_sites_ibd_table[formcomp == 'SS']
forminds = ss

xy <- as.data.frame(t(combn(filter(metadata, Form == 'S')$ind_id, 2)), ) #create combinations of all m form seq ids)
colnames(xy) = c('ida', 'idb')
xy$pairid = paste0(xy$ida,xy$idb) #paste to make unique pair id
length(unique(xy$ida)) #are they the same
length(unique(xy$idb))
forminds$pairid = paste0(forminds$ind1, forminds$ind2) #make pairid from ibd seg data
xy = left_join(xy, forminds, by = c('pairid' = 'pairid'))
xy$total_seglen[is.na(xy$total_seglen)] <- 0
xy$dist[is.na(xy$dist)] <- 0

gd = xy %>% dplyr::select(ida, idb, total_seglen) %>% pivot_wider(., names_from = ida, values_from = total_seglen ) 
gd = rbind(c('ind1',rep(0, ncol(gd))), gd)
gd = cbind(gd, rep(0, nrow(gd)))
gd = gd[,-1] #remove row ids
gd = as.dist(gd, upper = FALSE) #take lower triangle to dist matrix object


ed = xy %>% dplyr::select(ida, idb, dist) %>% pivot_wider(., names_from = ida, values_from = dist) 
ed = rbind(c('ind1',rep(0, ncol(ed))), ed)
ed = cbind(ed, rep(0, nrow(ed)))
ed = ed[,-1] #remove row ids
ed = as.dist(ed, upper = FALSE) #take lower triangle to dist matrix object
ed

ibd = mantel.randtest(ed, gd)
plot(ibd)
ibd
g = mantel.correlog(D.eco = ed, D.geo = gd, n.class = 200)
plot(g)



############
#with kinship
relframe = read.delim('~/Projects/MOVE/td_je_angam_2022/data/wg_res.res', header = F)
#specify column names 
colnames(relframe) = c('a',  'b',  'nSites',  'J9', 'J8', 'J7', 'J6', 'J5', 'J4', 'J3', 'J2', 'J1', 'rab', 'Fa', 'Fb', 'theta', 'inbred_relatedness_1_2', 'inbred_relatedness_2_1', 'fraternity', 'identity',  'zygosity',  '2of3IDB', 'FDiff', 'loglh', 'nIter','notsure', 'coverage',  '2dsfs', 'R0', 'R1', 'KING', '2dsfs_loglike', '2dsfsf_niter')
#join both ind sides with metadata
relframe = left_join(relframe, metadata, by = c('a' = 'bamorder')) %>% left_join(., metadata, by = c('b' = 'bamorder'))
#create within or between form level and get rid of between form comparisons
#get distances between sampling points
relframe$dist = distGeo(relframe[,c('Lat.x','Long.x')], relframe[,c('Lat.y','Long.y')])
relframe$formcomp = paste0(relframe$Form.x, relframe$Form.y)
relframe$dist[is.na(relframe$dist)] <- 0
#include only m:m comparisons (as we expect almost total reproductive isiolation)
res = relframe %>% filter(Form.x == 'S' & Form.y =='S')
#res[is.na(res$dist)] <- 0 #get rid of annoying na = 0
res = res[res$dist > 0,]
res$dist[is.na(res$dist)] <- 0
res$rab[is.na(res$rab)] <- 0
res$ida = paste0('ind', res$a)
res$idb = paste0('ind', res$b)
res$pairid = paste0(res$ida, res$idb)
xy <- as.data.frame(t(combn(filter(metadata, Form == 'S')$ind_id, 2)), ) #create combinations of all m form seq ids)
colnames(xy) = c('ida', 'idb')
xy$pairid = paste0(xy$ida,xy$idb) #paste to make unique pair id
xy = left_join(xy, res, by = c('pairid' = 'pairid'))
xy$dist[is.na(xy$dist)] <- 0

gd = xy %>% dplyr::select(ida.x, idb.x, rab) %>%  pivot_wider(., names_from = ida.x, values_from = rab) 
gd = rbind(c('ind1',rep(0, ncol(gd))), gd)
gd = cbind(gd, rep(0, nrow(gd)))
gd = gd[,-1] #remove row ids
gd = as.dist(gd, upper = FALSE) #take lower triangle to dist matrix object
gd[is.na(gd)] <- 0

ed = xy %>% dplyr::select(ida.x, idb.x, dist) %>% pivot_wider(., names_from = ida.x, values_from = dist) 
ed = rbind(c('ind1',rep(0, ncol(ed))), ed)
ed = cbind(ed, rep(0, nrow(ed)))
ed = ed[,-1] #remove row ids
ed = as.dist(ed, upper = FALSE) #take lower triangle to dist matrix object
ed[is.na(ed)] <- 0


ed

ibd = mantel.randtest(ed, gd)
ibd
plot(ibd)


g = mantel.correlog(D.eco = ed, D.geo = gd, n.class = 500)
plot(g)

#########################################################
#PLOT ISOLATION BY DISTANCE WITH KINSHIP AND IBD SEGS   #
#########################################################

#random effects
#try within or between species
#within or between ecozone

#try to visualise point density
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#get point density
summary_ibd_table$density <- get_density(summary_ibd_table$dist, summary_ibd_table$total_seglen, n = 1000)
relframe$density <- get_density(relframe$dist, relframe$rab, n = 1000)
relframe$dist
#plot seglen/dist with density as colour
ibdseg = ggplot(filter(summary_ibd_table, (formcomp == 'MM' | formcomp == 'SS') & dist > 0 ), aes(x=dist,y=total_seglen, color=density))+
  facet_wrap(~formcomp)+
  geom_point()+
  scale_colour_viridis()+
  theme_classic()+  labs(x='Geographic Distance (m)', y='Genome Shared In IBD Segments (bp)')+
  theme(legend.position =  "none")


#plot seglen/dist with density as colour
kin = ggplot(filter(relframe, ( formcomp == 'MM' | formcomp == 'SS') & dist > 0 ), aes(x=dist,y=rab, color=density))+
  facet_wrap(~formcomp)+
  geom_point()+
  scale_colour_viridis()+
  theme_classic()+
  labs(x='Geographic Distance (m)', y='Kinship Coefficient (rab)')+
  theme(legend.position =  "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank())


cowplot::plot_grid(kin, ibdseg, align = 'v', ncol=1)


.###take segs from gene of interest regions and do a pca with them/clustering


#####################
#Analysis of sitewise Fst
#####################

fst_site = fread('~/Projects/MOVE/td_je_angam_2022/data/fst_betweensites.csv')
left_join(fst_site, metadata, by=c('site_a' = 'Site')) %>% left_join(., metadata, by = c('site_b' = 'Site'))






#####################
#Analysis of downsampled segments
#Take each downsampled file (for chr3), load, bind, calculate summary stats  and isobd
#####################
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_chr_downsampled/')
dsample_list = list.files(path='/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_chr_downsampled/', pattern = '*.ibd', recursive=TRUE)
#colname vector
cnam=c('ind1', 'hapindex1', 'ind2', 'hapindex2', 'chr', 'start', 'end', 'lod')
#read all files into big df
ibd_df = rbindlist(lapply(seq_along(dsample_list), function(i){fread(dsample_list[[i]], col.names = cnam) %>% mutate(fn=dsample_list[[i]])})) #fread all, adding filename as col
ibd_df$fn = gsub(ibd_df$fn, pattern = '_form.core.', replacement = ',') #remove some rubbish
ibd_df$fn = gsub(ibd_df$fn, pattern = '.list.', replacement = ',')
ibd_df$fn = gsub(ibd_df$fn, pattern = '_phasedb4.0.ibd', replacement = '') 
ibd_df = separate(ibd_df, fn, sep=',', into = c('form', 'coverage', 'chrom')) #separate filename into grouping vars
ibd_df$len= ibd_df$end - ibd_df$start #find out the length
ibd_df$coverage_factor <- factor(ibd_df$coverage, levels = c("full", "10", "7.5", "5", "2.5", "1", "0.5"))

aggregate(ibd_df$len, by = list(ibd_df$chr, ibd_df$coverage), FUN = "median")

#calculate summary stats for downsampled ibd segs for length
ibd_df[ , .(skewness(len), as.numeric(median(len)), .N), by = .(chr, coverage)]
ibd_df[ , .(skewness(lod), as.numeric(median(lod)), .N), by = .(chr, coverage)]

#plot lod
ibd_df %>% filter(form=='m' & chr == 'AgamP4_3L') %>% 
  ggplot(aes(x=len, y=lod))+
  facet_wrap(~coverage_factor)+
  geom_point(size=1, alpha=0.1, width = .2)+
  theme_classic()

ibd_df %>% filter(form=='m' & chr == 'AgamP4_3L') %>% 
  ggplot(aes(x=coverage_factor, y=len))+
  geom_jitter(alpha=0.3)+
  geom_boxplot()+
  theme_classic()


