####IBD segment analysis
####Tristan Dennis 10.03.22

#load env and metadata
pkg = c("tidyverse", "gridExtra", "RColorBrewer", "moments", "viridis" , "data.table", "itertools", "geosphere", 'UpSetR', 'ggthemes', 'devtools', 'adegenet', 'vegan')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
source_url("https://raw.githubusercontent.com/sjmurdoch/fancyaxis/master/fancyaxis.R")
#load metadata and clean
metadata = read.csv('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/metadata/metadata_full_sequencedonly_062022.csv')

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
#summary of total and median seg length shared by individuals
summary_ibd_table = setDT(ibd_df)[,c('ind1', 'ind2', 'len')][ , .(total_seglen = sum(len), median_seglen = median(len), count_segs = .N), by = .(ind1, ind2)]
summary_ibd_table$dist[is.na(summary_ibd_table$dist)] <- 0

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
summary_ibd_table$prop_seg = summary_ibd_table$total_seglen / 265000000
#get rid of samesite comps
diff_sites_ibd_table = summary_ibd_table[dist >0]


form='S'
make_seg_matrices <- function(form) {
forminds = diff_sites_ibd_table[diff_sites_ibd_table$Form.x == form & diff_sites_ibd_table$Form.y == form]
xy <- as.data.frame(t(combn(filter(metadata, Form == form)$ind_id, 2)), ) #create combinations of all m form seq ids)
colnames(xy) = c('ida', 'idb')
xy$pairid = paste0(xy$ida,xy$idb) #paste to make unique pair id
length(unique(xy$ida)) #are they the same
length(unique(xy$idb))
forminds$pairid = paste0(forminds$ind1, forminds$ind2) #make pairid from ibd seg data
xy = left_join(xy, forminds, by = c('pairid' = 'pairid'))
xy$total_seglen[is.na(xy$total_seglen)] <- 0
xy$dist[is.na(xy$dist)] <- 0

gd = xy %>% dplyr::select(ida, idb, prop_seg) %>% pivot_wider(., names_from = ida, values_from = prop_seg ) 
gd = rbind(c('ind1',rep(0, ncol(gd))), gd)
gd = cbind(gd, rep(0, nrow(gd)))
gd = gd[,-1] #remove row ids
gd = as.dist(gd, upper = FALSE) #take lower triangle to dist matrix object


ed = xy %>% dplyr::select(ida, idb, dist) %>% pivot_wider(., names_from = ida, values_from = dist) 
ed = rbind(c('ind1',rep(0, ncol(ed))), ed)
ed = cbind(ed, rep(0, nrow(ed)))
ed = ed[,-1] #remove row ids
ed = as.dist(ed, upper = FALSE) #take lower triangle to dist matrix object
return(list(gd, ed))
}

m_seg_matrices = make_seg_matrices('M')
s_seg_matrices = make_seg_matrices('S')
plot(m_seg_matrices[1], m_seg_matrices[[2]])

plot(ed,gd)

m_seg_mantel = vegan::mantel(m_seg_matrices[[1]], m_seg_matrices[[2]], permutations = 999, na.rm = TRUE)
s_seg_mantel = vegan::mantel(s_seg_matrices[[1]], s_seg_matrices[[2]], permutations = 999, na.rm = TRUE)
mantel(ed, gd, na.rm = TRUE)

#plot mantel test output
tiff("~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/figures/sfig1_fst_isobd.tiff", width = 20, height = 25, units = 'cm', res = 800)
par(mfrow =c(3,2))
plot(as.vector(m_seg_matrices[[1]]), as.vector(m_seg_matrices[[2]]), ylab = 'Fst / (1 - Fst)', xlab = 'Geographic Distance (m)', col = '#1b9e77', pch = 16,mar = c(1.5,10,1.5,1.5), main="Genetic vs Geographic Distance for A. coluzzi")
plot(as.vector(s_seg_matrices[[1]]), as.vector(s_seg_matrices[[2]]), ylab = 'Fst / (1 - Fst)', xlab = 'Geographic Distance (m)', col = '#d95f02', pch = 16, main="Genetic vs Geographic Distance for A. gambiae")
m_corr = mantel.correlog(D.geo = m_matrices[[1]], D.eco = m_matrices[[2]], n.class = 20, mult = 'bonferroni')
s_corr = mantel.correlog(D.geo = s_matrices[[1]], D.eco = s_matrices[[2]], n.class = 20, mult = 'bonferroni')
plot(m_corr, main = 'Mantel Correlelogram  for A. coluzzi')
plot(s_corr, main = 'Mantel Correlelogram  for A. gambiae')
m_corr
dev.off()





############
#with kinship
relframe = read.delim('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/relatedness/wg_res.res', header = F)
#specify column names 
colnames(relframe) = c('a',  'b',  'nSites',  'J9', 'J8', 'J7', 'J6', 'J5', 'J4', 'J3', 'J2', 'J1', 'rab', 'Fa', 'Fb', 'theta', 'inbred_relatedness_1_2', 'inbred_relatedness_2_1', 'fraternity', 'identity',  'zygosity',  '2of3IDB', 'FDiff', 'loglh', 'nIter','notsure', 'coverage',  '2dsfs', 'R0', 'R1', 'KING', '2dsfs_loglike', '2dsfsf_niter')
#join both ind sides with metadata
relframe = left_join(relframe, metadata, by = c('a' = 'bamorder')) %>% left_join(., metadata, by = c('b' = 'bamorder'))
#sort
relframe <- relframe[
  order( relframe[,1], relframe[,2] ),
]

#create within or between form level and get rid of between form comparisons
#get distances between sampling points
relframe$dist = distGeo(relframe[,c('Lat.x','Long.x')], relframe[,c('Lat.y','Long.y')])
relframe$dist[is.na(relframe$dist)] <- 0

relframe$rab[relframe$dist == 0] <- NA
relframe$dist[relframe$dist == 0] <- NA

#res = relframe %>% filter(Form.x == 'M' & Form.y =='M')
#rabmat = pivot_wider(res[,c('a', 'b', 'rab')], names_from = a, values_from = rab)
#rabmat = rabmat[,-1] #remove row ids
#rabmat = as.dist(rabmat, upper = FALSE) #take lower triangle to dist matrix object
#rabmat[order(as.numeric(colnames(rabmat)))]
#plot isobd
relframe = mutate(relframe, spec_comp = case_when(
  Form.x == 'M' & Form.y == 'M' ~ 'A. coluzzi/A. coluzzi',
  Form.x == 'S' & Form.y == 'S' ~ 'A. gambiae/A. gambiae',
  Form.x == 'M' & Form.y == 'S' ~ 'A. coluzzi/A. gambiae',
  Form.x == 'S' & Form.y == 'M' ~ 'A. gambiae/A. coluzzi'
))


make_relmat_fromibs <- function(complevel){
#include only m:m comparisons (as we expect almost total reproductive isiolation)
 # complevel = 'A. coluzzi/A. coluzzi'
res = relframe[relframe$spec_comp == complevel,]
#res[is.na(res$dist)] <- 0 #get rid of annoying na = 0
#res = relframe[relframe$dist > 0,]
res$ida = paste0('ind', res$a)
res$idb = paste0('ind', res$b)
gd = res %>% dplyr::select(a, b, rab) %>%  pivot_wider(., names_from = a, values_from = rab) 
gd = rbind(c('ind1',rep(0, ncol(gd))), gd)
gd = cbind(gd, rep(0, nrow(gd)))
gd = gd[,-1] #remove row ids
gd = as.dist(gd, upper = FALSE) #take lower triangle to dist matrix object

ed = res %>% dplyr::select(a, b, dist) %>% pivot_wider(., names_from = a, values_from = dist) 
ed = rbind(c('ind1',rep(0, ncol(ed))), ed)
ed = cbind(ed, rep(0, nrow(ed)))
ed = ed[,-1] #remove row ids
ed = as.dist(ed, upper = FALSE) #take lower triangle to dist matrix object

return(list(ed, gd))
}



ggplot(relframe[relframe$dist > 0,],aes(x=dist, y=rab, color=spec_comp))+
  facet_wrap(~spec_comp)+
  geom_point()+
  geom_rug()+ 
  theme_classic()+
  scale_color_brewer(palette = 'Dark2')


m_kinmats <- make_relmat_fromibs('A. coluzzi/A. coluzzi')
s_kinmats <- make_relmat_fromibs('A. gambiae/A. gambiae')
m_ibd = vegan::mantel(m_kinmats[[1]], m_kinmats[[2]], permutations = 999, na.rm = TRUE)
s_ibd = vegan::mantel(s_kinmats[[1]], s_kinmats[[2]], permutations = 999, na.rm = TRUE)
plot(s_kinmats[[1]], s_kinmats[[2]])


m = mantel.correlog(D.eco = m_kinmats[[1]], D.geo = m_kinmats[[2]], n.class = 12)
warnings()


s = mantel.correlog(D.eco = s_kinmats[[1]], D.geo = s_kinmats[[2]], n.class = 500)
plot(m)
plot(m)


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
#Kinship matrix from pcrelate
#####################

m_meta <- metadata[metadata$Form == 'M',][c('ind_id','Lat', 'Long')]
s_meta <- metadata[metadata$Form == 'S',][c('ind_id','Lat', 'Long')]
m_kinmat = RcppCNPy::npyLoad('/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/kinship/mform_3L.pcangsdrelate.kinship.npy')
s_kinmat = RcppCNPy::npyLoad('/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/kinship/sform_3L.pcangsdrelate.kinship.npy')
m_xy <- as.data.frame(t(combn(m_meta$ind_id, 2)), ) 
s_xy <- as.data.frame(t(combn(s_meta$ind_id, 2)), ) 

m_geodat = merge(m_meta,merge(m_meta,m_xy,  by.x = 'ind_id', by.y = 'V1'),  by.x = 'ind_id', by.y = 'V2')
s_geodat = merge(s_meta,merge(s_meta,s_xy,  by.x = 'ind_id', by.y = 'V1'),  by.x = 'ind_id', by.y = 'V2')

m_geodat$dist = distGeo(m_geodat[,c('Lat.x','Long.x')], m_geodat[,c('Lat.y','Long.y')])
s_geodat$dist = distGeo(s_geodat[,c('Lat.x','Long.x')], s_geodat[,c('Lat.y','Long.y')])
sort(m_geodat)
pivot_wider(m_geodat[c('ind_id', 'ind_id.y', 'dist')], names_from = ind_id, values_from = 'dist')




#####################
#Analysis of sitewise Fst
#####################
global_fst <- fread("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/fst/betweensites/newfst_global.txt", col.names = c('site_form_a', 'site_form_b', 'unweighted', 'weighted'), header = T)
global_fst <- separate(global_fst, col = site_form_a, into = c('site_a', 'form_a'), sep = '_', remove = F)
global_fst <- separate(global_fst, col = site_form_b, into = c('site_b', 'form_b'), sep = '_', remove = F)
global_fst$weighted = global_fst$weighted / (1 - global_fst$weighted)
make_matrices_fst <- function(species) {
matfor= species
siteinfo = select(metadata, Site, Lat, Long) %>% unique()
fst_premat <- global_fst[global_fst$form_a == matfor & global_fst$form_b == matfor]
fst_premat$site_a <- gsub(' ','-',fst_premat$site_a)
fst_premat$site_b <- gsub(' ','-',fst_premat$site_b)
siteinfo$Site <- gsub(' ','-',siteinfo$Site )
fst_premat <- fst_premat %>% left_join(., siteinfo, by = c('site_a' = 'Site')) %>% left_join(., siteinfo, by = c('site_b' = 'Site'))
fst_premat$dist = distGeo(fst_premat[,c('Lat.x','Long.x')], fst_premat[,c('Lat.y','Long.y')])
fst_premat = fst_premat[dist > 0]

fd = fst_premat %>% dplyr::select(site_a, site_b, weighted) %>%  pivot_wider(., names_from = site_a, values_from = weighted) 
fd = fd[,-1] #remove row ids
fd = as.dist(fd, upper = FALSE) #take lower triangle to dist matrix object

#make geodist mat
siteinfo = unique(metadata[c('Site', 'Lat', 'Long')])

gd = fst_premat %>% dplyr::select(site_a, site_b, dist) %>%  pivot_wider(., names_from = site_a, values_from = dist) 
gd = gd[,-1] #remove row ids
gd = as.dist(gd, upper = FALSE) #take lower triangle to dist matrix object
return(list(gd, fd, fst_premat))
}
?distGeo
m_matrices = make_matrices_fst('M')
s_matrices = make_matrices_fst('S')
m_mantel <- mantel.randtest(m1 = m_matrices[[2]], m2 = m_matrices[[1]], nrepet = 999)
s_mantel <- mantel.randtest(m1 = s_matrices[[2]], m2 = s_matrices[[1]], nrepet = 999)

#plot mantel test output
tiff("~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/figures/sfig1_fst_isobd.tiff", width = 20, height = 25, units = 'cm', res = 800)
par(mfrow =c(3,2))
plot(as.vector(m_matrices[[1]]), as.vector(m_matrices[[2]]), ylab = 'Fst / (1 - Fst)', xlab = 'Geographic Distance (m)', col = '#1b9e77', pch = 16,mar = c(1.5,10,1.5,1.5), main="Genetic vs Geographic Distance for A. coluzzi")
plot(as.vector(s_matrices[[1]]), as.vector(s_matrices[[2]]), ylab = 'Fst / (1 - Fst)', xlab = 'Geographic Distance (m)', col = '#d95f02', pch = 16, main="Genetic vs Geographic Distance for A. gambiae")
plot(m_mantel, main = 'Mantel Test Histogram for A. coluzzi')
plot(s_mantel, main = 'Mantel Test Histogram for A. gambiae')
m_corr = mantel.correlog(D.geo = m_matrices[[1]], D.eco = m_matrices[[2]], n.class = 20, mult = 'bonferroni')
s_corr = mantel.correlog(D.geo = s_matrices[[1]], D.eco = s_matrices[[2]], n.class = 20, mult = 'bonferroni')
plot(m_corr, main = 'Mantel Correlelogram  for A. coluzzi')
plot(s_corr, main = 'Mantel Correlelogram  for A. gambiae')
dev.off()
s_corr
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

####using ibs mat


m_ibsmat <- as.matrix(fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/ibs/m_form.fullset.pop.AgamP4_3L.ibsMat'))
m_ibsmat <- m_ibsmat[,-160]
m_ibsmat = as.dist(m_ibsmat, upper = FALSE) #take lower triangle to dist matrix object

s_ibsmat <- fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/data/ibs/s_form.fullset.pop.AgamP4_3L.ibsMat')

  

m_xy <- as.data.frame(t(combn(metadata[metadata$Form == 'M',]$seq_id, 2)), ) 
m_xy = left_join(m_xy, metadata, by = c('V1' = 'seq_id')) %>% left_join(., metadata, by = c('V2' = 'seq_id')) %>% select(V1, V2, Lat.x, Long.x, Lat.y, Long.y)
m_xy$dist = distGeo(m_xy[,c('Lat.x','Long.x')], m_xy[,c('Lat.y','Long.y')])

gd = m_xy %>% dplyr::select(V1, V2, dist) %>%  pivot_wider(., names_from = V1, values_from = dist) 
gd = rbind(c('Nsa1',rep(0, ncol(gd))), gd)
gd = cbind(gd, rep(0, nrow(gd)))
gd = gd[,-1] #remove row ids
gd = as.dist(gd, upper = FALSE) #take lower triangle to dist matrix object


plot(gd, m_ibsmat)
?distGeo
m_xy$dist[is.na(m_xy$total_seglen)] <- 0
xy$dist[is.na(xy$dist)] <- 0
