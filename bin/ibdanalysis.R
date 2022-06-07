#####
#Analysis of IBD tract data
#tristan dennis 23.10.21
#####

#load our shit
#parse_pcangsd.py output
#run this in the same directory 
pkg = c('tidyverse', 'geosphere', 'data.table', 'qqman' ,'cowplot', 'parallel', 'Rcpp', 'RColorBrewer', 'UpSetR', 'sparseAHC', 'moments')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#remotes::install_github("khabbazian/sparseAHC")
#set num cores
numcores = 10
#setwd
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_fullset/ibd_analysis/')


#export list of pops for fst ibd analysis
source('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/bin/generate_metadata.R')
metadata = prepare_metadata()
#sitelist = metadata %>% group_by(Form, Agric_eco.activity) %>% count() 
#
#
#for (site in unique(sitelist$Site)) {
#  for (form in unique(metadata$Form)) {
#    s = metadata[metadata$Form == form,]
#    t = s[s$Site == site,]
#    x = t$bamid
#    write(x, paste0('~/Projects/MOVE/data/anopheles_03_21_DW/metadata/', site, '_', form, '.pop'))
#  }  
#}
?write.csv

#metadata = read.delim('~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/metadata_v1.txt')
#bamlist = read.csv('~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/bamlist.txt', header =F)
#bamlist$no = seq(0, (length(bamlist$V1) - 1))
#colnames(bamlist) = c('sample_id', 'bamname', 'samplenum')
#newmeta = left_join(metadata, bamlist, by=c('sample' = 'sample_id'))

#get our imputed files (vs unimputed ones)
ibdfiles = list.files(pattern='\\.ibd.gz$')
hbdfiles = list.files(pattern='\\.hbd.gz$')
nefiles = list.files(pattern='\\.ne$')


look = fread('sform_3L_phasednoref-0.01.ibd.gz')
look$sample_id = paste0(look$V1, '_', look$V3)


#calculate total ibd length sharing per individual pair and plot boxplots and over distance
get_ibdlen <- function(ibdfile) {
  ibdseqfile = fread(ibdfile, colClasses = c('character', 'numeric', 'character', 'numeric', 'character', 'numeric', 'numeric', 'numeric'))
  ibdseqfile$len = ibdseqfile[,7] - ibdseqfile[,6]
  totalibd = do.call(data.frame, aggregate(. ~ V1+V3+V5, data = ibdseqfile, FUN = function(ibdseqfile) c(sumlen = sum(ibdseqfile), n = length(ibdseqfile) ) ))
  totalibd = totalibd %>% select(V1, V3, len.sumlen, len.n)
  totalibd$fn = ibdfile
  totalibd = separate(totalibd, fn, into = c('form', 'chr', 'treatment'), sep = '_')
}

#get everything and stick together, change names to numbers for metadata
agg_ibd = do.call(rbind, mclapply(ibdfiles, get_ibdlen, mc.cores =10))

#faff around and clean up names, split tx col into cm and tx, make function next time
agg_ibd$cm = sub('.*\\-', '', agg_ibd$treatment)
agg_ibd$treatment = sub("-[^-]+$", "", agg_ibd$treatment)
agg_ibd$cm = sub('\\.ibd.gz', '', agg_ibd$cm)
agg_ibd$ida = gsub('ind', '', agg_ibd$V1)
agg_ibd$idb = gsub('ind', '', agg_ibd$V3)
agg_ibd$ida =  as.numeric(agg_ibd$ida) 
agg_ibd$idb =  as.numeric(agg_ibd$idb) 

#agg_ibd %>% select(cm) %>% distinct()

sform = filter(metadata, Form == 'S') %>% mutate(population = 's_form') %>% mutate(bam_key = seq(0, nrow(.)-1))
mform = filter(metadata, Form == 'M') %>% mutate(population = 'm_form') %>% mutate(bam_key = seq(0, nrow(.)-1))
#join to formspecific metadata
fullset_df = rbind(agg_ibd %>% filter(form == 'mform') %>% left_join(., mform, by = c('ida' = 'bam_key')) %>% left_join(., mform, by = c('idb' = 'bam_key')),
                   agg_ibd %>% filter(form == 'sform') %>% left_join(., sform, by = c('ida' = 'bam_key')) %>% left_join(., sform, by = c('idb' = 'bam_key')) )

#defin within/between ecozone pairing
fullset_df = fullset_df %>% mutate(withinorbetween = case_when(
  habitat.x == habitat.y ~ 'Same Ecozone',
  habitat.x != habitat.y ~ 'Different Ecozone'
))
#calc dist between samples and add as col to fullsetdf
dist = distVincentyEllipsoid(fullset_df[,c('Lat.x','long.x')], fullset_df[,c('Lat.y','long.y')])
fullset_df$pointdist = dist

#add chrlen
fullset_df = fullset_df %>% mutate(chrlen = case_when(
  chr == '3L' ~ 41946761,
  chr == '3R' ~ 53200684,
))

#calculate fraction of chromsome in ibd tracts
fullset_df$chrfrac = fullset_df$len.sumlen / fullset_df$chrlen

#remove gubbins
#plottingdf = fullset_df %>% select(ida, idb, len.n, len.sumlen, habitat.x, habitat.y, Lat.x, long.x, Lat.y, long.y, Insecticide.x, Insecticide.y, form, chr, cm, treatment, chrfrac, pointdist)

#plot ibd
ibd_iso_s = fullset_df %>% filter(form == 'sform') %>% ggplot(aes(x=pointdist,y=chrfrac))+facet_grid(rows = vars(cm), cols = vars(treatment)) + geom_point(colour= '#1B9E77')+theme_minimal()+labs(y='Shared fraction of Chr IBDesc', x='Distance (m)')+ylim(0,1)
ibd_iso_m = fullset_df %>% filter(form == 'mform') %>% ggplot(aes(x=pointdist,y=chrfrac))+facet_grid(rows = vars(cm), cols = vars(treatment)) + geom_point(colour= '#1B9E77')+theme_minimal()+labs(y='Shared fraction of Chr IBDesc', x='Distance (m)')+ylim(0,1)




#plot isoibd for specific pair
fullset_df %>% filter((treatment == 'imputed-popspecificref' | treatment == 'phasednoref') & cm == 0.1) %>% ggplot(aes(x=pointdist,y=chrfrac, color=form))+facet_wrap(~form+treatment)+geom_point()+theme_minimal()+scale_color_brewer(palette = 'Dark2')+labs(y='Shared fraction of Chr IBDesc', x='Distance (m)')

fullset_df %>% filter(( treatment == 'phasednoref' & cm == 0.05))  %>% 
  ggplot(aes(x=pointdist,y=chrfrac, color=form))+
  geom_point()+
  facet_wrap(~form)+
  theme_minimal()+
  scale_color_brewer(palette = 'Dark2')+
  labs(y='Shared fraction of Chr IBDesc', x='Distance (m)', colour = 'Form')


m = fullset_df %>% filter( treatment == 'phasednoref' & cm == 0.05 & form == 'mform')
mmod = glm(data = m, formula = chrfrac ~ pointdist, family = gaussian)
summary(mmod)
hist(m$chrfrac)


s = fullset_df %>% filter(cm == 0.1) 
unique(s$treatment)


#isolation by environment
ibedf = data.frame(fullset_df %>% select(chrfrac, habitat.x, habitat.y, form, chrfrac, withinorbetween) %>% group_by(habitat.x, form, habitat.y, withinorbetween) %>% summarise(medianfrac = median(chrfrac)) )
ggplot(ibedf, aes(x=habitat.x, habitat.y, fill=medianfrac))+geom_tile()+facet_wrap(~form+withinorbetween)+theme_minimal()

#newres = res %>% filter(species.x == 'S' & species.y == 'S')
#x = fullset_df %>% filter(form == 'm') %>% dplyr::select(ida, idb, dist)  %>% pivot_wider(., names_from = ida, values_from = dist)
#y = fullset_df %>% dplyr::select(a, b, V54) %>% pivot_wider(., names_from = a, values_from = V54)

fullset_df[is.na(fullset_df)] <- 0

sdf = fullset_df %>% filter(treatment == 'phasednoref' & form == 'mform' & cm == 0.05) %>% select(ida, idb, chrfrac, pointdist)



x = xtabs(sdf$chrfrac ~ sdf$idb + sdf$ida)
y = xtabs(sdf$pointdist ~ sdf$idb + sdf$ida)
s = as.matrix(x)
t = as.matrix(y)
x
view(x)


ibd = ade4::mantel.randtest(as.dist(s), as.dist(t), nrepet = 999)
plot(ibd)
ibd

#graphy clustery bits
#add transposed matrix (y?)
s = Matrix(s, sparse = TRUE)
A = s + t(s)
#calc adjacency graph
G = igraph::graph.adjacency(A, mode = "undirected", weighted=TRUE)
plot(G, edge.label=E(G)$weight, vertex.label=V(G)-1)
H <- sparseAHC(A, "average", TRUE)
plot(H)
####plot histograms of LOD scores

#function taking ibdseq output, reading, calculating per-pair ibd lengt

m_ibd_files = mclapply(mfiles, fread, mc.cores = 10)
s_ibd_files = mclapply(sfiles, fread, mc.cores = 10)
m_ibd_files = lapply(seq_along(m_ibd_files), function(i){m_ibd_files[[i]] %>% mutate(fn = mfiles[[i]])})
s_ibd_files = lapply(seq_along(s_ibd_files), function(i){s_ibd_files[[i]] %>% mutate(fn = sfiles[[i]])})
m_df = do.call(rbind, m_ibd_files)
s_df = do.call(rbind, s_ibd_files)
fulldf  = rbind(m_df, s_df)
fulldf = separate(fulldf, fn, into = c('impute', 'mincm', 'form', 'gub', 'gub2', 'chr'), sep = '_')
fulldf %>%  ggplot(aes(x=form, y=V8, fill=form)) +facet_wrap(~mincm+impute, ncol=2) +geom_violin()+theme_minimal()+theme(legend.position="none")+labs(x='LOD Score',y='Count')+scale_fill_brewer(palette='Dark2')
fulldf$len = fulldf$V7 - fulldf$V6


s_imputed = fulldf %>% filter(form=='s') %>% ggplot(aes(y=V8)) +facet_wrap(~mincm+impute, ncol=2)+geom_histogram( fill = '#1B9E77', colour = '#1B9E77')+theme_minimal()+theme(legend.position="none")+labs(x='LOD Score',y='Count')
m_imputed = fulldf %>% filter(form=='m') %>% ggplot(aes(y=V8)) +facet_wrap(~mincm+impute, ncol=2)+geom_histogram( fill = '#D95F02', colour = '#D95F02')+theme_minimal()+theme(legend.position="none")+labs(x='LOD Score',y='Count')

plot_grid(m_imputed, s_imputed, labels = c('S', 'M'))

fulldf %>% filter(impute == 'imputed' & mincm=0.01) %>% ggplot(aes(x=form,y=V8))+geom_jitter()
fulldf %>% group_by(form, mincm, impute) %>% count() %>% ggplot(aes(x=form,y=n, fill=form))+facet_wrap(~mincm+ impute, ncol=4)+geom_bar(stat='identity')+theme_minimal()+scale_fill_brewer(palette='Dark2')+labs(x='Form',y='Count IBD Segments')+theme(legend.position="none")


##hbd output
hbdfiles = list.files(pattern='*.hbd.gz')
hbdfiles
lapply(hbdfiles, data.table::fread)
test= fread(hbdfiles[[1]])
test$len=test$V7-test$V6
test %>% group_by(V3) %>% summarise(total_len = sum(len), count_frags = n()) %>% ggplot(aes(x=total_len, y=count_frags))+geom_point()


s_hbd_files = mclapply(hbdsfiles, fread, mc.cores = 10)
m_hbd_files = lapply(seq_along(m_hbd_files), function(i){m_hbd_files[[i]] %>% mutate(fn = hbdmfiles[[i]])})
s_hbd_files = lapply(seq_along(s_hbd_files), function(i){s_hbd_files[[i]] %>% mutate(fn = hbdsfiles[[i]])})
m_df = do.call(rbind, m_hbd_files)
s_df = do.call(rbind, s_hbd_files)
hbdf  = rbind(m_df, s_df)
hbdf = separate(hbdf, fn, into = c('impute', 'mincm', 'form', 'gub', 'gub2', 'chr'), sep = '_')
hbdf$len = hbdf$V7 - hbdf$V6
aggregate(hbdf, hbdf$len ~V1+ impute, FUN=sum)

hbdf %>% select(impute, V1) %>% group_by(impute, V1) %>% count()
hbdf %>% select(impute, len) %>% group_by(impute, len) %>% across(.funs = sum)


setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_fullset/ibd_analysis/')
nefiles = list.files(pattern='*.ne')
bootfiles = list.files(pattern='*.boot')
nefiles
nefile = fread(nefiles[[1]])
nefile %>% ggplot(aes(x=GEN,y=NE))+geom_line()+theme_minimal()+xlim(0,100)
nefiles






