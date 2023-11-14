#isolation by distance with relatedness and fst
pkg = c('tidyverse', 'geosphere', 'data.table', 'vegan' ,'cowplot', 'RColorBrewer')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#remotes::install_github("khabbazian/sparseAHC")
#set num cores
numcores = 10

#set wd and read our shit
setwd('~/Dropbox (LSTM)/td_je_angam_2022/')
fst_by_site = read.csv('~/Dropbox (LSTM)/td_je_angam_2022/data/fst_bysite/fst_betweensites.csv')
metadata = read.csv('metadata/sequenced_metadata.csv')

#read coluzzi data
m_res = fread('~/Dropbox (LSTM)/td_je_angam_2022/data/relatedness/m_form.res')
m_samples <- metadata[metadata$Form == 'M',]
m_samples$seq_id <- seq(0, nrow(m_samples)-1)

#read gambiae data
s_res = fread('~/Dropbox (LSTM)/td_je_angam_2022/data/relatedness/s_form.res')
s_samples <- metadata[metadata$Form == 'S',]
s_samples$seq_id <- seq(0, nrow(s_samples)-1)

####coluzzi
mres_meta = left_join(m_res, m_samples, by=c('a' = 'seq_id')) %>% left_join(., m_samples, by=c('b' = 'seq_id'))
mres_meta$pointdist = distVincentyEllipsoid(mres_meta[,c('Lat.x','long.x')], mres_meta[,c('Lat.y','long.y')])
mking <-mres_meta %>% mutate(sameordiff = ifelse(habitat.x ==habitat.y, 'Same','Diff')) %>% 
  ggplot(aes(y=KING, fill=sameordiff))+
  geom_histogram(bins = 100, colour='darkblue')+
  scale_fill_manual(values=c("dodgerblue4","dodgerblue1"))+
  theme_classic()+
  xlim(0, 2000)+
  ylim(-2,0.5)+
  labs(x='Count')+
  theme(legend.position = "none",
        axis.text.y=element_blank(), 
        axis.title.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.ticks.y=element_blank())

####gambiae
sres_meta = left_join(s_res, s_samples, by=c('a' = 'seq_id')) %>% left_join(., s_samples, by=c('b' = 'seq_id'))
sres_meta$pointdist = distVincentyEllipsoid(sres_meta[,c('Lat.x','long.x')], sres_meta[,c('Lat.y','long.y')])
sking<- sres_meta %>% mutate(sameordiff = ifelse(habitat.x ==habitat.y, 'Same','Diff')) %>% 
  ggplot(aes(y=KING, fill=sameordiff))+
  geom_histogram(bins = 100, colour='#013220')+
  scale_fill_manual(values=c('#037d50','#05c880'))+
  theme_classic()+
  labs(x='Count')+
  xlim(0, 2000)+
  ylim(-1.2,0.5)+
  theme(legend.position = "none",
        axis.text.y=element_blank(), 
        axis.line.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y=element_blank())

#gambiae geo/gendist
sisobd<- sres_meta %>% mutate(sameordiff = ifelse(habitat.x ==habitat.y, 'Same','Diff')) %>% 
  ggplot(., aes(x=pointdist, y=KING, colour=sameordiff))+
  geom_point(alpha=0.8, size=1.5)+
  scale_colour_manual(values=c('#037d50','#05c880'))+
  theme_classic()+
  ylim(-1.2,0.5)+
  xlim(0, 526267)+
  labs(x='Geographic Distance (m)')+
  theme(legend.position = "none")

#coluzzi geo/gendist
misobd<- mres_meta %>% mutate(sameordiff = ifelse(habitat.x ==habitat.y, 'Same','Diff')) %>% 
  ggplot(., aes(x=pointdist, y=rab, colour=sameordiff))+
  geom_point(alpha=0.8, size=1.5)+
  scale_colour_manual(values=c("dodgerblue4","dodgerblue1"))+
  theme_classic()+
  ylim(-2,0.5)+
  xlim(0, 526267)+
  labs(x='Geographic Distance (m)')+
  theme(legend.position = "none")

#plot final figure
sibfig <- cowplot::plot_grid(misobd, mking, sisobd,sking, labels=c('A','B','C','D'))

#mantel tests
m_rabmat = mres_meta %>% dplyr::select(KING, a, b) %>% pivot_wider(names_from = a, values_from=KING)
m_rabmat = m_rabmat[,-1] #remove idb column from matrix
#m_rabmat[m_rabmat<0] <- NA
m_dismat = mres_meta %>% dplyr::select(pointdist, a, b) %>% pivot_wider(names_from = a, values_from=pointdist)
m_dismat = m_dismat[,-1]
#m_dismat[m_dismat==0] <- NA
m_mantel_res = vegan::mantel(m_rabmat,m_dismat, na.rm = TRUE)
m_mantel_res

s_rabmat = sres_meta %>% dplyr::select(KING, a, b) %>% pivot_wider(names_from = a, values_from=KING)
s_rabmat = s_rabmat[,-1] #remove idb column from matrix
#s_rabmat[s_rabmat<0] <- NA
s_dismat = sres_meta %>% dplyr::select(pointdist, a, b) %>% pivot_wider(names_from = a, values_from=pointdist)
s_dismat = s_dismat[,-1]
#s_dismat[s_dismat==0] <- NA
s_mantel_res = vegan::mantel(s_rabmat,s_dismat, na.rm = TRUE)
s_mantel_res
