#isolation by distance with relatedness and fst

pkg = c('tidyverse', 'geosphere', 'data.table', 'vegan' ,'cowplot', 'sjPlot', 'Rcpp', 'RColorBrewer', 'UpSetR', 'sparseAHC', 'adegenet')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#remotes::install_github("khabbazian/sparseAHC")
#set num cores
numcores = 10

#set wd and read our shit
setwd('~/Projects/td_je_angam_2022/')
fst_by_site = read.csv('data/fst_bysite/fst_betweensites.csv')
metadata = read.csv('metadata/sequenced_metadata.csv')
relatedness = fread('data/relatedness/allsamples.res')


#fst first

#metadata subset to latlong
lalo = metadata[c('Site', 'long', 'Lat')]

#double join here (joins metadata to both sides of fst data (site a and site b))
fst_loc = merge(
  x = merge(x=fst_by_site, 
            y=lalo, 
            by.x="site_a", 
            by.y="Site")[],
  y=lalo, 
  by.x="site_b", 
  by.y = 'Site'
)[]

fst_loc = unique(fst_loc)

#make point dist between coords
fst_loc$pointdist = distVincentyEllipsoid(fst_loc[,c('Lat.x','long.x')], fst_loc[,c('Lat.y','long.y')])

#there's probably quite a bit of noise here (lots of small sample sizes etc), so let's filter our data down a bit

fst_loc$Form = fst_loc$form_a

#create a vector of sites/form that have N >= 8 (arbitrarily chosen)
filtsites = data.table(metadata)[, .N, by=.(Form, Site)][ N >= 5 ]
#create id list for filtering
site_filt_vec = paste0(filtsites$Site, '_', filtsites$Form)
fst_loc$ida = paste0(fst_loc$site_a, '_', fst_loc$Form)
fst_loc$idb = paste0(fst_loc$site_b, '_', fst_loc$Form)

#filter based on above
x = subset(fst_loc, idb %in% site_filt_vec)
y = subset(x, ida %in% site_filt_vec)
y = y[complete.cases(y), ]

formspec = y#[y$Form == c('M'),]

formspec$newfst = formspec$weighted_fst / (1-formspec$weighted_fst)
formspec$newdist = log10(formspec$pointdist)
plotly::plot_ly(data=formspec, x=~pointdist, y=~newfst, color=~Form)

filter(formspec) %>% group_by(Form) %>% summarise(min = min(newfst), max=max(newfst), mean = mean(newfst))

formspec %>% ggplot(aes(x=pointdist, y=newfst, colour=Form))+
  geom_point()+
  theme_classic()+
  scale_color_brewer(palette = 'Dark2')+
  facet_wrap(~Form)+
  labs(y='Fst / (1-Fst)', x='Geographic Distance (m)')+
#  geom_smooth(method = 'lm')+
  theme(legend.position = "none")

m = formspec[formspec$Form == "M", ]
s = formspec[formspec$Form == "S", ]

ibd_mod_m = glm(newdist ~ newfst, family = gaussian, data = s)
ibd_mod_s = glm(newdist ~ newfst, family = gaussian, data = m)

summary(ibd_mod_m)
summary(ibd_mod_s)

#insignificant, now try mantel test (also sites could be pseudoreplicates?)
#for coluzzi
m_fstmat =  xtabs(m$weighted_fst ~ m$site_b + m$site_a)
m_dstmat =  xtabs(m$pointdist ~ m$site_b + m$site_a)

#for gambiae
s_fstmat = as.dist(xtabs(s$weighted_fst ~ s$site_b + s$site_a))
s_dstmat = as.dist(xtabs(s$pointdist ~ s$site_b + s$site_a))

m_mantel_fst = vegan::mantel(m_fstmat,m_dstmat, )
s_mantel_fst = vegan::mantel(s_fstmat,s_dstmat)

m_mantel_fst
s_mantel_fst


unique(m$site_a)

#now for relatedness

allsamples_reslatedness = fread('~/Projects/td_je_angam_2022/data/relatedness/allsamples.res')
allsamples_reslatedness = left_join(allsamples_reslatedness, metadata, by=c('ida' = 'seq_id')) %>% left_join(., metadata, by=c('idb' = 'seq_id'))
allsamples_reslatedness$pointdist = distVincentyEllipsoid(allsamples_reslatedness[,c('Lat.x','long.x')], allsamples_reslatedness[,c('Lat.y','long.y')])


resM = filter(allsamples_reslatedness, Form.x == 'M' & Form.y == 'M')
resS = filter(allsamples_reslatedness, Form.x == 'S' & Form.y == 'S')


ggplot(filter(resM, KING > 0 & pointdist > 0), aes(x=pointdist, y=KING))+
  geom_point()+
  theme_classic()+
  labs(x='Geographic Distance (m)')



m_rabmat = resM %>% dplyr::select(KING, a, b) %>% pivot_wider(names_from = a, values_from=KING)
m_rabmat = m_rabmat[,-1]
m_rabmat[m_rabmat<0] <- NA
m_dismat = resM %>% dplyr::select(pointdist, a, b) %>% pivot_wider(names_from = a, values_from=pointdist)
m_dismat = m_dismat[,-1]
m_dismat[m_dismat==0] <- NA
m_mantel_res = vegan::mantel(m_rabmat,m_dismat, na.rm = TRUE)
summary(m_mantel_res)


filter(resS, pointdist > 0) %>% summarise(min = min(KING), max=max(KING), mean = mean(KING))


