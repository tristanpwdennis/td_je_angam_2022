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

#create poplists for angsd analysis, output as txt files for upload to hpc
#for (form in c('M', 'S')) {
#  for (site in unique(metadata$Site)) {
#    x = metadata[metadata$Site == site,]
#    t = x[x$Form == form,]$bamid
#    write(t, paste0('~/Projects/MOVE/data/anopheles_03_21_DW/metadata/', site, '_', form, '.pop'))
#  }
#}


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
filtsites = data.table(metadata)[, .N, by=.(Form, Site)][ N >= 10 ]
#create id list for filtering
site_filt_vec = paste0(filtsites$Site, '_', filtsites$Form)
fst_loc$ida = paste0(fst_loc$site_a, '_', fst_loc$Form)
fst_loc$idb = paste0(fst_loc$site_b, '_', fst_loc$Form)

#filter based on above
x = subset(fst_loc, idb %in% site_filt_vec)
y = subset(x, ida %in% site_filt_vec)
y = y[complete.cases(y), ]


formspec = y#[y$Form == c('M'),]



#formspec = y[y$Form == c('M'),]
formspec$newfst = formspec$weighted_fst / (1-formspec$weighted_fst)
formspec$newdist = log10(formspec$pointdist)

formspec %>% ggplot(aes(x=pointdist, y=newfst, colour=Form))+
  geom_point()+
  theme_classic()+
  scale_color_brewer(palette = 'Dark2')+
  facet_wrap(~Form)+
  labs(y='Fst / (1-Fst)', x='Geographic Distance')+
  geom_smooth(method = 'lm')+
  theme(legend.position = "none")

m = formspec[formspec$Form == "M", ]
s = formspec[formspec$Form == "S", ]

ibd_mod = glm(newdist ~ newfst, family = gaussian, data = m)
summary(ibd_mod)


e = as.dist(xtabs(formspec$weighted_fst ~ formspec$site_b + formspec$site_a))
f = as.dist(xtabs(formspec$newdist ~ formspec$site_b + formspec$site_a))




plot(e, f)
abline(lm(e~f), col="red",lty=2)

t = vegan::mantel(e,f)
t
m_test = mantel.randtest(e,f)
m_test
plot(m_test)

t = glm(newfst ~ newdist, data = formspec, family='gaussian')
summary(t)

library(MASS)
dens <- kde2d(e,f, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(e, f, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(e~f))
title("Isolation by distance plot")


ggplot(formspec, aes(y=newdist,x=newfst))+
  geom_point()+
  facet_wrap(~Form)
plot(e,f)
##create fst dist mat
##do ibd test
#fstdistmat = as.matrix(as.dist(xtabs(formspec$weighted_fst ~ formspec$site_a + formspec$site_b)))
#pdistmat = as.matrix(as.dist(xtabs(formspec$pointdist ~ formspec$site_a + formspec$site_b)))
#
#length(unique(formspec$site_b))
#length(unique(formspec$site_a))
#
#pd = formspec %>% select(site_a, site_b, weighted_fst) %>% pivot_wider(names_from = site_a, values_from = weighted_fst)
#pd
#as.matrix(pd)
#
#matrix(formspec$weighted_fst, nrow = 20, ncol = 20)
#
##ibd = mantel.randtest(as.dist(fstdistmat), as.dist(pdistmat), nrepet = 999)
##ibd
#ibd = mantel(fstdistmat, pdistmat, method = "spearman", permutations = 9999, na.rm = TRUE)
#summary(ibd)
#plot(abund_temp)
#fstdistmat
#
##plot the results on a map as a graph
#
#sites = metadata %>% select(Site, habitat) %>% distinct() #vertices
#distances = fst_loc %>% select(site_a, site_b, unweighted_fst, Form) #edges
#sites


