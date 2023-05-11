library(tidyverse)
library(data.table)

#wrangle
flist = list.files(path = '/Users/dennistpw/Projects/td_je_angam_2022/data/fst_by_ecotype_species/', pattern = 'window5step1.txt', full.names = T)
flistnam<-  str_remove(basename(flist), "_AgamP4.*")
flistnam<-  str_remove(basename(flistnam), "_forest")
names(flist) <- flistnam
fstlist <- map_dfr(flist, fread, .id = 'ID')
fstlist$species = gsub('_.*','',fstlist$ID)
fstlist$ecoregionb <- gsub('.*_','',fstlist$ID)
fstlist$ecoregiona <- gsub('_.*','',gsub('[sm]_','',fstlist$ID))

unique(fstlist$ID)

fstlist<-fstlist %>%
  mutate(comp = case_when( 
    ID == 'm_mangrove_decid' ~ "M:DF", 
    ID == 'm_mangrove_rainforest' ~ "M:RF", 
    ID == 'm_mangrove_savannah' ~ "M:CS", 
    ID == 'm_rainforest_decid' ~ "RF:DF", 
    ID == 'm_savannah_decid' ~ "CS:DF", 
    ID == 'm_savannah_rainforest' ~ "CS:RF", 
    ID == 's_rainforest_decid' ~ "RF:DF", 
    ID == 's_savannah_decid' ~ "CS:DF", 
    ID == 's_savannah_rainforest' ~ "CS:RF", 
    TRUE ~ NA
  ))



#fstlist = fstlist[,-2:-1]

fstlist %>% case_when(
  ecoregiona == 'decid' & ecoregionb == 'decid' ~ "M:DF")




write.csv(fstlist, '~/Projects/td_je_angam_2022/data/fst_by_ecotype_species/fst_window5step1_allecoregions_species.csv')

data = '~/Projects/td_je_angam_2022/data/fst_by_ecotype_species/fst_window5step1_allecoregions_species.txt'
fst = readr::read_delim(data)
ggplot(fst, aes(x=chr,))





par(mfrow = c(6, 5), mar = c(0.1, 0, 1.2, 0), oma=c(4,3,2,2))
#m mangrove decid forest
plot(x=fstlist[[1]]$chr, y=fstlist[[1]]$Nsites,type='l', ylim=c(0,0.5), xaxt="n",col='#8dd3c7',frame.plot=FALSE, las=2, main = "M:DF", adj = 0.2) #2L
plot(x=fstlist[[2]]$chr, y=fstlist[[2]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#8dd3c7') #2R
plot(x=fstlist[[3]]$chr, y=fstlist[[3]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#8dd3c7') #3L
plot(x=fstlist[[4]]$chr, y=fstlist[[4]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#8dd3c7') #3R
plot(x=fstlist[[5]]$chr, y=fstlist[[5]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#8dd3c7') #X
#m mangrove rainforest
plot(x=fstlist[[6]]$chr, y=fstlist[[6]]$Nsites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#feeea6',frame.plot=FALSE,las = 2, main = "M:RF", adj = 0.2) #2L
plot(x=fstlist[[7]]$chr, y=fstlist[[7]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#feeea6') #2R) #2R
plot(x=fstlist[[8]]$chr, y=fstlist[[8]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#feeea6') #3L) #3L
plot(x=fstlist[[9]]$chr, y=fstlist[[9]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#feeea6') #3R) #3R
plot(x=fstlist[[10]]$chr, y=fstlist[[10]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#feeea6') #XSE) #X
#m mangrove savannah
plot(x=fstlist[[11]]$chr, y=fstlist[[11]]$Nsites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#80b1d3',frame.plot=FALSE,las = 2, main = "M:CS", adj = 0.2) #2L
plot(x=fstlist[[12]]$chr, y=fstlist[[12]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#80b1d3') #2R) #2R
plot(x=fstlist[[13]]$chr, y=fstlist[[13]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#80b1d3') #3L) #3L
plot(x=fstlist[[14]]$chr, y=fstlist[[14]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#80b1d3') #3R) #3R
plot(x=fstlist[[15]]$chr, y=fstlist[[15]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#80b1d3') #XSE) #X
#m_rainforest_decid_forest
plot(x=fstlist[[11]]$chr, y=fstlist[[11]]$Nsites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#fb8072',frame.plot=FALSE,las = 2, main = "RF:DF", adj = 0.2) #2L
plot(x=fstlist[[17]]$chr, y=fstlist[[17]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#fb8072') #2R) #2R
plot(x=fstlist[[18]]$chr, y=fstlist[[18]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fb8072') #3L) #3L
plot(x=fstlist[[19]]$chr, y=fstlist[[19]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fb8072') #3R) #3R
plot(x=fstlist[[20]]$chr, y=fstlist[[20]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fb8072') #XSE) #X
#m_savannah_decid_forest
plot(x=fstlist[[11]]$chr, y=fstlist[[11]]$Nsites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#bebaea',frame.plot=FALSE, las = 2, main = "CS:DF", adj = 0.2) #2L
plot(x=fstlist[[22]]$chr, y=fstlist[[22]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#bebaea') #2R) #2R
plot(x=fstlist[[23]]$chr, y=fstlist[[23]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#bebaea') #3L) #3L
plot(x=fstlist[[24]]$chr, y=fstlist[[24]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#bebaea') #3R) #3R
plot(x=fstlist[[25]]$chr, y=fstlist[[25]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#bebaea') #XSE) #X
#m_savannah_rainforest
plot(x=fstlist[[26]]$chr, y=fstlist[[26]]$Nsites,type='l', ylim=c(0,0.5), axes= T, ann=TRUE, col='#fdb462',frame.plot=FALSE,las = 2, main = "CS:RF", adj = 0.2) #2L
plot(x=fstlist[[27]]$chr, y=fstlist[[27]]$Nsites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col = '#fdb462',frame.plot=FALSE,las = 2) #2R) #2R
plot(x=fstlist[[28]]$chr, y=fstlist[[28]]$Nsites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#fdb462',frame.plot=FALSE, las=2) #3L) #3L
plot(x=fstlist[[29]]$chr, y=fstlist[[29]]$Nsites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#fdb462',frame.plot=FALSE, las=2) #3R) #3R
plot(x=fstlist[[30]]$chr, y=fstlist[[30]]$Nsites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#fdb462',frame.plot=FALSE, las=2) #XSE) #X

par(mfrow = c(6, 5), mar = c(0.1, 0, 1.2, 0), oma=c(4,3,2,2))
#s_rainforest_decid_forest
plot(x=fstlist[[31]]$chr, y=fstlist[[31]]$Nsites,type='l', ylim=c(0,0.5), xaxt="n",col='#b3de69',frame.plot=FALSE, las=2, main = "RF:DF", adj = 0.2) #2L
plot(x=fstlist[[32]]$chr, y=fstlist[[32]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#b3de69') #2R
plot(x=fstlist[[33]]$chr, y=fstlist[[33]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#b3de69') #3L
plot(x=fstlist[[34]]$chr, y=fstlist[[34]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#b3de69') #3R
plot(x=fstlist[[35]]$chr, y=fstlist[[35]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#b3de69') #X
#s_savannah_decid_forest
plot(x=fstlist[[36]]$chr, y=fstlist[[36]]$Nsites,type='l', ylim=c(0,0.5), xaxt="n",col='#fccdf9',frame.plot=FALSE, las=2, main = "CS:DF", adj = 0.2) #2L
plot(x=fstlist[[37]]$chr, y=fstlist[[37]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#fccdf9') #2R) #2R
plot(x=fstlist[[38]]$chr, y=fstlist[[38]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fccdf9') #3L) #3L
plot(x=fstlist[[39]]$chr, y=fstlist[[39]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fccdf9') #3R) #3R
plot(x=fstlist[[40]]$chr, y=fstlist[[40]]$Nsites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fccdf9') #XSE) #X
#s_savannah_rainforest
plot(x=fstlist[[41]]$chr, y=fstlist[[41]]$Nsites,type='l', ylim=c(0,0.5), axes= T,col='#d9d9e9',frame.plot=FALSE, las=2, main = "CS:RF", adj = 0.2) #2L
plot(x=fstlist[[42]]$chr, y=fstlist[[42]]$Nsites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#d9d9e9',frame.plot=FALSE, las=2) #XSE) #X
plot(x=fstlist[[43]]$chr, y=fstlist[[43]]$Nsites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#d9d9e9',frame.plot=FALSE, las=2) #XSE) #X
plot(x=fstlist[[44]]$chr, y=fstlist[[44]]$Nsites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#d9d9e9',frame.plot=FALSE, las=2) #XSE) #X
plot(x=fstlist[[45]]$chr, y=fstlist[[45]]$Nsites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#d9d9e9',frame.plot=FALSE, las=2) #XSE) #X






