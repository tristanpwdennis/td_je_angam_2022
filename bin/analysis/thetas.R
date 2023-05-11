library(tidyverse)
library(data.table)


thetalist = list.files(path = '/Users/dennisTajimaw/Projects/td_je_angam_2022/data/thetas/', pattern = '.thetasWindow.gz.pesTajimaG', full.names = T)
thetas = lapply(thetalist, fread)
thetas$Pi <- thetas


par(mfrow = c(7, 5), mar = c(0.1, 0, 1.2, 0), oma=c(4,3,2,2))
#m mangrove decid forest
plot(x=thetas[[1]]$WinCenter, y=thetas[[1]]$tP/thetas[[1]]$nSites,type='l', ylim=c(0,0.5), xaxt="n",col='#7fc97f',frame.plot=FALSE, las=2, main = "Anc DF", adj = 0.2) #2L
plot(x=thetas[[2]]$WinCenter, y=thetas[[2]]$tP/thetas[[2]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#7fc97f') #2R
plot(x=thetas[[3]]$WinCenter, y=thetas[[3]]$tP/thetas[[3]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#7fc97f') #3L
plot(x=thetas[[4]]$WinCenter, y=thetas[[4]]$tP/thetas[[4]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#7fc97f') #3R
plot(x=thetas[[5]]$WinCenter, y=thetas[[5]]$tP/thetas[[5]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#7fc97f') #X
#m mangrove rainforest
plot(x=thetas[[6]]$WinCenter, y=thetas[[6]]$tP/thetas[[6]]$nSites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#beaed4',frame.plot=FALSE,las = 2, main = "Anc M", adj = 0.2) #2L
plot(x=thetas[[7]]$WinCenter, y=thetas[[7]]$tP/thetas[[7]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#beaed4') #2R) #2R
plot(x=thetas[[8]]$WinCenter, y=thetas[[8]]$tP/thetas[[8]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#beaed4') #3L) #3L
plot(x=thetas[[9]]$WinCenter, y=thetas[[9]]$tP/thetas[[9]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#beaed4') #3R) #3R
plot(x=thetas[[10]]$WinCenter, y=thetas[[10]]$tP/thetas[[10]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#beaed4') #XSE) #X
#m mangrove savannah
plot(x=thetas[[11]]$WinCenter, y=thetas[[11]]$tP/thetas[[11]]$nSites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#fdc086',frame.plot=FALSE,las = 2, main = "Anc RF", adj = 0.2) #2L
plot(x=thetas[[12]]$WinCenter, y=thetas[[12]]$tP/thetas[[12]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#fdc086') #2R) #2R
plot(x=thetas[[13]]$WinCenter, y=thetas[[13]]$tP/thetas[[13]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fdc086') #3L) #3L
plot(x=thetas[[14]]$WinCenter, y=thetas[[14]]$tP/thetas[[14]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fdc086') #3R) #3R
plot(x=thetas[[15]]$WinCenter, y=thetas[[15]]$tP/thetas[[15]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fdc086') #XSE) #X
#m_rainforest_decid_forest
plot(x=thetas[[16]]$WinCenter, y=thetas[[16]]$tP/thetas[[16]]$nSites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#fb8072',frame.plot=FALSE,las = 2, main = "Anc CS", adj = 0.2) #2L
plot(x=thetas[[17]]$WinCenter, y=thetas[[17]]$tP/thetas[[17]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#fb8072') #2R) #2R
plot(x=thetas[[18]]$WinCenter, y=thetas[[18]]$tP/thetas[[18]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fb8072') #3L) #3L
plot(x=thetas[[19]]$WinCenter, y=thetas[[19]]$tP/thetas[[19]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fb8072') #3R) #3R
plot(x=thetas[[20]]$WinCenter, y=thetas[[20]]$tP/thetas[[20]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#fb8072') #XSE) #X
#m_savannah_decid_forWinCenter
plot(x=thetas[[21]]$WinCenter, y=thetas[[21]]$tP/thetas[[21]]$nSites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#ffff99',frame.plot=FALSE, las = 2, main = "Ang DF", adj = 0.2) #2L
plot(x=thetas[[22]]$WinCenter, y=thetas[[22]]$tP/thetas[[22]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col = '#ffff99') #2R) #2R
plot(x=thetas[[23]]$WinCenter, y=thetas[[23]]$tP/thetas[[23]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#ffff99') #3L) #3L
plot(x=thetas[[24]]$WinCenter, y=thetas[[24]]$tP/thetas[[24]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#ffff99') #3R) #3R
plot(x=thetas[[25]]$WinCenter, y=thetas[[25]]$tP/thetas[[25]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=FALSE, col='#ffff99') #XSE) #X
#m_savannah_rainforesWinCenter
plot(x=thetas[[26]]$WinCenter, y=thetas[[26]]$tP/thetas[[26]]$nSites,type='l', ylim=c(0,0.5), xaxt="n", ann=TRUE, col='#f0027f',frame.plot=FALSE,las = 2, main = "Ang RF", adj = 0.2) #2L
plot(x=thetas[[27]]$WinCenter, y=thetas[[27]]$tP/thetas[[27]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=TRUE, col = '#f0027f',frame.plot=FALSE,las = 2) #2R) #2R
plot(x=thetas[[28]]$WinCenter, y=thetas[[28]]$tP/thetas[[28]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=TRUE, col='#f0027f',frame.plot=FALSE, las=2) #3L) #3L
plot(x=thetas[[29]]$WinCenter, y=thetas[[29]]$tP/thetas[[29]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=TRUE, col='#f0027f',frame.plot=FALSE, las=2) #3R) #3R
plot(x=thetas[[30]]$WinCenter, y=thetas[[30]]$tP/thetas[[30]]$nSites,type='l', ylim=c(0,0.5), axes = F, ann=TRUE, col='#f0027f',frame.plot=FALSE, las=2) #XSE) #X
#m_savannah_rainforest
plot(x=thetas[[31]]$WinCenter, y=thetas[[31]]$tP/thetas[[31]]$nSites,type='l', ylim=c(0,0.5), axes= T, ann=TRUE, col='#bf5b17',frame.plot=FALSE,las = 2, main = "Ang CS", adj = 0.2) #2L
plot(x=thetas[[32]]$WinCenter, y=thetas[[32]]$tP/thetas[[32]]$nSites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col = '#bf5b17',frame.plot=FALSE,las = 2) #2R) #2R
plot(x=thetas[[33]]$WinCenter, y=thetas[[33]]$tP/thetas[[33]]$nSites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#bf5b17',frame.plot=FALSE, las=2) #3L) #3L
plot(x=thetas[[34]]$WinCenter, y=thetas[[34]]$tP/thetas[[34]]$nSites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#bf5b17',frame.plot=FALSE, las=2) #3R) #3R
plot(x=thetas[[35]]$WinCenter, y=thetas[[35]]$tP/thetas[[35]]$nSites,type='l', ylim=c(0,0.5), yaxt="n", ann=TRUE, col='#bf5b17',frame.plot=FALSE, las=2) #XSE) #X


par(mfrow = c(7, 5), mar = c(0.1, 0, 1.2, 0), oma=c(4,3,2,2))
#m mangrove decid forest
plot(x=thetas[[1]]$WinCenter, y=thetas[[1]]$Tajima,type='l', ylim=c(-2.5,1), xaxt="n",col='#7fc97f',frame.plot=FALSE, las=2, main = "Anc DF", adj = 0.2) #2L
plot(x=thetas[[2]]$WinCenter, y=thetas[[2]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col = '#7fc97f') #2R
plot(x=thetas[[3]]$WinCenter, y=thetas[[3]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#7fc97f') #3L
plot(x=thetas[[4]]$WinCenter, y=thetas[[4]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#7fc97f') #3R
plot(x=thetas[[5]]$WinCenter, y=thetas[[5]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#7fc97f') #X
#m mangrove rainforest
plot(x=thetas[[6]]$WinCenter, y=thetas[[6]]$Tajima,type='l', ylim=c(-2.5,1), xaxt="n", ann=TRUE, col='#beaed4',frame.plot=FALSE,las = 2, main = "Anc M", adj = 0.2) #2L
plot(x=thetas[[7]]$WinCenter, y=thetas[[7]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col = '#beaed4') #2R) #2R
plot(x=thetas[[8]]$WinCenter, y=thetas[[8]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#beaed4') #3L) #3L
plot(x=thetas[[9]]$WinCenter, y=thetas[[9]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#beaed4') #3R) #3R
plot(x=thetas[[10]]$WinCenter, y=thetas[[10]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#beaed4') #XSE) #X
#m mangrove savannah
plot(x=thetas[[11]]$WinCenter, y=thetas[[11]]$Tajima,type='l', ylim=c(-2.5,1), xaxt="n", ann=TRUE, col='#fdc086',frame.plot=FALSE,las = 2, main = "Anc RF", adj = 0.2) #2L
plot(x=thetas[[12]]$WinCenter, y=thetas[[12]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col = '#fdc086') #2R) #2R
plot(x=thetas[[13]]$WinCenter, y=thetas[[13]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#fdc086') #3L) #3L
plot(x=thetas[[14]]$WinCenter, y=thetas[[14]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#fdc086') #3R) #3R
plot(x=thetas[[15]]$WinCenter, y=thetas[[15]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#fdc086') #XSE) #X
#m_rainforest_decid_forest
plot(x=thetas[[16]]$WinCenter, y=thetas[[16]]$Tajima,type='l', ylim=c(-2.5,1), xaxt="n", ann=TRUE, col='#fb8072',frame.plot=FALSE,las = 2, main = "Anc CS", adj = 0.2) #2L
plot(x=thetas[[17]]$WinCenter, y=thetas[[17]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col = '#fb8072') #2R) #2R
plot(x=thetas[[18]]$WinCenter, y=thetas[[18]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#fb8072') #3L) #3L
plot(x=thetas[[19]]$WinCenter, y=thetas[[19]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#fb8072') #3R) #3R
plot(x=thetas[[20]]$WinCenter, y=thetas[[20]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#fb8072') #XSE) #X
#m_savannah_decid_forWinCenter
plot(x=thetas[[21]]$WinCenter, y=thetas[[21]]$Tajima,type='l', ylim=c(-2.5,1), xaxt="n", ann=TRUE, col='#ffff99',frame.plot=FALSE, las = 2, main = "Ang DF", adj = 0.2) #2L
plot(x=thetas[[22]]$WinCenter, y=thetas[[22]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col = '#ffff99') #2R) #2R
plot(x=thetas[[23]]$WinCenter, y=thetas[[23]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#ffff99') #3L) #3L
plot(x=thetas[[24]]$WinCenter, y=thetas[[24]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#ffff99') #3R) #3R
plot(x=thetas[[25]]$WinCenter, y=thetas[[25]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=FALSE, col='#ffff99') #XSE) #X
#m_savannah_rainforesWinCenter
plot(x=thetas[[26]]$WinCenter, y=thetas[[26]]$Tajima,type='l', ylim=c(-2.5,1), xaxt="n", ann=TRUE, col='#f0027f',frame.plot=FALSE,las = 2, main = "Ang RF", adj = 0.2) #2L
plot(x=thetas[[27]]$WinCenter, y=thetas[[27]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=TRUE, col = '#f0027f',frame.plot=FALSE,las = 2) #2R) #2R
plot(x=thetas[[28]]$WinCenter, y=thetas[[28]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=TRUE, col='#f0027f',frame.plot=FALSE, las=2) #3L) #3L
plot(x=thetas[[29]]$WinCenter, y=thetas[[29]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=TRUE, col='#f0027f',frame.plot=FALSE, las=2) #3R) #3R
plot(x=thetas[[30]]$WinCenter, y=thetas[[30]]$Tajima,type='l', ylim=c(-2.5,1), axes = F, ann=TRUE, col='#f0027f',frame.plot=FALSE, las=2) #XSE) #X
#m_savannah_rainforest
plot(x=thetas[[31]]$WinCenter, y=thetas[[31]]$Tajima,type='l', ylim=c(-2.5,1), axes= T, ann=TRUE, col='#bf5b17',frame.plot=FALSE,las = 2, main = "Ang CS", adj = 0.2) #2L
plot(x=thetas[[32]]$WinCenter, y=thetas[[32]]$Tajima,type='l', ylim=c(-2.5,1), yaxt="n", ann=TRUE, col = '#bf5b17',frame.plot=FALSE,las = 2) #2R) #2R
plot(x=thetas[[33]]$WinCenter, y=thetas[[33]]$Tajima,type='l', ylim=c(-2.5,1), yaxt="n", ann=TRUE, col='#bf5b17',frame.plot=FALSE, las=2) #3L) #3L
plot(x=thetas[[34]]$WinCenter, y=thetas[[34]]$Tajima,type='l', ylim=c(-2.5,1), yaxt="n", ann=TRUE, col='#bf5b17',frame.plot=FALSE, las=2) #3R) #3R
plot(x=thetas[[35]]$WinCenter, y=thetas[[35]]$Tajima,type='l', ylim=c(-2.5,1), yaxt="n", ann=TRUE, col='#bf5b17',frame.plot=FALSE, las=2) #XSE) #X
