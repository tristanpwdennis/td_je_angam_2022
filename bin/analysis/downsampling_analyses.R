library(tidyverse)
library(data.table)
library(ggthemes)

metadata=fread('~/Projects/td_je_angam_2022/metadata/metadata_with_insecticide_info.csv')
corebamlist<- fread('/Users/dennistpw/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/td_je_angam/td_je_angam_2022/data/downsampling/pca/corebamlist.csv', header=T)
colnames(corebamlist) <- c('bamno1','seq_id')

corebamlist$bamno <- seq(0, nrow(corebamlist) -1)
metadata <- left_join(corebamlist, metadata, by=c('seq_id' ='seq_id' ))

fstlist = list.files(path='~/Projects/td_je_angam_2022/data/downsampling/', pattern = '*window*', full.names = T)
fstf <- lapply(fstlist, function(file_path) {
  # Use fread to read in the file
  file_data <- fread(file_path)
  
  # Add a column with the filename
  file_data[, filename := basename(file_path)]
  
  # Return the modified data.table
  return(file_data)
})


par(mfrow = c(6, 1), mar = c(0.1, 0, 0.8, 0), oma=c(4,3,3,2), adj = 0.1)
#10x
#plot(x=fstf[[11]]$chr, y=fstf[[11]]$Nsites,type='l', ylim=c(0,1), xaxt="n",col='#0868ac',frame.plot=FALSE, las=2, main = "10X") #2L
#plot(x=fstf[[12]]$chr, y=fstf[[12]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col = '#0868ac') #2R
#plot(x=fstf[[13]]$chr, y=fstf[[13]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#0868ac') #3L
#plot(x=fstf[[14]]$chr, y=fstf[[14]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#0868ac') #3R
plot(x=fstf[[15]]$chr, y=fstf[[15]]$Nsites,type='l', ylim=c(0,1), xaxt="n", col='#0868ac',frame.plot=FALSE, las=2, main = "10X") #X
#7.5
#plot(x=fstf[[26]]$chr, y=fstf[[26]]$Nsites,type='l', ylim=c(0,1), xaxt="n",col='#2b8cbe',frame.plot=FALSE, las=2, main = "7.5X") #2L
#plot(x=fstf[[27]]$chr, y=fstf[[27]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col = '#2b8cbe') #2R
#plot(x=fstf[[28]]$chr, y=fstf[[28]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#2b8cbe') #3L
#plot(x=fstf[[29]]$chr, y=fstf[[29]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#2b8cbe') #3R
#plot(x=fstf[[30]]$chr, y=fstf[[30]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#2b8cbe') #X
plot(x=fstf[[30]]$chr, y=fstf[[30]]$Nsites,type='l', ylim=c(0,1), xaxt="n", col='#2b8cbe',frame.plot=FALSE, las=2, main = "7.5X") #X

#5
#plot(x=fstf[[21]]$chr, y=fstf[[21]]$Nsites,type='l', ylim=c(0,1), xaxt="n",col='#4eb3d3',frame.plot=FALSE, las=2, main = "5X") #2L
#plot(x=fstf[[22]]$chr, y=fstf[[22]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col = '#4eb3d3') #2R
#plot(x=fstf[[23]]$chr, y=fstf[[23]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#4eb3d3') #3L
#plot(x=fstf[[24]]$chr, y=fstf[[24]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#4eb3d3') #3R
#plot(x=fstf[[25]]$chr, y=fstf[[25]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#4eb3d3') #X
plot(x=fstf[[25]]$chr, y=fstf[[25]]$Nsites,type='l', ylim=c(0,1), xaxt="n", col='#4eb3d3',frame.plot=FALSE, las=2, main = "5X") #X

#2.5
#plot(x=fstf[[16]]$chr, y=fstf[[16]]$Nsites,type='l', ylim=c(0,1), xaxt="n",col='#7bccc4',frame.plot=FALSE, las=2, main = "2.5X") #2L
#plot(x=fstf[[17]]$chr, y=fstf[[17]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col = '#7bccc4') #2R
#plot(x=fstf[[18]]$chr, y=fstf[[18]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#7bccc4') #3L
#plot(x=fstf[[19]]$chr, y=fstf[[19]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#7bccc4') #3R
#plot(x=fstf[[20]]$chr, y=fstf[[10]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#7bccc4') #X
plot(x=fstf[[20]]$chr, y=fstf[[20]]$Nsites,type='l', ylim=c(0,1), xaxt="n", col='#7bccc4',frame.plot=FALSE, las=2, main = "2.5X") #X
#1
#plot(x=fstf[[6]]$chr, y=fstf[[6]]$Nsites,type='l', ylim=c(0,1), xaxt="n",col='#a8ddb5',frame.plot=FALSE, las=2, main = "1X") #2L
#plot(x=fstf[[7]]$chr, y=fstf[[7]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col = '#a8ddb5') #2R
#plot(x=fstf[[8]]$chr, y=fstf[[8]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#a8ddb5') #3L
#plot(x=fstf[[9]]$chr, y=fstf[[9]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#a8ddb5') #3R
#plot(x=fstf[[9]]$chr, y=fstf[[9]]$Nsites,type='l', ylim=c(0,1), axes = F, ann=FALSE, col='#a8ddb5') #X
plot(x=fstf[[10]]$chr, y=fstf[[10]]$Nsites,type='l', ylim=c(0,1), xaxt="n", col='#a8ddb5',frame.plot=FALSE, las=2, main = "1X") #X
#0.5
#plot(x=fstf[[1]]$chr, y=fstf[[1]]$Nsites,type='l', ylim=c(0,0.5), axes= T, ann=TRUE, col='#ccebc5',frame.plot=FALSE,las = 2, main = "0.5X") #2L
#plot(x=fstf[[2]]$chr, y=fstf[[2]]$Nsites,type='l', ylim=c(0,1), yaxt="n", ann=FALSE, col = '#ccebc5',frame.plot=FALSE, las=2) #2R
#plot(x=fstf[[3]]$chr, y=fstf[[3]]$Nsites,type='l', ylim=c(0,1), yaxt="n", ann=FALSE, col='#ccebc5',frame.plot=FALSE, las=2) #3L
#plot(x=fstf[[4]]$chr, y=fstf[[4]]$Nsites,type='l', ylim=c(0,1), yaxt="n", ann=FALSE, col='#ccebc5',frame.plot=FALSE, las=2) #3R
#plot(x=fstf[[5]]$chr, y=fstf[[5]]$Nsites,type='l', ylim=c(0,1), yaxt="n", ann=FALSE, col='#ccebc5',frame.plot=FALSE, las=2) #X
plot(x=fstf[[5]]$chr, y=fstf[[5]]$Nsites,type='l', ylim=c(0,1), xaxt="n", col='#ccebc5',frame.plot=FALSE, las=2, main = "0.5X") #X


xfst <-do.call(rbind, list(
fstf[[15]],
fstf[[30]],
fstf[[25]],
fstf[[20]],
fstf[[10]],
fstf[[5]]
))

xfst$filename <- gsub('m_form_s_form_','',xfst$filename)
xfst$filename <- gsub('_AgamP4_Xwindow5step1.txt','',xfst$filename)
xfst$filename <- factor(xfst$filename, levels = c("10", "7.5", "5", "2.5", "1", "0.5"))

fst <- ggplot(xfst, aes(x=chr,y=Nsites,colour=filename))+
  geom_line()+
  facet_wrap(~filename, ncol = 1)+
  scale_color_manual(values=pal)+
  theme_classic()+
  labs(x="Position", y='Fst')+
  theme(legend.position = 'none')
  


##relkatednes
reslist = list.files(path='~/Projects/td_je_angam_2022/data/downsampling/', pattern = '*.res', full.names = T)
resf <- lapply(reslist, function(file_path) {
  # Use fread to read in the file
  file_data <- fread(file_path)
  
  # Add a column with the filename
  file_data[, filename := basename(file_path)]
  
  # Return the modified data.table
  return(file_data)
})
pal =   c("#0868ac", "#2b8cbe", "#4eb3d3", "#7bccc4", "#a8ddb5", "#ccebc5")

resf <- do.call(rbind, resf)
length(unique(resf$a))
resf$doc <- as.numeric(gsub("bamlist.*","",resf$filename))
resf$doc <- factor(resf$doc, levels = c("10", "7.5", "5", "2.5", "1", "0.5"))

resf  = left_join(resf, metadata, by=c('a' = 'bamno')) %>% left_join(., metadata, by=c('b' = 'bamno'))
#filtresf = resf[resf$Site.x == resf$Site.y]
filtresf = resf[resf$Form.x != resf$Form.y]

colres <- resf[resf$Form.x == 'M' & resf$Form.y == 'M']

rab <- ggplot(colres, aes(x=doc, y=rab, colour=as.factor(doc)))+
  geom_jitter(alpha=0.7)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  theme_classic()+
  labs(x="Mean Depth-of-Coverage",y="Rxy")+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")

KING <- ggplot(colres, aes(x=doc, y=KING, colour=as.factor(doc)))+
  geom_jitter(alpha=0.7)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  theme_classic()+
  labs(x="Mean Depth-of-Coverage",y="KING")+
  scale_color_manual(values=pal)+
  theme(legend.position = "none")


#now pca


pcadown<-fread('/Users/dennistpw/Library/Mobile Documents/com~apple~CloudDocs/Projects/td_je_angam/td_je_angam_2022/data/downsampling/pca/pca_downsampling.csv', header=F)
pcadown <- pcadown[-1,]
pcadown$V220 <- as.numeric(pcadown$V220)
pcadown = left_join(pcadown, metadata, by=c('V220'='bamno'))

pcadown$V221 <- gsub('AgamP4_3L_glf2_downsampled_','',pcadown$V221)
pcadown$coverage <- gsub('0','0.5',pcadown$V221)
pcadown$coverage <- str_replace_all(pcadown$coverage, c("7" = "7.5", "2" = "2.5", "10.5"="10"), )
pcadown$species = str_replace_all(pcadown$Form, c("M"="An coluzzi", "S" = 'An gambiae'))

pcadown$coverage <- factor(pcadown$coverage, levels=c('10','7.5','5','2.5','1','0.5'))
pcad = ggplot(pcadown[!is.na(pcadown$species)], aes(x=as.numeric(V2),y=as.numeric(V3), colour=coverage))+
  facet_wrap(~coverage)+
  geom_point()+
  scale_color_manual(values=pal)+
  theme_classic()+
  labs(x='PC1',y='PC2')+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

prowa <- cowplot::plot_grid(rab, pcad, labels = c('A','B'),ncol = 2)
prowb <- cowplot::plot_grid(prowa, fst, labels = c('','C'), nrow = 2, rel_heights = c(0.6, 1))

prowb


##pca
covlist = list.files(path='~/Projects/td_je_angam_2022/data/downsampling/', pattern = '*.cov*', full.names = T)
pcaf <- lapply(covlist, function(file_path) {
  # Use fread to read in the file
  file_data <- fread(file_path)
  
  # Add a column with the filename
  file_data[, filename := basename(file_path)]
  
  # Return the modified data.table
  return(file_data)
})
plot(pcaf[[1]]$V1, pcaf[[1]]$V2)
covlist
