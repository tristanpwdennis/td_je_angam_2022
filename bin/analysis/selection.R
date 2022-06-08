#parse_pcangsd.py output
#run this in the same directory 
pkg = c('tidyverse', 'RcppCNPy', 'data.table', 'qqman')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

#
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled')


selection = list.files(path = ".", pattern = ".selection.npy")
admix = list.files(path = ".", pattern = ".admix.")

#load selection
s = RcppCNPy::npyLoad('AgamP4_2L_glf2_minorfilts_angsd_samtools_downsampled_1.selection.npy')
p =  fread("AgamP4_2L_glf2_minorfilts_angsd_samtools_downsampled_1.pos.gz")[ ,1:2]


## function for QQplot
qqchi<-function(x,...){
  lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
  legend("topleft",paste("lambda=",lambda))
}

#make 11plot to qc the test stat
qqchi(s)

# convert test statistic to p-value
pval<-1-pchisq(s,1)

names(p)<-c("chr","pos")

colnames(pval) <-c('pc1')

selchr  = cbind(p, data.frame(pval))

ggplot(selchr, aes(x=pos,y=-log10(pc1)))+
  geom_point()+
  theme_minimal()







