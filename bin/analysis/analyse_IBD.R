#####
#Analysis of ROH tract data
#tristan dennis 23.10.21
#####

#load our shit
#parse_pcangsd.py output
#run this in the same directory 
pkg = c('tidyverse', 'RcppCNPy', 'data.table', 'qqman' ,'cowplot', 'parallel', 'ragg', 'RColorBrewer', 'UpSetR')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

#set num cores
numcores = 10
#setwd
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_full_cov')

