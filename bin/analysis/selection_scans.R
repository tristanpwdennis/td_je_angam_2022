library(pcadapt)
library(tidyverse)
library(data.table)
library(parallel)






pcadapt <- as.matrix(fread('~/Projects/larvae_paper_stash/pcangsd_selection/m_mangrove.pop.AgamP4_X.beagle.gz.pcangsd.pcadapt.zscores.npy.txt'))
D <- fread('~/Projects/larvae_paper_stash/pcangsd_selection/m_mangrove.pop.AgamP4_X.beagle.gz.pcangsd.selection.npy.txt')
Dsite <- fread('~/Projects/larvae_paper_stash/pcangsd_selection/m_mangrove.pop.AgamP4_X.beagle.gz.site.gz')
Dsite <- separate(Dsite, marker, into = c('piece','chrom','pos'))
p <- pchisq(D$V1, 1, lower.tail=FALSE)
p<-cbind(Dsite,p)
p$logp <- -log10(p$p)
p$pos <- as.numeric(p$pos)

qval <- qvalue(p$p)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)

hist(-log10(p$p))


ggplot(p, aes(x=pos, y=logp))+
  geom_point()


#now let's look at pcadapt
library(bigutilsr)

K <- ncol(pcadapt)

# For one component only
if (K == 1) {
  d2 <- (pcadapt - median(pcadapt))^2
} else {
  d2 <- dist_ogk(pcadapt)
}
d2p <- pchisq(d2, df=K, lower.tail=F)

pcadaptp<-cbind(Dsite,d2p)
pcadaptp$pos <- as.numeric(pcadaptp$pos)

ggplot(pcadaptp, aes(x=pos, y=-log10(V1)))+
  geom_point()


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


get_outlier_threshold <- function(column){
  outlierserthreshold <- getmode(column) + (3*(getmode(column) - min(column)))
  return(outlierserthreshold)
}



get_outlier_threshold(pcadaptp$V1)
get_outlier_threshold(-log10(pcadaptp$V1))

hist(-log10(pcadaptp$V1), breaks=100)

pcadaptp %>% mutate(modal_value = get_outlier_threshold(-log10(V1)),
       new_column = ifelse(-log10(V1) > modal_value, "Outlier", "NotOutlier")) %>% 
  ggplot(aes(x=pos, y=-log10(V1), colour=new_column))+
    geom_point()


library(qvalue)
qval <- qvalue(pcadaptp$V1)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)

pcadaptp$padj <- p.adjust(-log10(pcadaptp$V1),method="bonferroni")
alpha <- 0.1
pcadaptp %>% mutate(outlier = ifelse(V1 < 0.01, "Outlier", "NotOutlier")) %>% 
  ggplot(aes(x=pos,y=-log10(V1), colour=outlier))+
    geom_point()

summary(-log10(pcadaptp$V1))

hist(pcadaptp$V1, breaks=50)

pcadaptp$padj < alpha
length(outliers)



pcadaptp[outliers,]

#for each ecoregion and species
#read for each chromosome the pcadapt file and sites file
#do dist normalisation
#generate pvalues
#generale histograms and outlier position lists
#bind pcadapt and sites file
#bind  all chr/sites together
#plot with outliers marked

#list of genes in outlier positions

