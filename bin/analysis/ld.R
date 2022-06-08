#plotting LD blocks
#ls pkg
pkg = c("data.table", "tidyverse", 'gtools','LDheatmap', 'reshape2')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

#read ld data
ld = read.table('~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/AgamP4_2R_glf2_minorfilts_angsd_samtools_ngsLD')
r = ld %>% select(V1, V3, V5, V6, V7, V8, V9) 
colnames(r) = c('snp1', 'snp2', 'dist', 'r2p', 'D', 'Dp', 'r2')


#id <- unique(mixedsort(c(r[,"snp1"],r[,"snp2"])))
id = unique(cbind(r$snp1, r$snp2))


mixedsort(id)


posStart <- head(id,1)
posEnd <- tail(id,1)
r <- rbind(r, c(posStart,posStart,0,NA,NA,NA,NA), c(posEnd,posEnd,0,NA,NA,NA,NA))
r = r %>% select(snp1, snp2, r2)
m = acast(r, snp1 ~ snp2, drop=FALSE)




r$r2 = as.numeric(r$r2)
r = r %>% select(snp1, snp2, r2) %>% pivot_wider(names_from = snp1, values_from = r2)
m = as.matrix(r)
row.names(m) <- m[,1]
m = m[,-1]
t<- m[mixedorder(rownames(m)),mixedorder(colnames(m))]
id = rownames(m)
dist <- as.numeric(sub(".*:","",id))

LDheatmap(m, genetic.distances=dist, SNP.name=id, geneMapLabelX=0.75, geneMapLabelY=0.25, color="blueToRed", LDmeasure=ld)





for (ld in c("r2","Dp")) {
  m <- apply(acast(r, snp1 ~ snp2, value.var=ld, drop=FALSE),2,as.numeric)
  rownames(m) <- colnames(m)
  m <- m[mixedorder(rownames(m)),mixedorder(colnames(m))]
  id <- rownames(m)
  dist <- as.numeric(sub(".*:","",id))
  # Save plot
  pdf(paste("LD_blocks", ld,"pdf", sep="."), width=10, height=10)
  LDheatmap(m, genetic.distances=dist, SNP.name=id, geneMapLabelX=0.75, geneMapLabelY=0.25, color="blueToRed", LDmeasure=ld)
  x <- dev.off()
}
