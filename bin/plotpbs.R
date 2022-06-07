####analysis of
#-between form and ecotype fst
#-between form fst - downsampled
#-thetas per form and per ecotype
#thetas downsampled

#install.packages('tidygenomics')
#ls pkg
pkg = c("tidyverse", "cowplot", "RColorBrewer", "ggrepel", "viridis" , "data.table", "ragg", "parallel")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

###let's start with the by ecotype thetas
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_ecotype/')

colpal = c('AgamP4_2L', 'AgamP4_2R', 'AgamP4_3L', 'AgamP4_3R', 'AgamP4_X')

colpal = RColorBrewer::brewer.pal(5, "Set2")

setwd('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/genome_scans/insecticide_usage/by_ecotype_decidlow_s/')


pbsdffile = '/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/genome_scans/insecticide_usage/by_ecotype_decidlow_s/rainforest_savannah_decidlow_s.winpbs'


##-----------------------------------------------------------------------------------------------------------##
##                                           Anopheles coluzzi PBS                                           ##
##-----------------------------------------------------------------------------------------------------------##

#function to plot pbs
make_pbs_df <- function(pbsdffile) {
  pbs_df = fread(pbsdffile)
  colnames(pbs_df) <- c('Angsd_gubbins', 'chr', 'pos', 'nsites', 'fst01', 'fst02', 'fst12', 'PBS0', 'PBS1', 'PBS2')
  pbs_df = pbs_df[grep("Agam", chr)]
  pbs_df = pbs_df[!grep("UNK", chr)]
  pbs_df = pbs_df[!grep("unplaced", chr)]
  chrlens = aggregate(pbs_df$pos, by = list(pbs_df$chr), max) #get chromosome lengths (mx bp)
  colnames(chrlens) = c('chr', 'chrlen') #rename
  chrlens$cumsum = cumsum(chrlens$chrlen) - chrlens$chrlen #calculate cumulative positions of each chr start
  pbs_df = pbs_df[chrlens, on = .(chr)] #join to thetas data
  pbs_df$BPcum = pbs_df$pos + pbs_df$cumsum #get cumulative position in genome
  #pbs_df <- sample_n(pbs_df, size = nrow(pbs_df) * 0.5)
  return(pbs_df)
}


#ancol (figure)
mform_mh <- list.files('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/genome_scans/insecticide_usage/by_ecotype_and_high_isecusage_mangrove/', full.names = T)
hilist <- mclapply(mform_mh, make_pbs_df, mc.cores = 6)

#make chr colour palette
pal = c('#d73027', '#fc8d59', '#91bfdb', '#4575b4', '#fee090')

#df: mangrove high
ac_df_mh = ggplot(hilist[[1]], aes(x=BPcum, y=fst02, colour = chr))+
  geom_line()+
  ylim(0, 1)+
  theme_classic()+
  scale_color_manual(values = pal)+
  theme(legend.position = "none", axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        plot.margin = margin(0,0,0,1, "cm"))+  labs(x='Position', y='Fst')
#rf: mangrove high
ac_rf_mh = ggplot(hilist[[1]], aes(x=BPcum, y=fst12, colour = chr))+
  geom_line()+
  ylim(0, 1)+
  theme_classic()+
  scale_color_manual(values = pal)+
  theme(legend.position = "none", axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        plot.margin = margin(0,0,0,1, "cm"))+
  labs(x='Position', y='Fst')

#rf: mangrove high
ac_cs_mh = ggplot(hilist[[2]], aes(x=BPcum, y=fst12, colour = chr))+
  geom_line()+
  ylim(0, 1)+
  theme_classic()+
  scale_color_manual(values = pal)+
  theme(legend.position = "none", axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        plot.margin = margin(0,0,0,1, "cm"))+  labs(x='Position', y='Fst')

ac_fig <- plot_grid(ac_df_mh, ac_rf_mh, ac_cs_mh, ncol =1, labels = "AUTO")
#ac_fig = egg::ggarrange(ac_df_mh, ac_rf_mh, ac_cs_mh, ncol = 1, nrow = 3, labels=c('A', 'B', 'C'))
#ggsave(plot=ac_fig,filename = '~/Projects/MOVE/td_je_angam_2022/figures/fig2_ancol_fst.tiff', device = grDevices::tiff, width = 20, height = 10, units = 'cm')
tiff("~/Projects/MOVE/td_je_angam_2022/figures/fig2_ancol_fst.tiff")
print(ac_fig)
dev.off()
rm(hilist)

#angam (figure)
mform_mh <- list.files('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/genome_scans/insecticide_usage/by_ecotype_and_low_isecusage_mangrove/', full.names = T)
hilist <- mclapply(mform_mh[[7]], make_pbs_df, mc.cores = 6)

#df: rf and dfl
ag_rf_dfl = ggplot(hilist[[1]], aes(x=BPcum, y=fst02, colour = chr))+
  geom_line()+
  ylim(0, 1)+
  theme_classic()+
  scale_color_manual(values = pal)+
  theme(legend.position = "none", axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        plot.margin = margin(0,0,0,1, "cm"))+  labs(x='Position', y='Fst')
#rf: cs and dfl
ag_cs_dfl = ggplot(hilist[[1]], aes(x=BPcum, y=fst12, colour = chr))+
  geom_line()+
  ylim(0, 1)+
  theme_classic()+
  scale_color_manual(values = pal)+
  theme(legend.position = "none", axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        plot.margin = margin(0,0,0,1, "cm"))+
  labs(x='Position', y='Fst')

ag_fig <- plot_grid(ag_rf_dfl, ag_cs_dfl, ncol =1, labels = "AUTO")

#cs: decid low
#rf: decid low


mform_mh <- list.files('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/genome_scans/insecticide_usage/by_sec_and_ecotype/', full.names = T)


#plot chromosomes


chrom.sizes <- c('2L' = 49364325, '2R' = 61545105,  '3L' = 41963435, '3R' = 53200684, 'X' = 24393108)
gaps=2000000
ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
# get a vector of start points for each chromosome (ie: cumulative sizes + gaps)
cs <- ce - chrom.sizes
chrom.offset = -1
chrom.cex = 1.3
gene.cex = 0.9
gene.col = 'grey20'
point.cex = 1.2
point.col = 'grey30'
chrom.col = NULL
chrom.cex = 1.4
chrom.offset = 0
plot(c(cs[1], ce[5]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
cs
cs['2L']
# show the kdr region
kdr.region.mean <- cs['2L'] + mean(c(2358158, 2431617))
points(kdr.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
text(kdr.region.mean, -1.7, 'Vgsc', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
# show the Rdl region
rdl.region.mean <- cs['2L'] + mean(c(25363652, 25434556))
points(rdl.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
text(rdl.region.mean, -1.7, 'Rdl', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
# show the Carboxylesterase region
coeae.region.mean <- cs['2L'] + mean(c(28548433, 28550748))
points(coeae.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
text(coeae.region.mean, -1.7, 'Coeae2f', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
# show the Ace1 gene region
ace1.region.mean <- cs['2R'] + mean(c(3484107, 3495790))
points(ace1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
text(ace1.region.mean, -1.7, 'Ace1', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
# show the CYP6 region
cyp6.region.mean <- cs['2R'] + mean(c(28463000, 28568000))
points(cyp6.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
text(cyp6.region.mean, -1.7, 'Cyp6aa1-\nCyp6p2 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
# show the GST region
gst.region.mean <- cs['3R'] + mean(c(28580000, 28605000))
points(gst.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
text(gst.region.mean, -1.7, 'Gstu4-\nGste3 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
# show the CYP6M2-Z1 region
cyp6m2.z1.region.mean <- cs['3R'] + mean(c(6900000, 7030000))
points(cyp6m2.z1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
text(cyp6m2.z1.region.mean, -1.7, 'Cyp6m2-\nCyp6z1 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
# show the CYP9K1 region
cyp9k1.region.mean <- cs['X'] + mean(c(15222000, 15257000))
points(cyp9k1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
text(cyp9k1.region.mean, -1.7, 'Cyp9k1', srt = 45, adj = 1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
# Plot the outline of the chromosomes
chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2L' = 'black', '2R' = 'black', '3L' = 'black', '3R' = 'black', 'X' = 'black')
chrom.y <- -8.5 - chrom.offset

lines(c(ce['2L'], ce['2L'] - gaps/2, cs['2L'], cs['2L'], ce['2L'] - gaps/2, ce['2L']), 
      c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['2L']) #c(-0.2, -1, -1, 1, 1, 0.2)
text((cs['2L'] + ce['2L'])/2, chrom.y, '2L', adj = 0.5, xpd = NA, cex = chrom.cex)

lines(c(cs['2R'], cs['2R'] + gaps/2, ce['2R'], ce['2R'], cs['2R'] + gaps/2, cs['2R']), 
      c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['2R'])
text((cs['2R'] + ce['2R'])/2, chrom.y, '2R', adj = 0.5, xpd = NA, cex = chrom.cex)

lines(c(ce['3L'], ce['3L'] - gaps/2, cs['3L'], cs['3L'], ce['3L'] - gaps/2, ce['3L']), 
      c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['3L'])
text((cs['3L'] + ce['3L'])/2, chrom.y, '3L', adj = 0.5, xpd = NA, cex = chrom.cex)


lines(c(cs['3R'], cs['3R'] + gaps/2, ce['3R'], ce['3R'], cs['3R'] + gaps/2, cs['3R']), 
      c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['3R'])
text((cs['3R'] + ce['3R'])/2, chrom.y, '3R', adj = 0.5, xpd = NA, cex = chrom.cex)



lines(c(cs['X'], cs['X'], ce['X'] - gaps/2, ce['X'], ce['X'], ce['X'] - gaps/2, cs['X']), 
      c(-1, 1, 1, 0.2, -0.2, -1, -1), lwd = 2, col = chrom.col['X'])
text((cs['X'] + ce['X'])/2, chrom.y, 'X', adj = 0.5, xpd = NA, cex = chrom.cex)

plot(1,2)






