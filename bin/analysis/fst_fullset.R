####analysis of
#-between form and ecotype fst
#-between form fst - downsampled
#-thetas per form and per ecotype
#thetas downsampled

#install.packages('tidygenomics')
#ls pkg
pkg = c("tidyverse", "cowplot", "RColorBrewer", "ggrepel", "viridis" , "data.table", "plotly", "RColorBrewer")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

###let's start with the by ecotype thetas
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_ecotype/')

genome_length = 280e6


##############################################
#between ecotype and isec fststs for all forms and levels
#almost same as above but data need to be cleaned a bit differently
############
# 
setwd('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/genome_scans/insecticide_usage/by_sec_and_ecotype//')


#select files, read, as usual
windowed_fst_eco_files = list.files(pattern = '*winfst')
names = gsub('window1step0p2.txt', "",  windowed_fst_eco_files)
winfst = lapply(seq_along(windowed_fst_eco_files), function(i){fread(windowed_fst_eco_files[[i]]) %>% mutate(fn=windowed_fst_eco_files[[i]])})
winfst = do.call(rbind, winfst)
colnames(winfst) = c('angsd_rubbish', 'chr', 'pos', 'nsites', 'fst','piece')
winfst$piece = gsub("_Agam.*","",winfst$piece)
#winfst$piece = gsub("decid_forest","decid-forest",winfst$piece)
#winfst = separate(winfst,col=piece, into=c('form','ecoa', 'ecob'), sep='_', remove = FALSE)

#get outliers (the data are quite noisy with lots of little local peaks so I am just going to select the most extreme windows)
my_threshold = quantile(winfst$fst, 0.95, na.rm = T)
winfst = winfst %>% mutate(outlier = ifelse(fst > my_threshold, "outlier", "background"))
#winfst = mutate(winfst, id = paste0(chr, '_', pos,'_',form, '_',piece)) 
#winfst$form = gsub("_.*","",winfst$piece)
#chrlens = aggregate(winfst$pos, by = list(winfst$chr), max) #get chromosome lengths (mx bp)
#colnames(chrlens) = c('chr', 'chrlen') #rename
#chrlens$cumsum = cumsum(chrlens$chrlen) - chrlens$chrlen #calculate cumulative positions of each chr start
#winfst = winfst[chrlens, on = .(chr)] #join to thetas data
#winfst$BPcum = winfst$pos + winfst$cumsum #get cumulative position in genome
#s = winfst[chr == 'AgamP4_X' & piece == 'm_savannah_rainforest' & outlier == 'outlier' & pos == 15240100]
#
#winfst = goi_df[winfst, on = .(id)]





####################################################################
# #plot genomewide fsts between ecotype, labelling with highly differentiated gene annotations
# #granular plots by chr and condition
# ####################################################################
# 
#prepare bed format 
s = winfst %>% filter(outlier == 'outlier') #get only outliers > 99.9th %ile
s = mutate(s, start = pos -500, end = pos + 500) #make start stop coords for bed file (from window midpoint)
s = select(s, chr, start, end, piece, fst, form, id)#select bed cols
s = s[order(chr, start, end)] #order our windows

##read annotations and intersect
genes = fread('grep protein_coding_gene /Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/annotation/VectorBase-54_AgambiaePEST.gff') #read gff file
genebed = select(genes, V1, V4,V5, V9) #convert to bed (see UCSC documentation)
genebed = genebed[order(V1, V4),] #sort just in case
colnames(genebed) = c('chr', 'start', 'end', 'name') #add colnames
s_genes = tidygenomics::genome_intersect(s, genebed) #intersect

s_genes = data.table(s_genes)
max_windows = s_genes[s_genes[, .I[which.max(fst)], by=list(piece, name)]$V1]
max_windows = max_windows[- grep("unspecified", max_windows$name),] 
winfst = max_windows[winfst, on = 'id']

winfst$id = gsub(";.*","",winfst$name)
winfst$description = gsub(".*=","",winfst$name)
winfst$description = gsub("\\[.*","",winfst$description)
#winfst$text <- paste("Gene:", winfst$description,"\nChromosome: ", winfst$i.chr, "\nWindow: ", winfst$pos, "\nFst:", winfst$i.fst %>% round(2), sep="")

#write final fst df to a compressed csv
#write.csv(x, file=gzfile('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/windowed_stats_1kb/fsts_windowed_byecotype.csv.gz'))
#winfst= fread('fsts_windowed_byecotype.csv.gz')
winfst = data.table(winfst)[,c('angsd_rubbish', 'text', 'ecoa', 'ecob', 'i.chr', 'i.piece', 'i.fst', 'i.form', 'i.name', 'i.start'):=NULL]
#write.csv(winfst, file=gzfile('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/windowed_stats_1kb/fsts_windowed_byecotype.csv.gz'))

#winfst[, c("foo","bar"):=NULL]


plot_labelled_manhattan = function(spec, condition, chrom){
  sub = x[complete.cases(x$name), ]
  sub = sub[i.chr == chrom & i.piece == condition]
  ptitle = paste0(chrom,'-', condition)
  pname = paste0(ptitle,'.tiff')
  xplot = x %>% filter(i.form == spec & i.chr == chrom & i.piece == condition) %>% 
        ggplot(aes(x=BPcum, y=i.fst, colour=outlier))+
        geom_point(size=0.3)+
        scale_color_brewer(palette='Set2')+
        theme_minimal()+
        #geom_label_repel(data=sub, aes(label=description), max.overlaps = Inf) +
      # labs(x='Fst',y='Position',colour='Chromosome')+
        labs(x='Position', y='Fst', colour = '')+
        ggtitle(pname)
    #ggplotly(xplot, tooltip="text")
  ggsave(pname, xplot, width = 15, height = 10, units = 'in', bg='white', device = grDevices::tiff, path = '~/Projects/MOVE/td_je_angam_2022/data/genome_scans/fst/')
}

#get combinations of conditions to apply over
combs = unique(as.data.frame(max_windows)[c('form', 'piece', 'chr')])

#plot
t = mapply(plot_labelled_manhattan, combs$form, combs$piece, combs$chr)

t = x %>% filter(outlier == 'outlier')
####################################################################
#plot all of our genomewide fsts
####################################################################



goi_df$goi_name
goi_df$goi_cumstart

s = x %>% filter(i.form=='m' & i.chr == 'AgamP4_X' & i.piece == 'm_savannah_rainforest' & i.fst > 0.2)

m_fst_all = winfst %>% filter(form=='s') %>% 
  ggplot(aes(x=BPcum, y=fst, colour= chr))+
  scale_color_brewer(palette='Set2')+
  geom_point(size=0.3)+
  theme_classic()+
  ylim(-0.1,1)+
  facet_wrap(~piece, ncol = 1)+
  labs(x='Cumulative Position (bp)', y = 'Fst',colour =  'Chromosome')+
  theme(legend.position="none")
ggsave('s_fst.tiff', m_fst_all, width = 12, height = 10, units = 'in',bg='white',device = grDevices::tiff)
m_fst_all


m_fst = winfst %>% filter(form == 'm') %>% ggplot(aes(x=BPcum, y=fst, colour=chr))+geom_point(size=0.3)+scale_color_brewer(palette='Set2')+theme_minimal()+labs(x='Fst',y='Position',colour='Chromosome')+facet_wrap(~comp, nrow = 6, ncol = 1)
s_fst = winfst %>% filter(form == 's') %>% ggplot(aes(x=BPcum, y=fst, colour=chr))+geom_point(size=0.3)+scale_color_brewer(palette='Set2')+theme_minimal()+labs(x='Fst',y='Position',colour='Chromosome')+facet_wrap(~comp, nrow = 6, ncol = 1)
ggsave('m_fst.tiff', m_fst, width = 12, height = 10, units = 'in',device = grDevices::tiff)
ggsave('s_fst.tiff', s_fst, width = 12, height = 10, units = 'in',device = grDevices::tiff)

########################################
#global fst between ecotypes
#####plot heatmaps of pairwise global fsts by chr for ecotype:ecotype comparison
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_ecotype/') #setwd
global_fst_fullset_files = list.files(pattern = '*globalstat.txt') 
global_fst_fullset_files = global_fst_fullset_files[2:length(global_fst_fullset_files)] 
glob_fst = lapply(seq_along(global_fst_fullset_files), function(i){fread(global_fst_fullset_files[[i]]) %>% mutate(fn=global_fst_fullset_files[[i]])}) #read files
glob_fst= do.call(rbind,glob_fst) #merge into big df
glob_fst$fn = gsub('decid_forest', "decidforest",  glob_fst$fn) #some cleaning
glob_fst = separate(glob_fst,col=fn, into=c('form','ecoone', 'ecotwo' ,'gub', 'chrgub', 'chr'), sep='_') #parse chr out of filename
glob_fst$chr = gsub('globalstat.txt', "",  glob_fst$chrgub) #clean
glob_fst = glob_fst[glob_fst$ecoone != 'form']

#set row and col orders vectors
m_column_order = c('mangrove', 'savannah', 'rainforest')
m_row_order = c('decidforest', 'rainforest', 'savannah')
glob_fst = as.data.frame(glob_fst) #coerce from dt to df
cowplot::plot_grid(ggplot(aes(x=form, y=V2, color=form), data=filter(glob_fst, chr == '3L'))+
  geom_point(size=4)+
  theme_minimal()+
  labs(x='Species', y='Weighted Fst')+
  scale_color_brewer(palette = 'Dark2')+
  theme_classic()+
  theme(legend.position = 'none',
        text = element_text(size = 15))+
  theme(legend.position = 'none')
  , labels = 'D')


#make heatmap funxn
make_fst_heatmap <- function(chrom, spec){
  f = glob_fst[(glob_fst$form == spec & glob_fst$chr == chrom),]
  f = as.data.frame(f)
  f = f[c('ecoone', 'ecotwo', 'V2')]
  f$ecoone = factor(f$ecoone, levels = m_column_order)
  f$ecotwo = factor(f$ecotwo, levels = m_row_order)
  return(f)
  #ggplot(f, aes(x=ecoone, y=ecotwo, fill=V2))+geom_tile()+labs(x='',y='')+scale_fill_viridis(discrete=FALSE)+theme_minimal()+theme(legend.position = 'NA')+geom_text(aes(label=V2))
  }

#make plots
list_m = list(make_fst_heatmap('3L', 's')) #,make_fst_heatmap('3R', 's'),make_fst_heatmap('2L', 's'),make_fst_heatmap('2R', 's'), make_fst_heatmap('X', 's'))
list_s = list(make_fst_heatmap('3L', 'm'))#,make_fst_heatmap('3R', 'm'),make_fst_heatmap('2L', 'm'),make_fst_heatmap('2R', 'm'), make_fst_heatmap('X', 'm'))
do.call(rbind, list_m)


#plot plots
m_ibe = cowplot::plot_grid(plotlist = list_m, labels = c('3L', '3R', '2L', '2R', 'X'))
s_ibe = cowplot::plot_grid(plotlist = list_s, labels = c('3L', '3R', '2L', '2R', 'X'))

s_ibe

unique(max_windows$piece)
max_windows$id = gsub(";.*","",max_windows$name)
max_windows$id = gsub("ID=","",max_windows$id)
write_csv(max_windows, file = '/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/genome_scans/fst_genomewide/genes.txt')

gen = max_windows[piece == 'ma_savannah_decid-forest']$id
as.data.frame(gen)

m_savannah_rainforest = 'insig'
m_mangrove_savannah = 'insig'

setwd('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/genome_scans/insecticide_usage/by_ecotype_and_high_isecusage_mangrove')

plot_pbs <- function(pbsdffile) {
  pbs_df = fread(pbsdffile)
  pbs_df = pbs_df[grep("Agam", chr)]
  pbs_df = pbs_df[!grep("UNK", chr)]
  pbs_df = pbs_df[!grep("unplaced", chr)]
  chrlens = aggregate(pbs_df$midPos, by = list(pbs_df$chr), max) #get chromosome lengths (mx bp)
  colnames(chrlens) = c('chr', 'chrlen') #rename
  chrlens$cumsum = cumsum(chrlens$chrlen) - chrlens$chrlen #calculate cumulative positions of each chr start
  pbs_df = pbs_df[chrlens, on = .(chr)] #join to thetas data
  pbs_df$BPcum = pbs_df$midPos + pbs_df$cumsum #get cumulative position in genome
  ggplot(pbs_df, aes(x=BPcum, y=PBS2, color=chr))+
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank())+
    geom_line(size=0.3)+
    theme_classic()+
    labs(x='', y = '', colour = '')+
    scale_color_brewer(palette='Set2')
}

filelist = list.files(pattern = '*winpbs')
filelist
plotlist = lapply(filelist, plot_pbs)

plotlist[[1]]+ theme(legend.position="none")+theme(axis.text.x = element_blank())

prow <- plot_grid(
  plotlist[[1]] + theme(legend.position="none")+ylim(0,1)+theme(axis.text.x = element_blank()),
  plotlist[[2]]  + theme(legend.position="none")+ylim(0,1)+theme(axis.text.x = element_blank()),
  plotlist[[4]]  + theme(legend.position="none")+ylim(-0.6,1)+theme(axis.text.x = element_blank()),
  plotlist[[3]]  + theme(legend.position="none")+ylim(0,1)+theme(axis.text.x = element_blank()),
  labels = c('  DF, RF, MH', '  DF, CS, MH', '  RF, CS, MH','  DF, CS, RF'),
ncol = 1
)

prow


# Write a function to draw the chromosomes on an existing plotting device
add.chromosomes <- function(gaps, cs, ce, gene.cex = 0.9, gene.col = 'grey20', point.cex = 1.2, point.col = 'grey30', chrom.col = NULL, chrom.cex = 1.4, chrom.offset = 0){
  plot(c(cs[1], ce[5]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
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
  chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
  chrom.y <- -8.5 - chrom.offset

  lines(c(cs['2L'], cs['2L'] + gaps/2, ce['2L'], ce['2L'], cs['2L'] + gaps/2, cs['2L']), 
        c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['2L'])
  text((cs['2L'] + ce['2L'])/2, chrom.y, '2L', adj = 0.5, xpd = NA, cex = chrom.cex)
  lines(c(ce['2R'], ce['2R'] - gaps/2, cs['2R'], cs['2R'], ce['2R'] - gaps/2, ce['2R']), 
        c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['2R'])
  text((cs['2R'] + ce['2R'])/2, chrom.y, '2R', adj = 0.5, xpd = NA, cex = chrom.cex)
  lines(c(cs['3L'], cs['3L'] + gaps/2, ce['3L'], ce['3L'], cs['3L'] + gaps/2, cs['3L']), 
        c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['3L'])
  text((cs['3L'] + ce['3L'])/2, chrom.y, '3L', adj = 0.5, xpd = NA, cex = chrom.cex)
  lines(c(ce['3R'], ce['3R'] - gaps/2, cs['3R'], cs['3R'], ce['3R'] - gaps/2, ce['3R']), 
        c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['3R'])
  text((cs['3R'] + ce['3R'])/2, chrom.y, '3R', adj = 0.5, xpd = NA, cex = chrom.cex)
  lines(c(cs['X'], cs['X'], ce['X'] - gaps/2, ce['X'], ce['X'], ce['X'] - gaps/2, cs['X']), 
        c(-1, 1, 1, 0.2, -0.2, -1, -1), lwd = 2, col = chrom.col['X'])
  text((cs['X'] + ce['X'])/2, chrom.y, 'X', adj = 0.5, xpd = NA, cex = chrom.cex)
}
add.chromosomes(gaps = gaps, cs = cs, ce = ce, chrom.offset = -1, chrom.cex = 1.3)




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
