library(data.table)

setwd('~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_chr_downsampled/')
metadata = read.csv('~/Projects//MOVE/td_je_angam_2022/metadata/metadata_with_insecticide_info.csv')
metadata %>% dplyr::filter(INSECTICIDE_USE == 'high' & Form == 'M' & habitat == 'mangrove_swamp')

path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_chr_downsampled/'
?fread




mform_sform_fst = rbind(fread(paste0(path,'m_form_s_form_full_AgamP4_2Lwindow5step1.txt')),
fread(paste0(path,'m_form_s_form_full_AgamP4_2Rwindow5step1.txt')),
fread(paste0(path,'m_form_s_form_full_AgamP4_3Lwindow5step1.txt')),
fread(paste0(path,'m_form_s_form_full_AgamP4_3Rwindow5step1.txt')),
fread(paste0(path,'m_form_s_form_full_AgamP4_Xwindow5step1.txt'))
)
colnames(mform_sform_fst) <- c('angsd_gubbins', 'chr', 'pos', 'nsites', 'fst')
chrlens = aggregate(mform_sform_fst$pos, by = list(mform_sform_fst$chr), max) #get chromosome lengths (mx bp)
colnames(chrlens) = c('chr', 'chrlen') #rename
chrlens$cumsum = cumsum(chrlens$chrlen) - chrlens$chrlen #calculate cumulative positions of each chr start
mform_sform_fst = mform_sform_fst[chrlens, on = .(chr)] #join to thetas data
mform_sform_fst$BPcum = mform_sform_fst$pos + mform_sform_fst$cumsum #get cumulative position in genome
ggplot(mform_sform_fst, aes(x=BPcum, y=fst, colour=chr))+
  theme_classic()+
  scale_color_brewer(palette='Set2')+
  geom_point(size=0.8)+
  ylim(-0.05, 1)+
  theme(legend.position = 'none',
        text = element_text(size = 15))+
  labs(x='', y='Fst')

#zoomed in onto chrx
chrxonly = ggplot(mform_sform_fst[chr == 'AgamP4_X'], aes(x=pos, y=fst, color="#a6d854"))+
  theme_classic()+
  geom_line(size=0.8, colour='#a6d854')+
  annotate("text", x=15251572, y=1, label = 'Cyp9k1', angle = 0)+
  theme(legend.position = "none")+
  labs(x="Position", y='Fst', size = 12)+
  theme(axis.title.x = element_blank())+
  theme(text = element_text(size = 15))  +
  theme(axis.title.x = element_blank())

221306572-206055000

RColorBrewer::brewer.pal(name="Set2", n=5)
