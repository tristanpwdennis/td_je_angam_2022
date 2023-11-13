#####genomescans
library(tidyverse) #tidyverse 
library(data.table) #datatable
library(RColorBrewer) #nice colours
library(kableExtra)
library(parallel)

####
#fst and dxy

#read fst data
fstfiles <- list.files(path='~/Projects/td_je_angam_2022/data/fst_by_ecotype_species/', pattern='*win10kstp5k*', full.names = TRUE)
fstdata  <- mclapply(fstfiles, function(x) data.table(fread(x), 'name' = basename(x))) #read data
fstdata <- do.call(rbind, fstdata)
colnames(fstdata) <- c('bit','chrom','midpos','nSites','fst','name')

#sort comp levels for fsts
fstdata$name <-gsub('form','',fstdata$name)
fstdata <- fstdata %>% separate(name, c('species_a','ecozone_a','species_b', 'ecozone_b', 'chrpref','chrom','suffix'))

#read dxy data
dxyfiles <- list.files(path='~/Projects/td_je_angam_2022/data/dxy/', pattern='*windowed.txt', full.names = TRUE)
dxydata  <- mclapply(dxyfiles, function(x) data.table(fread(x), 'name' = basename(x))) #read data
dxydata <- do.call(rbind, dxydata)
colnames(dxydata) <- c('bit','chrom','midpos','dxy','name')

#sort comp levels for fsts
dxydata$name <-gsub('.dxy.windowed.txt','',dxydata$name)
dxydata$name <-gsub('decid_forest','decidforest',dxydata$name)
unique(dxydata$name)
dxydata <- dxydata %>% separate(name, c('species_a','ecozone_a','species_b', 'ecozone_b'))


#fiuddle with names for df merge
dxydatasub <- dxydata[c('chrom','midpos','dxy','species_a','species_b', 'ecozone_a','ecozone_b')]
colnames(dxydatasub) <- c('chrom','midpos','stat','species_a','species_b', 'ecozone_a','ecozone_b')
fstdatasub <- fstdata[c('chrom','midpos','fst','species_a','species_b', 'ecozone_a','ecozone_b')]
colnames(fstdatasub) <- c('chrom','midpos','stat','species_a','species_b', 'ecozone_a','ecozone_b')
mean(fstdata$nSites)

#need ro swap values for ecozone a and b for species comp gam and col
# Condition to select rows for which you want to swap values
condition <- fstdatasub$species_a == 'm' & fstdatasub$species_b == 's'

# Create a temporary variable to hold col1 values
temp_values <- fstdatasub$ecozone_a[condition]

# Swap values in col1 and col2 for rows where the condition is met
fstdatasub$ecozone_a[condition] <- fstdatasub$ecozone_b[condition]
fstdatasub$ecozone_b[condition] <- temp_values

#add stat type to dxy and fst
dxydatasub$stat_type <- 'dxy'
fstdatasub$stat_type <- 'fst'

#temporary keep only gam
#bind
specstats <- rbind(dxydatasub, fstdatasub)

#chrom id tidy
specstats$chrom <- gsub('AgamP4_','',specstats$chrom) #tidy chrid


#plottable ecozone names
specstats$ecozone_a_s <- case_when(
  specstats$ecozone_a == 'rainforest'~ 'RF',
  specstats$ecozone_a == 'decidforest'~ 'DF',
  specstats$ecozone_a == 'mangrove' ~ 'MS',
  specstats$ecozone_a == 'savannah' ~ 'CS'
)

#plottable ecozone names
specstats$ecozone_b_s <- case_when(
  specstats$ecozone_b == 'rainforest'~ 'RF',
  specstats$ecozone_b == 'decidforest'~ 'DF',
  specstats$ecozone_b == 'mangrove' ~ 'MS',
  specstats$ecozone_b == 'savannah' ~ 'CS'
)

#plottable species names
specstats$species_a <- case_when(
  specstats$species_a == 'm'~ 'An. col',
  specstats$species_a == 's'~ 'An. gam s.s.',
)
specstats$species_b <- case_when(
  specstats$species_b == 'm'~ 'An. col',
  specstats$species_b == 's'~ 'An. gam s.s.',
)

#chrom ordering
specstats$chrom <- factor(specstats$chrom, levels=c('2R','2L','3R','3L','X'))

#spec_comp 
specstats$spec_comp <- paste0(specstats$species_a, ':', specstats$species_b)
#habitat_comp 
specstats$hab_comp <- paste0(specstats$ecozone_a_s, ':', specstats$ecozone_b_s)

#################################################
#identify outliers in data
#################################################

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


get_outlier_threshold <- function(column){
  outlierserthreshold <- getmode(column) + (3*(getmode(column) - min(column)))
}

specstats <- specstats %>% 
  filter(stat_type == 'fst') %>% 
  group_by(spec_comp, hab_comp) %>% 
  mutate(modal_value = get_outlier_threshold(stat),
         new_column = ifelse(stat > modal_value, "Outlier", "NotOutlier"))

#################################################
#plot fsts
#################################################

m_fsts <- specstats %>% filter(stat_type=='fst' & hab_comp != 'CS:DF' & hab_comp != 'RF:CS' & hab_comp != 'RF:DF' & species_a == 'An. col' & species_b == 'An. col') %>% 
  ggplot(aes(x=midpos,y=stat, colour=new_column))+
  geom_point(alpha=0.7, size=0.2)+
    scale_y_continuous(n.breaks=3)+
    scale_color_brewer(palette = 'Set2')+
    facet_grid(rows=vars(hab_comp),cols=vars(chrom),scales = "free_x", space = "free_x")+
    theme_classic()+
    labs(x='Position',y='Fst', colour=c('Species:Species'))+
    theme(
    #strip.text.y = element_blank(),
    axis.line.x = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.1, 0, 0, 0), "cm"))

s_fsts <- specstats %>% filter(stat_type=='fst' & hab_comp != 'CS:DF' & hab_comp != 'RF:CS' & hab_comp != 'RF:DF' & species_a != 'An. col' & species_b != 'An. col') %>% 
  ggplot(aes(x=midpos,y=stat, colour=new_column))+
  geom_point(alpha=0.7, size=0.2)+
  scale_y_continuous(n.breaks=3)+
  scale_color_brewer(palette = 'Set2')+
  facet_grid(rows=vars(hab_comp),cols=vars(chrom),scales = "free_x", space = "free_x")+
  theme_classic()+
  labs(x='Position',y='Fst', colour=c('Species:Species'))+
  theme(
    #strip.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    plot.margin = unit(c(0.1, 0, 0, 0), "cm"))


#################################################
#plot thetas
#################################################

thetas <- list.files(path='~/Projects/td_je_angam_2022/data/thetas/', pattern='.*indow*', full.names = TRUE)
ecozones <- gsub('_nomaffilter.thetasWindow.gz.pestPG','',basename(thetas))
#function plot thetas per ecozone and species
thetafilelist <- mclapply(thetas, function(x) data.table(fread(x), 'ecozone' = x)) #read data
thetafilelist <- mclapply(thetafilelist, function(x) x[x$Chr %like% "AgamP4_2" | x$Chr %like% "AgamP4_3" | x$Chr %like% "AgamP4_X", ]) #subset to chroms of interest
thetafilelist <- do.call(rbind, thetafilelist)


thetafilelist$ecozone <- gsub('_nomaffilter.thetasWindow.gz.pestPG','',basename(thetafilelist$ecozone))

thetafilelist$ecozone <- gsub('_nomaf.*','',basename(thetafilelist$ecozone))

thetafilelist$ecozone <- gsub('decid_forest','decidforest',basename(thetafilelist$ecozone))
thetafilelist$species <- gsub('_.*','',basename(thetafilelist$ecozone))
thetafilelist$ecozone <- gsub('*._','',basename(thetafilelist$ecozone))
thetafilelist$ecozone <- gsub('*._','',basename(thetafilelist$ecozone))

#plottable ecozone names
thetafilelist$ecozone <- case_when(
  thetafilelist$ecozone == 'rainforest'~ 'RF',
  thetafilelist$ecozone == 'decidforest'~ 'DF',
  thetafilelist$ecozone == 'mangrove' ~ 'MS',
  thetafilelist$ecozone == 'savannah' ~ 'CS'
)

thetafilelist$species <- case_when(
  thetafilelist$species == 'm'~ 'An. col',
  thetafilelist$species == 's'~ 'An. gam s.s.',
)

  #order ecozones and chgroms
thetafilelist$ecozone <- factor(thetafilelist$ecozone, levels=c('DF','RF','CS', 'MS'))
thetafilelist$Chr <- factor(thetafilelist$Chr, levels=c('AgamP4_2R','AgamP4_2L','AgamP4_3R','AgamP4_3L','AgamP4_X'))
thetafilelist <- thetafilelist[thetafilelist$Chr != 'Chr'] #remove accidential headers lol
  
  #make coords for plotting x axis facet sizes (so x axes limits per chrom facet are lengths proportionate to actual sizes)
maxchromlengths <- aggregate(WinCenter ~ Chr, thetafilelist, max)

#calculate pi
thetafilelist$Pi <- thetafilelist$tP / thetafilelist$nSites

m_thetas <- thetafilelist %>% filter(species == 'An. col') %>% 
  pivot_longer(cols=c('Tajima','Pi')) %>% 
  ggplot(aes(x=WinCenter,y=value,colour=ecozone))+
  geom_line(alpha=0.5, linewidth=0.5)+
  scale_x_continuous(n.breaks=3)+
    facet_grid(rows=vars(species,name), cols=vars(Chr), scales='free', space = "free_x")+
    theme_classic()+
  labs(x='Position (bp)',y='Stat', colour='Ecoregion')+
  scale_color_brewer(palette = "Set2")+
  theme(strip.text.x = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.7))

thetafilelist$theta <- thetafilelist$tW / thetafilelist$nSites
s_thetas <- thetafilelist %>% filter(species != 'An. col') %>% 
  pivot_longer(cols=c('Tajima','Pi')) %>% 
  ggplot(aes(x=WinCenter,y=value,colour=ecozone))+
  geom_line(alpha=0.5, linewidth=0.5)+
  labs(x='Position (bp)',y='Stat', colour='Ecoregion')+
  scale_x_continuous(n.breaks=3)+
  facet_grid(rows=vars(species,name), cols=vars(Chr), scales='free', space = "free_x")+
  theme_classic()+
  scale_color_brewer(palette = "Set2")+
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.7))
s_thetas
mthetaslegend <- cowplot::get_legend(m_thetas)
sthetaslegend <- cowplot::get_legend(s_thetas)
plot(mthetaslegend)
mplots <- cowplot::plot_grid(m_fsts, m_thetas+theme(legend.position = "none"),nrow=2, ncol=1,align = 'v', axis='bt', labels = c('A','B'))
ggsave('~/Projects/td_je_angam_2022/figures/fst_coluzzi_byecology.tiff', mplots, width = 18, height = 20, units = 'cm')

splots <- cowplot::plot_grid(s_fsts, s_thetas+theme(legend.position = "none"),nrow=2, ncol=1,align = 'v', axis='bt', labels = c('A','B'))
ggsave('~/Projects/td_je_angam_2022/figures/fst_gambiae_byecology.tiff', splots, width = 18, height = 15, units = 'cm')


####################################################################################################
#now let's make nice thetas table to summarise data
####################################################################################################


thetafilelist$pi <- thetafilelist$tP/thetafilelist$nSites
thetafilelist$Ne <- thetafilelist$tW / 3.5e-09 / thetafilelist$nSites / 4 # calculate effective population size
thetatab <- thetafilelist[,c('Chr','ecozone', 'species', 'Tajima','pi', 'Ne')]
thetatab <- thetatab %>% group_by(ecozone,species) %>% summarise(π=mean(pi), Tajima = mean(Tajima), Ne = mean(Ne))
thetatab$species <- case_when(
  thetatab$species == 'm'~ 'An. col',
  thetatab$species == 's'~ 'An. gam s.s.',
)

colnames(thetatab) <- c('Ecozone','Species','π','Tajima\'s D','Ne')
#make dataframe of syumamry stats
thetatab %>% arrange(Ecozone) %>% 
  kable("html") %>% 
  remove_column(1) %>% 
  kable_styling("striped") %>%
  column_spec(1, bold = TRUE, italic = TRUE,border_right = TRUE) %>%
  pack_rows("DF", 1,2) %>% 
  pack_rows("RF", 3,4) %>% 
  pack_rows("CS", 5,6) %>% 
  pack_rows("MS", 7,7)

####################################################################################################
#now let's make nice fst table to summarise fsts
####################################################################################################

specstats %>% filter(stat_type=='fst' & hab_comp != 'CS:DF' & hab_comp != 'RF:CS' & hab_comp != 'RF:DF') %>% 
 # filter(species_a =='An. col' & species_b =='An. col') %>% 
  group_by(species_a, species_b, ecozone_a_s, ecozone_b_s) %>% summarise(Fst = mean(stat)) 

specstats %>% filter(stat_type=='fst' & hab_comp != 'CS:DF' & hab_comp != 'RF:CS' & hab_comp != 'RF:DF') %>% 
   filter(species_a !='An. col' & species_b !='An. col') %>% 
  group_by(ecozone_a_s, ecozone_b_s) %>% summarise(Fst = mean(stat)) %>% 
  pivot_wider(names_from = ecozone_a_s, values_from =Fst) %>% 
  kable("html") %>% 
  column_spec(1, bold = TRUE,border_right = TRUE) %>%
  kable_styling("striped") %>% 
  add_header_above(c("Ecozone A" = 1,"Ecozone B"=2))
?add_header_above  


####################################################################################################
#intersect outlier fst windows with gene data for supplementary tables
####################################################################################################
#prepare bed format 
s = winfst %>% filter(outlier == 'outlier') #get only outliers > 99.9th %ile
s = mutate(s, start = pos -500, end = pos + 500) #make start stop coords for bed file (from window midpoint)
s = select(s, chr, start, end, piece, fst, form, id)#select bed cols
s = s[order(chr, start, end)] #order our windows


beddf <- specstats %>% filter(stat_type=='fst' & new_column == 'Outlier' & hab_comp != 'CS:DF' & hab_comp != 'RF:CS' & hab_comp != 'RF:DF' & spec_comp != 'An. col:An. gam s.s.') %>% 
  mutate(start = midpos - 5000) %>% 
  mutate(end = midpos + 5000) %>% 
  dplyr::select(chrom, start, end, spec_comp, hab_comp, new_column, stat) %>% 
  arrange(chrom, start, end)

##read annotations and intersect
genes = fread('grep protein_coding_gene ~/Projects/td_je_angam_2022/annotations/files/VectorBase-55_AgambiaePEST.gff') #read gff file
genebed = dplyr::select(genes, V1, V4,V5, V9) #convert to bed (see UCSC documentation)
genebed = genebed[order(V1, V4),] #sort just in case
colnames(genebed) = c('chrom', 'start', 'end', 'name') #add colnames
genebed$chrom <- gsub('AgamP4_' , '',genebed$chrom)

outlierwindow_genes = tidygenomics::genome_intersect(beddf, genebed, c("chrom"="chrom", "start"="start", "end"="end")) #intersect
write_csv(outlierwindow_genes, '~/Projects/td_je_angam_2022/data/annotation/outlier_window_genes.csv')



outlierwindow_genes = data.table(outlierwindow_genes)

outlierwindow_genes %>% filter(chrom == '2L' & start > 20528089 & end <42165182)

max_windows = outlierwindow_genes[outlierwindow_genes[, .I[which.max(stat)], by=list(name)]$chrom]


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
