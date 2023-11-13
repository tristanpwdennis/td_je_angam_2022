pkg = c("tidyverse", 
        "rnaturalearth", 
        "rnaturalearthdata", 
        "sf", 
        "cowplot", 
        "parallel", 
        "gridExtra", 
        "grid", 
        "data.table",
        'phangorn', 
        'phytools', 
        'pegas',
        'vcfR', 
        'plotly',
        'RColorBrewer', 
        'rehh')

#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
knitr::opts_chunk$set(echo = TRUE)
metadata<- fread('/Users/dennistpw/Library/Mobile Documents/com~apple~CloudDocs/Projects/td_je_angam/td_je_angam_2022/metadata/metadata_full_sequencedonly_062022.csv', header=T)


#make plottable factor levels
metadata$Species = ifelse(metadata$Form=='M', 'A. coluzzi', 'A. gambiae')
metadata$plothab <- case_when(
          metadata$habitat == 'rainforest'~ 'RF',
          metadata$habitat == 'deciduous_forest'~ 'DF',
          metadata$habitat == 'mangrove_swamp' ~ 'MS',
          metadata$habitat == 'coastal_savannah' ~ 'CS'
          )

#### nice map figure

#load geog objects
world <- ne_countries(scale = 'large', returnclass = "sf", continent = 'africa')
lakes <- ne_download(scale = 'large', type = 'lakes', category = 'physical', returnclass = "sf")
rivers <- ne_download(scale = 'large', type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
ocean <- ne_download(scale = 'large', type = 'ocean', category = 'physical', returnclass = "sf")
#get site only data and coilour by habitat (darken for border)
sites = metadata %>% select(Site, Long, Lat, habitat) %>% unique() #%>% full_join(colesites)
#darkerhabitatcolours3 = colorspace::darken(habitatocolours, amount=0.3)
#plot sites coloured by isec and eco, and species 
maps = ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  geom_point(data=metadata, aes(y=Lat,x=Long, color=plothab, fill=plothab), size = 4, alpha=0.7)+
  scale_color_brewer(palette = 'Dark2')+
  scale_fill_brewer(palette = 'Dark2')+
  #ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-3.6, 1.5), ylim = c(4.5, 11.5), expand = FALSE)+
  labs(x='Longitude', y='Latitude', color='Ecozone',fill = 'Ecozone', shape='Insecticide Usage Pattern')+
  #guides(fill=guide_legend(override.aes=list(colour=habitatocolours)))+
  theme_classic()+
  facet_wrap(~Species)
#
ggsave(plot=maps,filename = '~/Projects/td_je_angam_2022/figures/fig1_maps.tiff', device = grDevices::tiff, width = 30, height = 15, units = 'cm', dpi = 1000)

##pca and admix for both species
#fix spcies name
metadata$Species <- gsub('A. gambiae','A. gambiae s.s.',metadata$Species)

bothspecies <- fread('~/Projects/td_je_angam_2022/data/structure/3_euch_pruned.cov')
pca = eigen(bothspecies)
pca$values
pca_bothtaxa <- cbind(pca$vectors, metadata)

pcaplotbothspecies <- ggplot(pca_bothtaxa,aes(x=V1, y=V2, colour=Species))+
  geom_point(size=3, alpha=0.6)+
  scale_colour_manual(values=c('#1f78b4','#e31a1c'))+
  theme_classic()+
  labs(x='PC1 (20.23%)', y='PC2 (5.73%)', colour = 'Species')

#--------------------------------------------------#
#pca                                               #
#--------------------------------------------------#
mf = 'M'
mp = '~/Projects/td_je_angam_2022/data/pca/ancol.3L.pcangsd.cov'

sf = 'S'
sp = '~/Projects/td_je_angam_2022/data/pca/angam.3L.pcangsd.cov'

#function for plotting pcas

plot_pca <- function(form, path, pc1, pc2) {
  x = fread(path)
  pca = eigen(x)
  me = metadata[metadata$Form == form,]
  me$plothab <- factor(me$plothab, levels=c('DF','RF','CS', 'MS'))
  
  f_me <- cbind(pca$vectors, me)
  #how much variation explained by each PC?
  varexp <- pca$values / sum(pca$values)
  #plot pca
  
  
  pc1_2 <- ggplot(f_me,aes(x=V1, y=V2, colour=plothab))+
    geom_point(size=2, alpha=0.6)+
    scale_colour_brewer(palette='Dark2')+
    theme_classic()+
    labs(x='PC1', y='PC2', colour = 'Ecozone')+
    theme(legend.position = 'none')
  
  pc3_4 <- ggplot(f_me,aes(x=V3, y=V4, colour=plothab))+
    geom_point(size=2, alpha=0.6)+
    scale_colour_brewer(palette='Dark2')+
    theme_classic()+
    labs(x='PC3', y='PC4', colour = 'Ecozone')
  
  
  legend <- get_legend(pc3_4)
  pcas <- plot_grid(pc1_2, pc3_4+theme(legend.position = 'none'), labels = c('A','B'))
  return(pcas)
  #return(list(legend, pcas))
}



gamK = 3
gampath = '~/Projects/td_je_angam_2022/data/pca/angam.3L.pcangsd.admix.3.Q'
gampal = c('#7fc97f', '#beaed4', '#fdc086')

colK = 2
colpath = '~/Projects/td_je_angam_2022/data/pca/ancol.3L.pcangsd.admix.2.Q'
colpal = c('#ffff99', '#386cb0')

plot_admix <- function(K, path, palcols) {
  
  q = fread(path)
  
  q$id <- row.names(q)
  metadata$id = seq(1: nrow(metadata))
  metadata$id = as.numeric(metadata$id)
  q$id = as.numeric(q$id)
  popqmat = left_join(q, metadata)
  
  #get max index
  # Function to find the index of the maximum value in a row
  find_max_index <- function(row) {
    return(which.max(row[1:K]))
  }
  
  #get max prob per row:
  popqmat$max_prob <- apply(popqmat, 1, find_max_index)
  
  long_sort = pivot_longer(popqmat, cols = sprintf("V%s",seq(1:K)))
  
  long_sort %>% mutate(id = forcats::fct_reorder(as.factor(id), as.numeric(max_prob))) %>% 
    ggplot(., aes(factor(id), value, fill = factor(name))) +
    geom_col(color = "gray28", size = 0.1) +
    facet_grid(~plothab, switch = "x", scales = "free", space = "free") +
    theme_minimal() + labs(x = NULL, title = paste0("K=",K), y = "Ancestry") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expand_scale(add = 1)) +
    scale_fill_manual(values=palcols)+
    theme(
      axis.text.x = element_blank(),
      panel.spacing.x = unit(0.1, "lines"),
      panel.grid = element_blank(),
      legend.position = "none",
    ) 
  
}

col_pca <-plot_pca(mf, mp)
admixcol <- plot_admix(colK, colpath, colpal)


gam_pca <- plot_pca(sf, sp)
admixgam <- plot_admix(gamK, gampath, gampal)


colpcaplot <- col_pca
plot_grid(col_pca, admixcol, ncol=1, nrow=2, labels = c('','C'))
plot_grid(gam_pca, admixgam, ncol=1, nrow=2, labels = c('','C'))

