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
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_samples_downsampled/')

#plot_admix_data

#function to order admix data and join metadata to it
admix_data_prep <- function(admixdf) {
  admixture_plot_data = as.data.frame(admixdf, .id="id") %>% 
    mutate(id = row_number()) %>%
    pivot_longer(cols=c(1:(length(admixdf))), names_to = 'pop', values_to = 'prob') %>% 
    select(pop, prob, id) %>% 
    group_by(id) %>% 
    mutate(main_proportion = pop[which.max(prob)], assignment_prob = max(prob)) %>% 
    arrange(main_proportion, desc(assignment_prob)) %>%  
    ungroup() %>% 
    mutate(id = forcats::fct_inorder(factor(id))) 
  return(admixture_plot_data)
}

#function to plot admix plots
plot_admix_data <- function(i) {
  ggplot(prepped_admix_data[[i]], aes(id, prob, fill = pop)) +
    geom_bar(position="fill", stat="identity", width = 1) +
    theme_minimal() +
    labs(x='Individual', y= 'Admixture Proportion')+
    facet_grid(~main_proportion, scales = 'free_x', space = "free_x")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(
      legend.position = 'none',
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
      plot.title = element_text(size = 15, face = "bold"),
      strip.text.y = element_text(angle = 270, face = "bold"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.placement = "outside",
      strip.text = element_blank(),
      panel.grid.major.y = element_blank())
}

#get all admix qopt files
filelist = list.files(path = '.', pattern = '*.Q') 
filelist = filelist[grep('3L', filelist)]
#read into list of qopt dfs
data_frame_list <- lapply(filelist, read.delim, header = F, sep = ' ')
#prep (sort) admix data for each df
prepped_admix_data  = mclapply(data_frame_list, admix_data_prep, mc.cores = numcores)
#get names
tmpnames = gsub("\\.admix.*","",filelist)

s = mclapply(X = tmpnames, FUN = function(t) gsub(pattern = "downsampled", replacement = "", x = t, fixed = TRUE), mc.cores = numcores)

#set names and plot
list_of_admix_plots = mclapply(seq_along(prepped_admix_data), plot_admix_data, mc.cores = numcores)
list_of_admix_plots = setNames(as.list(list_of_admix_plots), tmpnames)
covlabs = c('Full', '10', '7.5', '5', '2.5', '1', '0.5')
admixgrid = cowplot::plot_grid(list_of_admix_plots[[7]]+ggtitle('Full'),
                               list_of_admix_plots[[1]]+ggtitle('10'),
                               list_of_admix_plots[[6]]+ggtitle('7.5'),
                               list_of_admix_plots[[5]]+ggtitle('5'),
                               list_of_admix_plots[[4]]+ggtitle('2.5'),
                               list_of_admix_plots[[3]]+ggtitle('1'),
                               list_of_admix_plots[[2]]+ggtitle('0.5'))
admixgrid

#plot everything by chr and cov and pop
coverages = list('full', '10', '7.5', '5', '2.5', '1', '0.5')
populations = list('rainforest', 'mangrove', 'decid_forest', 'savannah', 's_form', 'm_form')
chrom=list('AgamP4_2L', 'AgamP4_2R', 'AgamP4_3L', 'AgamP4_3R', 'AgamP4_X')

#nested loop that iterates over populations, for each pop iterate over coverages and apply a function that extracts a list of plots per pop
#passes list to cowplot plotgrid, saves plot grid as pdf
for (population in populations) {
  for (coverage in coverages) {
    ggsave(
      plot_grid(plotlist = list_of_admix_plots[unlist(as.vector(lapply(chrom, function(d) {str_glue(population,".core.", coverage,".list.", d)})))]),
      filename = str_glue(population, "_by_chr_full_",coverage,"_admix.pdf"), device = 'pdf', path = '~/OneDrive - University of Glasgow/anopheles_lcwgs/data/by_pop_downsampling/admix/', width = 15, height = 10
    )
  }
}





############################################################################
#plot PCA results from full set and downsampled core set

#get covariance matrices, filter out 'tree.cov
yesfiles = list.files(path = '.', pattern = '.cov')
nofiles = list.files(path='.', pattern = 'tree.cov')
cofiles = setdiff(yesfiles, nofiles)
cofiles = cofiles[grep('3L', cofiles)]

#read everything
covlist = mclapply(cofiles, read.delim, header = F, sep = ' ', mc.cores = numcores)
#eigenval and eigenvec compute of covariance matrices
eigenlist = mclapply(covlist, eigen, mc.cores = numcores)
#rmove rubbish from names
tmpnames = gsub(cofiles, pattern = "_minorfilts_angsd_samtools", replacement = "")
tmpnames = gsub(tmpnames, pattern = "_core", replacement = "")
pcaplotnames = gsub(tmpnames, pattern = "\\.cov", replacement = "")
#set klist names to plot names
eigenlist = setNames(as.list(eigenlist), pcaplotnames)
#coerce covariance matrix to df with first two pcs
eigenframes = lapply(seq_along(eigenlist), function(i){as.data.frame(eigenlist[[i]]$vectors)})
#add column with bam number as data key
eigenframes2 = lapply(seq_along(eigenframes), function(i){eigenframes[[i]] %>% mutate(bam_key = seq(1, nrow(eigenframes[[i]])))})
#add ecozone as column by subbing out rubbish from dataframe title
eigenframes3 = lapply(seq_along(eigenframes2), function(i){cbind(eigenframes2[[i]],population=gsub("\\..*","", pcaplotnames[[i]]))})
#make list of pop metadatas
metadatalist = list(decid_forest = decid_forest, s_form = sform, m_form = mform, savannah = savannah, rainforest=rainforest, mangrove=mangrove)
##join to corresponding dataframe

#neweigenframes = lapply(seq_along(eigenframes3), function(i) {left_join(metadatalist[[unique(eigenframes3[[i]]$population)]], eigenframes3[[i]], by = c('bamorder' = 'bam_key'))})
metadata = read.csv('~/Projects//MOVE/td_je_angam_2022/metadata/metadata_with_insecticide_info.csv')
corebams = fread('~/Projects//MOVE/data/anopheles_03_21_DW/metadata/core_bams.list', header=F)
coremeta = metadata[metadata$seq_id %in% corebams$V1 ,]

neweigenframes = lapply(seq_along(eigenframes3), function(i) {left_join(coremeta, eigenframes3[[i]], by = c('bamorder' = 'bam_key'))})



s = left_join(eigenframes3[[1]], coremeta, by =c('bam_key'= 'bamorder'))
s

ggplot(s, aes(x=V1, y=V2, colour=Form))+geom_point()
t
#function to plot pcas
plotpca <- function(i, colourfactor) {
  ggplot(as.data.frame(neweigenframes[[i]]), aes(x=V1, y=V2))+
    geom_point()+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
    #theme(legend.position="none")+
}





#plot list of eigen dfs
pcaplots = mclapply(seq_along(neweigenframes), plotpca, habitat, mc.cores = numcores * 1.5)

#setnames
pcaplots = setNames(as.list(pcaplots), pcaplotnames)


cowplot::plot_grid( pcaplots[[7]]+ggtitle('Full'),
                    pcaplots[[1]]+ggtitle('10'),
                    pcaplots[[6]]+ggtitle('7.5'),
                    pcaplots[[5]]+ggtitle('5'),
                    pcaplots[[4]]+ggtitle('2.5'),
                    pcaplots[[3]]+ggtitle('1'),
                    pcaplots[[2]]+ggtitle('0.5'))



#nested loop that iterates over populations, for each pop iterate over coverages and apply a function that extracts a list of plots per pop
#passes list to cowplot plotgrid, saves plot grid as pdf
for (population in populations) {
  for (coverage in coverages) {
    ggsave(
      plot_grid(plotlist = pcaplots[unlist(as.vector(lapply(chrom, function(d) {str_glue(population,".core.", coverage,".list.", d)})))], legend = get_legend(baseplot)),
      filename = str_glue(population, "_by_chr_full_",coverage,"_pca.pdf"), device = 'pdf', path = '~/OneDrive - University of Glasgow/anopheles_lcwgs/data/by_pop_downsampling/full_sample_set_fullcov/', width = 15, height = 10
    )
  }
}


setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_chr_downsampled')
fstlist = list.files(pattern = 'window5step1')
fstlist = fstlist[grep(pattern = 'X',fstlist)]

readfst <- function(x) {
  fil = fread(fstlist[[x]], header=F)
  setnames(fil, c('region', 'chr', 'midpoint', 'numsites', 'fst'))
}
#create plot labels from filenames
plotlabs = gsub("m_form_s_form", "mform_sform", fstlist)
plotlabs = stringr::str_extract(plotlabs, "[^_]*_[^_]*")
plotlabs
#read in our shit
s=mclapply(seq_along(fstlist), readfst, mc.cores = 8)
#add 'pop' row in df for ggplot colouring
s= mclapply(seq_along(s), function(i){s[[i]] %>% mutate(pop = plotlabs[[i]])}, mc.cores=10)
s= setNames(as.list(s), fstlist)
#make plots
fstplots = mclapply(seq_along(s), mc.cores = 6, function(i)
  {
    ggplot(s[[i]], aes(x=midpoint, y=fst, fill=as.factor(pop)))+
    scale_color_brewer(palette = "Set2")+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())+
    geom_line(color='#fde725')+
    ylim(0,1)})
  


thetalist = list.files(pattern = '.windowout.pestPG')
thetalist = thetalist[grep(pattern = 'X',thetalist)]
thetalist = thetalist[grep(pattern = 'm_form',thetalist)]
thetafiles = lapply(thetalist, fread)
thetafiles = setNames(thetafiles, thetalist)
thetanames = gsub(".windowout.pestPG", "", thetalist)
thetanames = gsub("m_form_s_form", "mform_sform", thetanames)
pops = stringr::str_extract(thetanames, "[^_]*")

tajimaplots = lapply(seq_along(thetafiles), function(i){
  ggplot(thetafiles[[i]], aes(WinCenter, y=Tajima))+
    geom_line()+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  })
tajimaplots

piplots = lapply(seq_along(thetafiles), function(i){
  ggplot(thetafiles[[i]], aes(WinCenter, y=tP))+
    geom_line(color='#440154')+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
})


pigrid = plot_grid(piplots[[7]]+ggtitle(''),
                   piplots[[3]]+ggtitle(''),
                   piplots[[6]]+ggtitle(''),
                   piplots[[5]]+ggtitle(''),
                   piplots[[4]]+ggtitle(''),
                   piplots[[2]]+ggtitle(''),
                   piplots[[1]]+ggtitle(''),
                   ncol=1)

fstgrid = plot_grid(fstplots[[7]]+ggtitle('Full'),
                    fstplots[[3]]+ggtitle('10'),
                    fstplots[[6]]+ggtitle('7.5'),
                    fstplots[[5]]+ggtitle('5'),
                    fstplots[[4]]+ggtitle('2.5'),
                    fstplots[[2]]+ggtitle('1'),
                    fstplots[[1]]+ggtitle('0.5'),
                    ncol=1)


fstgrid
pigrid
plot_grid(fstgrid, pigrid, align = 'v', labels = c('A', 'B'), hjust = -25)


resfiles = list.files(path = ".", pattern  = "res")
resdflist = lapply(resfiles, fread)
resdflist = setNames(resdflist, covvec)
covvec = c(10, 1, 2.5, 5, 7.5)
resdflist = lapply(seq_along(resdflist), function(i){mutate(resdflist[[i]], meandoc = covvec[[i]])})
resdflist = lapply(resdflist, function(resdf) {mutate(resdf, relid = paste0(a,"_",b))})
relidlist = resdflist[[1]] %>%  filter(KING > 0.2) %>% select(relid)
downsampleinds = lapply(resdflist, function(resdf) {filter(resdf, relid %in% relidlist$relid)})
downsampleinds = do.call("rbind", downsampleinds)
a
a = downsampleinds %>% ggplot(aes(x=as.factor(meandoc), y=KING))+
  geom_point()+
  geom_line(aes(group=relid), alpha=0.5) +
  theme(axis.title.x = element_blank())+
  theme_classic()+
  labs(x='Mean DOC', y= 'KING kinship')

b = downsampleinds %>% ggplot(aes(x=as.factor(meandoc), y=R0))+
  geom_point()+
  geom_line(aes(group=relid), alpha=0.5) +
  theme(axis.title.x = element_blank())+
  theme_classic()+
  labs(x='Mean DOC', y= 'R0')

c = downsampleinds %>% ggplot(aes(x=as.factor(meandoc), y=R0))+
  geom_point()+
  geom_line(aes(group=relid), alpha=0.5) +
  theme(axis.title.x = element_blank())+
  theme_classic()+
  labs(x='Mean DOC', y= 'R1')

d = downsampleinds %>% ggplot(aes(x=as.factor(meandoc), y=rab))+
  geom_point()+
  theme(axis.text.x  = element_blank())+
  geom_line(aes(group=relid), alpha=0.5) +
  theme_minimal()+
  labs(x='Mean DOC', y= 'rab')
d

reldownsample = plot_grid(a, b, c, align = 'hv')
reldownsample




