#parse_pcangsd.py output
#run this in the same directory 
pkg = c('tidyverse', 'RcppCNPy', 'data.table', 'qqman' ,'cowplot', 'parallel', 'ragg', 'flextable', 'UpSetR')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

#set num cores
numcores = 10
#setwd
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_samples_downsampled/')

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
    ggtitle(names(prepped_admix_data[i]))+
    theme(
      legend.position = 'none',
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
      plot.title = element_text(size = 15, face = "bold"),
      strip.text.y = element_text(angle = 270, face = "bold"),
      strip.placement = "outside",
      panel.grid.major.y = element_blank())
}

#function to get p value from chisq dist selection stats
dochisq <- function(selection_array) {
  apply(selection_array, 1:2, pchisq, 1)
}

############################################################################
#plot admixture results from full set and downsampled core set

#get all admix qopt files
filelist = list.files(path = '.', pattern = '.Q') 
#read into list of qopt dfs
data_frame_list <- lapply(filelist, read.delim, header = F, sep = ' ')
#prep (sort) admix data for each df
prepped_admix_data  = mclapply(data_frame_list, admix_data_prep, mc.cores = numcores)
#get list of plotnames, remove some of filename junk
tmpnames = mclapply(X = filelist, FUN = function(t) gsub(pattern = "_glf2_minorfilts_angsd_samtools", replacement = "", x = t, fixed = TRUE), mc.cores = numcores)

s = mclapply(X = tmpnames, FUN = function(t) gsub(pattern = "downsampled", replacement = "", x = t, fixed = TRUE), mc.cores = numcores)

#horribly gsub out a load of rubbish from names
tmpnames = gsub(filelist, pattern = 'glf2_minorfilts_angsd_samtools_', replacement = "")
tmpnames = gsub(tmpnames, pattern = '_glf2_minorfilts_angsd_samtools', replacement = "")
tmpnames = gsub(tmpnames, pattern = '_downsampled', replacement = "")
tmpnames = gsub(tmpnames, pattern = '_core', replacement = "")
plotnames = gsub(tmpnames, pattern = 'admix.*', replacement = "")
#set names
prepped_admix_data = setNames(as.list(prepped_admix_data), plotnames)
#plot into list of plots
list_of_admix_plots = mclapply(seq_along(prepped_admix_data), plot_admix_data, mc.cores = numcores)
#setnames again
list_of_admix_plots = setNames(as.list(list_of_admix_plots), plotnames)



#plot chrs
AgamP4_2L = plot_grid(list_of_admix_plots[[7]],
                      list_of_admix_plots[[1]], 
                      list_of_admix_plots[[6]],
                      list_of_admix_plots[[5]],
                      list_of_admix_plots[[4]], 
                      list_of_admix_plots[[3]], 
                      list_of_admix_plots[[2]], 
                      ncols = 2)

tmpnames

ggsave('admix_AgamP4_2L_core_downsampled.pdf', plot = AgamP4_2L, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_samples_downsampled/')




AgamP4_2R = plot_grid(list_of_admix_plots[[14]],
                      list_of_admix_plots[[8]],
                      list_of_admix_plots[[13]],
                      list_of_admix_plots[[12]],
                      list_of_admix_plots[[11]], 
                      list_of_admix_plots[[10]], 
                      list_of_admix_plots[[9]], 
                      ncols = 2)
#AgamP4_2R
ggsave('admix_AgamP4_2R_core_downsampled.pdf', plot = AgamP4_2R, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_samples_downsampled/')

AgamP4_3L = plot_grid(list_of_admix_plots[[21]],
                      list_of_admix_plots[[15]], 
                      list_of_admix_plots[[20]],
                      list_of_admix_plots[[19]],
                      list_of_admix_plots[[18]],
                      list_of_admix_plots[[17]],
                      list_of_admix_plots[[16]],
                      ncols = 2)
#AgamP4_3L
ggsave('admix_AgamP4_3L_core_downsampled.pdf', plot = AgamP4_3L, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_samples_downsampled/')

AgamP4_3R  = plot_grid(list_of_admix_plots[[28]],
                       list_of_admix_plots[[22]],
                       list_of_admix_plots[[27]], 
                       list_of_admix_plots[[26]],
                       list_of_admix_plots[[25]], 
                       list_of_admix_plots[[24]],
                       list_of_admix_plots[[23]], 
                       ncols = 2)
#AgamP4_3R
ggsave('admix_AgamP4_3R_core_downsampled.pdf', plot = AgamP4_3R, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_samples_downsampled/')

AgamP4_X  = plot_grid(list_of_admix_plots[[29]],
                      list_of_admix_plots[[34]], 
                      list_of_admix_plots[[33]],
                      list_of_admix_plots[[32]], 
                      list_of_admix_plots[[31]],
                      list_of_admix_plots[[30]], 
                      ncols = 2)

ggsave('admix_AgamP4_X_core_downsampled.pdf', plot = AgamP4_X, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_samples_downsampled/')

############################################################################
#plot PCA results from full set and downsampled core set

#get covariance matrices
cofiles = list.files(path = '.', pattern = '.cov') 
covlist = mclapply(cofiles, read.delim, header = F, sep = ' ', mc.cores = numcores)
#eigenval and eigenvec compute of covariance matrices
eigenlist = mclapply(covlist, eigen, mc.cores = numcores)
#rmove rubbish from names
tmpnames = gsub(cofiles, pattern = "_minorfilts_angsd_samtools", replacement = "")
tmpnames = gsub(tmpnames, pattern = "_core", replacement = "")
plotnames = gsub(tmpnames, pattern = "\\.cov", replacement = "")
#set klist names to plot names
eigenlist = setNames(as.list(eigenlist), plotnames)

eigenlist[[1]]$vectors


#function to plot pcas
plotpca <- function(i) {
  ggplot(as.data.frame(eigenlist[[i]]$vectors), aes(x=V1, V2))+
    geom_point()+
    theme_minimal()+
    theme(plot.title = element_text(size = 7, face = "bold"))+
     ggtitle(plotnames[[i]])
}




#plot list of eigen dfs
pcaplots = mclapply(seq_along(eigenlist), plotpca, mc.cores = numcores * 1.5)


#plot chrs downsampled
pca_AgamP4_2L = plot_grid(pcaplots[[7]], 
                          pcaplots[[1]],
                          pcaplots[[6]], 
                          pcaplots[[5]], 
                          pcaplots[[4]], 
                          pcaplots[[3]],
                          pcaplots[[2]]
)
pca_AgamP4_2L
#pca_AgamP4_2L
ggsave('pca_AgamP4_2L.pdf', plot = pca_AgamP4_2L, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')

pca_AgamP4_2R = plot_grid(pcaplots[[14]], 
                          pcaplots[[8]],
                          pcaplots[[13]], 
                          pcaplots[[12]], 
                          pcaplots[[11]], 
                          pcaplots[[10]], 
                          pcaplots[[9]]
)
pcaplots[[14]]
ggsave('pca_AgamP4_2R.pdf', plot = pca_AgamP4_2R, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')

pca_AgamP4_3L = plot_grid(pcaplots[[21]], 
                          pcaplots[[15]],
                          pcaplots[[20]],
                          pcaplots[[19]], 
                          pcaplots[[18]], 
                          pcaplots[[17]], 
                          pcaplots[[16]])
pca_AgamP4_3L
ggsave('pca_AgamP4_3L.pdf', plot = pca_AgamP4_3L, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')

pca_AgamP4_3R  = plot_grid(pcaplots[[28]],
                           pcaplots[[22]],
                           pcaplots[[27]], 
                           pcaplots[[26]], 
                           pcaplots[[25]], 
                           pcaplots[[24]], 
                           pcaplots[[23]])
pca_AgamP4_3R   
ggsave('pca_AgamP4_3R.pdf', plot = pca_AgamP4_3R, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')

pca_AgamP4_X  = plot_grid( pcaplots[[29]],
                           pcaplots[[34]],
                           pcaplots[[33]],
                           pcaplots[[32]], 
                           pcaplots[[31]], 
                           pcaplots[[30]]
) 
#pca_AgamP4_X                         
ggsave('pca_AgamP4_X.pdf', plot = pca_AgamP4_X, device = 'pdf', path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')


#############################
#check out all the outliers
#after looking at all the plots, select outlier clusters and concat into a single df to look at later
outliersamples= rbind(
  list_of_admix_plots[[14]]$data %>% filter(main_proportion != 'V2' & main_proportion != 'V4') %>% select(id , main_proportion) %>% unique() %>% mutate(chrom = 'agamp4_2r'),
  list_of_admix_plots[[21]]$data %>% filter(main_proportion != 'V3' & main_proportion != 'V4') %>% select(id , main_proportion) %>% unique() %>% mutate(chrom = 'agamp4_3l'),
  list_of_admix_plots[[28]]$data %>% filter(main_proportion != 'V2' & main_proportion != 'V5') %>% select(id , main_proportion) %>% unique() %>% mutate(chrom = 'agamp4_3r')
)

list_of_admix_plots[[8]]


metadata$bam_key = as.numeric(metadata$bam_key)
outliersamples$id = as.numeric(outliersamples$id)
t = outliersamples %>% left_join(., metadata, by = c('id' = 'bamorder'))
s = t %>% select(chrom, ID, bam_key, Site)


#############################
#PLOT selection stats
list_of_pval_dfs = list()
#load files and calc p val for each site
selfiles = list.files(path = '.', pattern = '.selection') 
list_of_pval_dfs <- mclapply(mclapply(selfiles, RcppCNPy::npyLoad, mc.cores = numcores),dochisq, mc.cores = numcores)
list_of_pval_dfs <- mclapply(seq_along(list_of_pval_dfs), function(i){-log10(1-list_of_pval_dfs[[i]])}, mc.cores = numcores)
list_of_pval_dfs[[2]]

#remove rubbish from names and name df list
tmpnames = gsub(selfiles, pattern = "_minorfilts_angsd_samtools", replacement = "")
tmpnames = gsub(tmpnames, pattern = "_core", replacement = "")
plotnames = gsub(tmpnames, pattern = "\\.selection.npy", replacement = "")
list_of_pval_dfs = setNames(as.list(list_of_pval_dfs), plotnames)
#get all 'pos' files
posfiles = list.files(path = '.', pattern = '.pos.gz') 
#read corresponding pos file and bind 'pos' column to selection df list
list_of_pval_dfs = mclapply(seq_along(list_of_pval_dfs), function(i){cbind(fread(posfiles[[i]])[,2], list_of_pval_dfs[[i]])}, mc.cores = numcores)


test_dfs = lapply(list_of_pval_dfs, sample_n, 100000)

#plot pc1 of selection plots
plot_selection <- function(i, pcs) {
  pcnum = pcs
  pc_col = paste0("V", pcnum)
  pnam = paste0(plotnames[[i]], "_PC",as.character(pcnum))
  ggplot(test_dfs[[i]], aes(x=.data[["pos"]], y=.data[[pc_col]]))+
    geom_point()+
    theme_minimal()+
    labs(title = pnam)
}


test_dfs[[2]]
pc1plots = lapply(seq_along(test_dfs), plot_selection, 1)
pc1plots[[15]]
list_of_pval_dfs[[15]] %>% filter(V1 < 0.05) %>% count()

selection_pc2_2l <- plot_grid(pc1plots[[2]],
                              pc1plots[[3]],
                              pc1plots[[4]],
                              pc1plots[[5]],
                              pc1plots[[6]],
                              pc1plots[[1]],
                              ncol =1
)
selection_pc2_2l
save_plot('selection_pc1_2l.pdf', selection_pc1_2l, base_width = 10, base_height = 20, path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')

selection_pc1_2r <- plot_grid(pc1plots[[8]],
                              pc1plots[[9]],
                              pc1plots[[10]],
                              pc1plots[[11]],
                              pc1plots[[12]],
                              pc1plots[[7]],
                              ncol =1)
save_plot('selection_pc1_2r.pdf', selection_pc1_2r, base_width = 10, base_height = 20, path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')

selection_pc1_3l <- plot_grid(pc1plots[[14]],
                              pc1plots[[15]],
                              pc1plots[[16]],
                              pc1plots[[17]],
                              pc1plots[[18]],
                              pc1plots[[13]],
                              ncol =1)
save_plot('selection_pc1_3l.pdf', selection_pc1_3l, base_width = 10, base_height = 20, path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')

selection_pc1_3r <- plot_grid(pc1plots[[20]],
                              pc1plots[[21]],
                              pc1plots[[22]],
                              pc1plots[[23]],
                              pc1plots[[24]],
                              pc1plots[[19]],
                              ncol =1)
save_plot('selection_pc1_3r.pdf', selection_pc1_3r, base_width = 10, base_height = 20, path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')

selection_pc1_x <-  plot_grid(pc1plots[[26]],
                              pc1plots[[27]],
                              pc1plots[[28]],
                              pc1plots[[29]],
                              pc1plots[[30]],
                              pc1plots[[25]],
                              ncol =1)
save_plot('selection_pc1_x.pdf', selection_pc1_x, base_width = 10, base_height = 20, path = '~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_chr_downsampled/')


common_sites <- Reduce(intersect, list(test_dfs[[13]][, 1], test_dfs[[18]][, 1], test_dfs[[17]][, 1], test_dfs[[16]][, 1], test_dfs[[15]][, 1]));


selsites = rbind(
test_dfs[[18]]  %>% filter(pos %in% common_sites$pos) %>% select(pos, V1) %>% mutate(doc = 7.5),
test_dfs[[13]]  %>% filter(pos %in% common_sites$pos) %>% select(pos, V1) %>% mutate(doc = 10))

selsites

selsites %>% ggplot(aes(x=as.factor(doc), y=V1))+
  geom_point()+
  geom_line(aes(group = pos))

