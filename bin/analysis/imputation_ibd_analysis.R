#imputation and IBD analysis
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_chr_downsampled/')

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



x=list.files(pattern = '*.txt')
n = gsub("imputed_","",gsub(".fullset.*","",x))

x

imputedf = do.call(rbind, lapply(seq_along(x), function(i){fread(x[[i]]) %>% mutate(form = n[[i]])}))
colnames(imputedf) = c('chr', 'pos', 'R2', 'AF', 'IMP', 'form')
imputedf$IMP = gsub('\\.', '0', imputedf$IMP)


imputedf %>% group_by(form, chr, IMP) %>% count() %>% ggplot(aes(x=chr,y=n,fill=IMP))+geom_bar(stat='identity') +facet_wrap(~form)+theme_minimal()+labs(x='Chromosome', y='No. Sites')

imputedf %>% group_by(form, chr, IMP) %>% count()

simpute = imputedf %>% filter(form == 's_form' & IMP == 1) %>% sample_n(18650367/100)%>% ggplot(aes(x=R2, y=AF)) + facet_wrap(~chr)+theme_minimal()+geom_point()
mimpute = imputedf %>% filter(form == 'm_form' & IMP == 1) %>% sample_n(18650367/100)%>% ggplot(aes(x=R2, y=AF)) + facet_wrap(~chr)+theme_minimal()+geom_point()


plot_grid(simpute, mimpute, align = 'v',ncol=1,labels = c('S', 'M'))

imputedf %>% filter(R2 >= 0.4) %>% group_by(form, chr, IMP) %>% count()

#ok, now let's look at the downsampled data
rm(imputedf)
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_form_by_chr_downsampled/')
x = list.files(pattern = '_3R.bg5_imputestats.txt')
n = gsub(".bg5_imputestats.txt","", x)
n
impute_stats_downsampled_dfs = do.call(rbind, mclapply(seq_along(x), function(i){fread(x[[i]]) %>% mutate(form = n[[i]]) %>% sample_n(nrow(.)/100)}, mc.cores = 10))
impute_stats_downsampled_dfs = separate(impute_stats_downsampled_dfs, form, into = c('form', 'gub', 'cov', 'chr'), sep = '_')
colnames(impute_stats_downsampled_dfs) = c('chr', 'pos', 'R2', 'AF', 'IMP', 'form', 'gub', 'cov', 'chrom')
impute_stats_downsampled_dfs$cov = ordered(impute_stats_downsampled_dfs$cov, levels = c("full", "10", "7.5", "5", "2.5","1","0.5"))
impute_stats_downsampled_dfs %>% filter(form=='s' & IMP == 1) %>% ggplot(aes(x=R2,y=AF))+geom_violin()+facet_wrap(~cov)


#plot imputed summary data
impute_sum  = read.csv('imputation_results.csv')
colnames(impute_sum) = c('Form', 'Total Imputed Sites', 'Imputed Sites DR > 0.4', 'Original Sites', 'Coverage')
impute_sum$Coverage = ordered(impute_sum$Coverage, levels = c("0.5", "1", "2.5", "5", "7.5", "10", "full"))

plots = lapply(list('m_form','s_form'), function(i){
    impute_sum %>% pivot_longer(cols = 2:4, names_to = "Site Type") %>% 
    filter(Form == i) %>% 
    ggplot(aes(x=Coverage, y=value, fill=`Site Type`))+
    facet_wrap(~Coverage)+
    geom_line()
    theme_minimal()+
    labs(x='Site Type', y='Count')+
    theme(axis.text.x=element_blank())+
    scale_fill_brewer(palette = "Set2")
})


impute_sum %>% select(-`Total Imputed Sites`) %>% pivot_longer(cols = 2:3, names_to = "Site Type") %>% 
  ggplot(aes(x=Coverage, y=value, group=`Site Type`, colour=`Site Type`))+
  geom_line()+
  geom_point()+
  facet_wrap(~Form)+
  theme_minimal()+
  labs(x='Coverage', y='Count')+
  scale_color_brewer(palette = "Set1")


legend = get_legend(plots[[1]])
plot_grid((plots[[1]]+theme(legend.position="none")), (plots[[2]]+theme(legend.position="none")), legend, rel_widths = c(1,1, .5), labels=c('M','S'))

          