########
#parse_downsampled_ngsrelate_output


#parse_pcangsd.py output
#run this in the same directory 
pkg = c('tidyverse', 'RcppCNPy', 'data.table', 'qqman' ,'cowplot', 'parallel', 'ragg')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

#set num cores
numcores = 10
#setwd
setwd('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_samples_wg_downsampled/')
bamlist = fread('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/core_sample_bamlist', header = F, col.names = 'sample_id')
bamlist$samplenum = seq(0, length(bamlist$sample_id)-1)
res = fread('10bamlist.txtdownsampled.res')
res = left_join(res, bamlist, by = c('a' = 'samplenum')) %>% left_join(., bamlist, by = c('b' = 'samplenum'))
res = res %>% mutate(KING = replace(KING, KING < 0, 0)) %>% mutate(R0 = replace(R0, R0 < 0, 0)) %>% mutate(R1 = replace(R1, R1 < 0, 0))
kingr1 = res %>% filter(R1 < 1) %>% 
  ggplot(aes(x=R1, y=KING))+
  geom_point()+
  theme_minimal()
kingr1

resfiles = list.files(path = ".", pattern  = "res")
resdflist = lapply(resfiles, fread)
resdflist = setNames(resdflist, covvec)
covvec = c(10, 1, 2.5, 5, 7.5)

plotr1king <- function(i, bamlist) {
  resdf = resdflist[[i]] %>% mutate(KING = replace(KING, KING < 0, 0)) %>% mutate(R0 = replace(R0, R0 < 0, 0)) %>% mutate(R1 = replace(R1, R1 < 0, 0))
  resdf %>% filter(R1 < 1) %>% 
    ggplot(aes(x=R1, y=KING))+
    geom_point()+
    ggtitle(covvec[[i]])+
    theme_minimal()
}

plotr1r0 <- function(i, bamlist) {
  resdf = resdflist[[i]] %>% mutate(KING = replace(KING, KING < 0, 0)) %>% mutate(R0 = replace(R0, R0 < 0, 0)) %>% mutate(R1 = replace(R1, R1 < 0, 0))
  resdf %>% filter(R1 < 1) %>% 
    ggplot(aes(x=R1, y=R0))+
    geom_point()+
    ggtitle(covvec[[i]])+
    theme_minimal()
}

plotrelhist <- function(i, bamlist) {
  resdf = resdflist[[i]] %>% mutate(KING = replace(KING, KING < 0, 0)) %>% mutate(R0 = replace(R0, R0 < 0, 0)) %>% mutate(R1 = replace(R1, R1 < 0, 0))
  resdf %>% filter(R1 < 1) %>% 
    ggplot(aes(x=rab))+
    geom_histogram()+
    ggtitle(covvec[[i]])+
    theme_minimal()
}

plotkinghist <- function(i, bamlist) {
  resdf = resdflist[[i]] %>% mutate(KING = replace(KING, KING < 0, 0)) %>% mutate(R0 = replace(R0, R0 < 0, 0)) %>% mutate(R1 = replace(R1, R1 < 0, 0))
  resdf %>% filter(R1 < 1) %>% 
    ggplot(aes(x=KING))+
    geom_histogram()+
    ggtitle(covvec[[i]])+
    theme_minimal()
}

relhistlist = lapply(seq_along(resdflist), plotrelhist)
kinghistlist = lapply(seq_along(resdflist), plotkinghist)
r1r0list = lapply(seq_along(resdflist), plotr1r0)
kingr1list = lapply(seq_along(resdflist), plotr1king)

rabhist = plot_grid(relhistlist[[2]],
                    relhistlist[[5]],
                    relhistlist[[4]],
                    relhistlist[[3]],
                    relhistlist[[1]])
ggsave('rabhist.pdf', plot=last_plot())


kinghist =  plot_grid(kinghistlist[[2]],
                      kinghistlist[[5]],
                      kinghistlist[[4]],
                      kinghistlist[[3]],
                      kinghistlist[[1]])
ggsave('kinghist.pdf', plot=last_plot())

r1r0plots = plot_grid(
          r1r0list[[1]],
          r1r0list[[5]],
          r1r0list[[4]],
          r1r0list[[3]],
          r1r0list[[2]]
          )
r1r0plots

ggsave('r1r0plots.pdf', plot=last_plot())

kinggr1plots = plot_grid(
  kingr1list[[1]],
  kingr1list[[5]],
  kingr1list[[4]],
  kingr1list[[3]],
  kingr1list[[2]])
ggsave('kinggr1plots.pdf', plot=last_plot())
resfiles = list.files(path = ".", pattern  = "res")
resdflist = lapply(resfiles, fread)
resdflist = setNames(resdflist, covvec)
covvec = c(10, 1, 2.5, 5, 7.5)
resdflist = lapply(seq_along(resdflist), function(i){mutate(resdflist[[i]], meandoc = covvec[[i]])})
resdflist = lapply(resdflist, function(resdf) {mutate(resdf, relid = paste0(a,"_",b))})
relidlist = resdflist[[1]] %>%  filter(KING > 0.2) %>% select(relid)
downsampleinds = lapply(resdflist, function(resdf) {filter(resdf, relid %in% relidlist$relid)})
downsampleinds = do.call("rbind", downsampleinds)


a = downsampleinds %>% ggplot(aes(x=as.factor(meandoc), y=KING))+
  geom_point()+
  geom_line(aes(group=relid), alpha=0.5) +
  theme_classic()+
  labs(x='Mean DOC', y= 'KING kinship')

b = downsampleinds %>% ggplot(aes(x=as.factor(meandoc), y=R0))+
  geom_point()+
  geom_line(aes(group=relid), alpha=0.5) +
  theme_classic()+
  labs(x='Mean DOC', y= 'R0')

c = downsampleinds %>% ggplot(aes(x=as.factor(meandoc), y=R0))+
  geom_point()+
  geom_line(aes(group=relid), alpha=0.5) +
  theme_classic()+
  labs(x='Mean DOC', y= 'R1')

d = downsampleinds %>% ggplot(aes(x=as.factor(meandoc), y=rab))+
  geom_point()+
  geom_line(aes(group=relid), alpha=0.5) +
  theme_minimal()+
  labs(x='Mean DOC', y= 'rab')


reldownsample = plot_grid(a, b, c, align = 'hv')
reldownsample
ggsave('reldownsample.pdf')

filter(downsampleinds, KING < -0.50)
