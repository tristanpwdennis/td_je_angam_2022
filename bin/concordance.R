gatk_concordance = fread('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/concordance_results.txt')
gatk_concordance

gatk_concordance$V1 = gsub('.txt:SNP', '', gatk_concordance$V1)
gatk_concordance$V1 = gsub(".*/","",gatk_concordance$V1)
gatk_concordance = separate(gatk_concordance, V1, into = c('ind', 'cov'), sep = '_')
gatk_concordance = separate(gatk_concordance, cov, into = c('Coverage', 'Filtered or Unfiltered'), sep = '-')

colnames(gatk_concordance) = c('ind', 'Coverage', 'Filtering', 'TP', 'FP', 'FN', 'Precision', 'Recall')


t = gatk_concordance %>% pivot_longer(colsß = 7:8) 
t
t %>% filter(name=='Precision' | name == 'Recall') %>% ggplot(., aes(x = Coverage,y=value, color=name))+geom_line()+theme_minimal()+labs(y='Value',color='Precision and Recall')+geom_point()
gatk_concordance %>% ggplot(aes(x=Precision, y=Recall, color=Coverage))+geom_point()

t$Filtering = gsub('filteval', 'Imputing', t$Filtering)
t$Filtering = gsub('phasedonly', 'Phasing', t$Filtering)

t %>% filter(Filtering != 'unImputing') %>% 
  ggplot(., aes(x=as.numeric(Coverage), y=value, color=name))+
  facet_grid(cols = vars(ind), rows=vars(Filtering))+
  geom_point()+
  geom_line()+
  labs(x='Coverage', y='Proportion', colour = 'Precision and Recall')+
  theme_classic()+
  scale_shape_discrete(labels = c('Filtered Imputed', 'Phased'))




gatk_concordance %>% ggplot(aes(x=Precision, y=Recall))+geom_line()+facet_grid(cols=vars(ind), rows = vars(Coverage))
cnam = c('sample', 'V1', 'V2', 'V3', 'V4', 'MISSING_ENTRY_truth-MISSING_ENTRY_test','MISSING_ENTRY_truth-MISSING_GT_test','MISSING_ENTRY_truth-REF_test','MISSING_ENTRY_truth-ALT_1_test','MISSING_ENTRY_truth-ALT_2_test', 'MISSING_GT_truth-MISSING_ENTRY_test',	'MISSING_GT_truth-MISSING_GT_test', 'MISSING_GT_truth-REF_test','MISSING_GT_truth-ALT_1_test', 'MISSING_GT_truth-ALT_2_test',	'REF_test-MISSING_ENTRY_truth',	'REF_test-MISSING_GT_truth','REF_truth-REF_test','REF_truth-ALT_1_test','REF_truth-ALT_2_test', 'ALT_1_truth-MISSING_ENTRY_test',	'ALT_1_truth-MISSING_GT_test',	'ALT_1_truth-REF_test',	'ALT_1_truth-ALT_1_test','ALT_1_truth-ALT_2_test','ALT_2_truth-MISSING_ENTRY_test',	'ALT_2_truth-MISSING_GT_test',	'ALT_2_truth-REF_test',	'ALT_2_truth-ALT_1_test',	'ALT_2_truth-ALT_2_test','ERROR')
setwd('/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/data/haplotype_concordance/')
concfilelist = list.files('.', pattern = 'by_sample.txt')
#read our shit
concdata =do.call(rbind, lapply(seq_along(concfilelist), function(i){fread(concfilelist[[i]], col.names = cnam) %>% mutate(fn=concfilelist[[i]])}))
#make columns for coverage and test file out of the filenames (hacky regex)
concdata$ind = sapply(str_split(concdata$fn, pattern = '_'),function(x) x[2])
fnstr = sapply(str_split(concdata$fn, pattern = 'AgamP4-'),function(x) x[3])
fnstr = gsub(pattern = '3L_', '', fnstr)
fnstr = gsub(pattern = '.by_sample.txt', '', fnstr)
concdata$level = fnstr
concdata$test_cov = str_match(gsub('^.+?Agam(.*)', "\\1",concdata$fn), "mform_\\s*(.*?)\\s*_AgamP4")[,2]
concdata$`MISSING_GT_truth-ALT_2_test`
concdata$`ALT_1_truth-ALT_2_test`

