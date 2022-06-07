library(tidyverse)
#define groups of samples for per-pop analysis

metadata = read.csv('~/Projects/MOVE/data/anopheles_03_21_DW/metadata/dw_complete_metadata_clusters_corrected.csv')
key = read.delim('~/Projects/MOVE/data/anopheles_03_21_DW/metadata/bam_name_to_seqid_key.txt')
metadata = left_join(metadata, key, by = c('ID' = 'seq_key'))
bamorder = read.csv('~/Projects/MOVE/data/anopheles_03_21_DW/metadata/bamorder_seqid.csv', header=T)
sequenced_keys = left_join(bamorder, key, by =c('seq_id' = 'bam_key'))
sequenced_samples = left_join(sequenced_keys, metadata, by = c('seq_id' = 'bam_key'))
nam = sequenced_samples %>% select(seq_id)

outgroups = c('Aka22','Ash12','Dzo14','Mam9','NeA.18','NeA.23','Nsa1','Nsa16','Nsa18','Nsa4','OBP09','OBP18','OBP14','Tub8')

sequenced_samples = sequenced_samples[ ! sequenced_samples$seq_id %in% outgroups, ]

sequenced_samples = sequenced_samples %>% mutate(bamid = paste0('/export/projects/III-data/lamberton/tdd3v/anopheles/bam/', seq_id,'.srt.dp.rg.bam' ))


sequenced_samples$Lat = as.numeric(gsub("N", "", sequenced_samples$Lat))
sequenced_samples$long = as.numeric(gsub("W", "", sequenced_samples$long))

sequenced_samples

sform = filter(sequenced_samples, Form == 'S') %>% mutate(population = 's_form') %>% mutate(bam_key = seq(0, nrow(.)-1))
mform = filter(sequenced_samples, Form == 'M') %>% mutate(population = 'm_form') %>% mutate(bam_key = seq(0, nrow(.)-1))

Sites = sform %>% select(Site) %>% group_by(Site) %>% count() %>% filter(n > 9)

savannah = sequenced_samples %>% filter(habitat == 'coastal_savannah') %>% mutate(bam_key = seq(1, nrow(.))) %>% mutate(population = 'savannah')
decid_forest = sequenced_samples %>% filter(habitat == 'deciduous_forest') %>% mutate(bam_key = seq(1, nrow(.))) %>% mutate(population = 'decid_forest')
mangrove = sequenced_samples %>% filter(habitat == 'mangrove_swamp') %>% mutate(bam_key = seq(1, nrow(.))) %>% mutate(population = 'mangrove')
rainforest = sequenced_samples %>% filter(habitat == 'rainforest') %>% mutate(bam_key = seq(1, nrow(.))) %>% mutate(population = 'rainforest')

savannah %>% filter(Form == 'S')

mangrove = sequenced_samples %>% filter(habitat == 'mangrove_swamp')

metadata %>% filter(habitat == 'mangrove_swamp') %>% group_by(Form) %>% count()


#for (site in sequenced_samples$Site %>% unique()) {
#  t = sequenced_samples %>% filter(Site == site) %>% select(seq_id)
#  write_csv(t, paste0('/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/',site,".pop"))
#}

mangrove %>% group_by(Form) %>% count()

x = sequenced_samples %>% group_by(Site) %>% count()
x %>% filter(n > 4)
#
#sform %>% select(seq_id) %>% write_csv(., col_names = F, file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/s_form.pop')
#mform %>% select(seq_id) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/m_form.pop')
#mangrove %>% select(seq_id) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/mangrove.pop')
#decid_forest %>% select(seq_id) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/decid_forest.pop')
#savannah %>% select(seq_id) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/savannah.pop')
#rainforest %>% select(seq_id) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/rainforest.pop')

mangrove %>% filter(Form == 'M') %>% select(bamid) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/m_mangrove.pop')
mangrove %>% filter(Form == 'S') %>% select(bamid) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/s_mangrove.pop')

decid_forest %>% filter(Form == 'M') %>% select(bamid) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/m_decid_forest.pop')
decid_forest %>% filter(Form == 'S') %>% select(bamid) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/s_decid_forest.pop')

savannah  %>%   filter(Form == 'M') %>% select(bamid) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/m_savannah.pop')
savannah  %>%   filter(Form == 'S') %>% select(bamid) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/s_savannah.pop')

rainforest %>% filter(Form == 'M') %>% select(bamid) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/m_rainforest.pop')
rainforest %>% filter(Form == 'S') %>% select(bamid) %>% write_csv(., col_names = F,file = '/Users/tristanpwdennis/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/by_pop_downsampled/s_rainforest.pop')

