library(tidyverse)
#define groups of samples for per-pop analysis

metadata_basedir = '/Users/tristanpwdennis/Projects/MOVE/td_je_angam_2022/metadata/'

prepare_metadata <- function(basedir) {
  #specify our files
  metadata_path = paste0(metadata_basedir, 'dw_complete_metadata_clusters_corrected.csv')
  bam_name_to_key_path = paste0(metadata_basedir, 'bam_name_to_seqid_key.txt')
  bam_name_to_id_path = paste0(metadata_basedir, 'bamorder_seqid.csv')
  #read our shit
  bamorder = read.csv(bam_name_to_id_path, header=T)
  metadata = read.csv(metadata_path)
  key = read.delim(bam_name_to_key_path)
  #compile metadata for sequenced samples
  metadata = left_join(metadata, key, by = c('ID' = 'seq_key'))
  sequenced_keys = left_join(bamorder, key, by =c('seq_id' = 'bam_key'))
  sequenced_samples = left_join(sequenced_keys, metadata, by = c('seq_id' = 'bam_key'))
  #list of outgroupsm, get rid of them
  outgroups = c('Aka22','Ash12','Dzo14','Mam9','NeA.18','NeA.23','Nsa1','Nsa16','Nsa18','Nsa4','OBP09','OBP18','OBP14','Tub8')
  sequenced_samples = sequenced_samples[ ! sequenced_samples$seq_id %in% outgroups, ]
  #add bam fullname
  sequenced_samples = sequenced_samples %>% mutate(bamid = paste0('/export/projects/III-data/lamberton/tdd3v/anopheles/bam/', seq_id,'.srt.dp.rg.bam' ))
  #prepare lat long
  sequenced_samples$Lat = as.numeric(gsub("N", "", sequenced_samples$Lat))
  sequenced_samples$long = as.numeric(gsub("W", "", sequenced_samples$long))
  return(sequenced_samples)
}

#function to import metadata with new keys for ordered bamlists on the hpc
import_formpops <- function(sampdf) {
  sform = filter(sequenced_samples, Form == 'S') %>% mutate(population = 's_form') %>% mutate(bam_key = seq(0, nrow(.)-1))
  mform = filter(sequenced_samples, Form == 'M') %>% mutate(population = 'm_form') %>% mutate(bam_key = seq(0, nrow(.)-1))
}

#import_ecopops <- function(sampdf) {
#  metadata = prepare_metadata(metadata_basedir)
#  habitat = unique(metadata$habitat)
#  species = unique(metadata$species)
#}

#to add functions to import ecopops and write all pops to csvs as an if statement to their function calls












