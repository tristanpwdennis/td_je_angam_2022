#----------------------
#Haplotype analysis
#Tristan Dennis May 2022
#----------------------

#----------------------
#HSet up environment and clean crap metadata
#----------------------

#ls pkg
pkg = c("network", "GGally", 'phangorn', 'phytools', 'pegas' ,'vcfR', 'RColorBrewer', 'data.table', 'tidyverse', "rnaturalearth", "rnaturalearthdata", "sf", 'rehh')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

fread('')
metadata=fread('~/OneDrive - University of Glasgow/MOVE/manuscripts/td_je_angam_2022/metadata/metadata_full_sequencedonly_062022.csv')

metadata %>% filter(Form == 'S' & habitat == 'deciduous_forest' & INSECTICIDE_USE == 'high')
#create haplotype specific metadata (duplicate ind rows and add hap1/hap2 designation)
meta_h1 = metadata
meta_h2 = metadata
meta_h1$hapid = 'h1'
meta_h2$hapid = 'h2'
hapmeta = rbind(meta_h1, meta_h2)
hapmeta$indid = paste0(hapmeta$seq_id, '_', hapmeta$hapid)


#load geog data
world <- ne_countries(scale = 'large', returnclass = "sf", continent = 'africa')
lakes <- ne_download(scale = 'large', type = 'lakes', category = 'physical', returnclass = "sf")
rivers <- ne_download(scale = 'large', type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
ocean <- ne_download(scale = 'large', type = 'ocean', category = 'physical', returnclass = "sf")

#load alignment and create hap membership tables

#load alignment, clean, add clean names
cypfas = read.FASTA('~/Projects/MOVE/td_je_angam_2022/data/cyp9k_hap_alignment/from_samtools/af_001.grouped.rn.fasta')
acefas = read.FASTA('~/Projects/MOVE/td_je_angam_2022/data/ace1_alignment/ace1_af_001.grouped.fasta')

names = gsub(pattern = "af_001/" , '',names(cypfas))
names = gsub(pattern = "geneonly" , '',names)
names(cypfas) = names
cypfas

#load phylogenies and make tip data
cyptreefile = '~/Projects/MOVE/td_je_angam_2022/data/cyp9k_hap_alignment/from_samtools/af_001.grouped.rn.fasta.treefile'
acetreefile = '~/Projects/MOVE/td_je_angam_2022/data/ace1_alignment/ace1_af_001.grouped.fasta.treefile'
treedata = dplyr::select(hapmeta, indid, hapid, habitat, Form, Site) #make tree data

#----------------------
#Plot CYP9K1 minimum spanning network
#right now the network layout is broken as fuck and I will only publish the phylogeny, but it's still useful to see
#----------------------
#get haplotypes from alignment (clustering sequences basically) and make network
cyphaps = haplotype(cypfas)
net = haploNet(cyphaps)
#

#make df of individuals and respective haplotype group memberships
allind.hap<-with(
  stack(setNames(attr(cyphaps, "index"), rownames(cyphaps))),
  table(hap=ind, individuals=names(cypfas)[values]))
hapmembership = data.frame(allind.hap)
hapmembership = hapmembership[hapmembership$Freq > 0,]
hapmeta = left_join(hapmeta, hapmembership, by=c('indid' = 'individuals'))

#make network metadata for pies
treedata$hapnetnames = paste0(hapmeta$seq_id, '_',hapmeta$hapid)
networkdata = treedata[match(names(cypfas), treedata$hapnetnames), ]

#change fasta names to form id
names(cypfas) <- networkdata$Form

#make table of now many individuals from each form belong to each cluster
ind.hap<-with(
  stack(setNames(attr(cyphaps, "index"), rownames(cyphaps))),
  table(hap=ind, individuals=names(cypfas)[values]))

x = plot(net, attr(net, "freq"), fast = T, show.mutation = 2, scale.ratio = 10, cex = 0.5, labels=FALSE, threshold = 0,  pie=ind.hap)

#----------------------
#Plot CYP9K1 and ACE1 phylogeny with and without outgroups
#----------------------

plot_phy <- function(treefile) {
  phy = ape::read.tree(treefile)
  phy$tip.label <- gsub(pattern = 'geneonly', replacement = '', phy$tip.label)
  #read tree, root
  phy_p = midpoint.root(phy)
  #create ordered treedata
  treedata = treedata[match(phy$tip.label, treedata$indid), ] 
  #configure tip labels and set lengths and colours
  phy_p$species = treedata$Form 
  phy_p$site  =treedata$Site
  phy_p$habitat = treedata$habitat
  phy_p$tip.colors = c('#D95F02', "#1B9E77")
  phy_p$edge.length[phy_p$edge.length >  0.005] = 0.005
  phy_p$edge.length[phy_p$edge.length == 0]    = 5e-5
  spec <- factor(phy_p$species)
  mycol <- c('#D95F02', "#1B9E77")[spec]
  spec
  #pliot tree and tiplabs
  par(mar=c(0.5,0.5,0.5,0.5))
  plot.phylo(phy_p, type = 'unr',
             use.edge.length = T, 
             show.tip.label = F,
             lab4ut = "axial",
             edge.color = "slategray3", font=0.5, edge.width = 0.5, node.depth = 1)
  tiplabels(pch=16,col=mycol)
  add.scale.bar(cex=0.7, "topleft")
}

#make final phylogeny plot object
par(mfcol = c(2, 1))
layout.matrix <- matrix(c(0, 1), nrow = 2, ncol = 1)
plot_phy(cyptreefile)
title("cyp9k1                                                           ")
plot_phy(acetreefile)
title("ace1                                                           ")
legend("bottomleft", legend = c("An. gambiae", "An. coluzzi"), pch = 22,
       pt.bg = c('#D95F02', "#1B9E77"), pt.cex = 2.5, bty="n", box.lty=1)

?par
#----------------------
#Plot CYP9K1 variant frequency table by form and site
#----------------------

#get 'm' seq ids
m = metadata[metadata$species == 'M',]$seq_id
s = metadata[metadata$species == 'S',]$seq_id

cyp_vcf = '~/Projects/MOVE/td_je_angam_2022/data/cyp9k_hap_alignment/from_samtools/cyp9k_samtools.minmaf0.01.b5phased.finfix.eff'
ace_vcf = '~/Projects/MOVE/td_je_angam_2022/data/ace1_alignment/ace1_samtools.minmaf0.01.b5phased.finfix.eff'

make_variant_table <- function(hapvcf){
  vcf = read.vcfR(paste0(ace_vcf, '.vcf.gz')) #load vcf
  #get formwise maf
  mfreq = as.data.frame(maf(vcf[samples=m]))
  sfreq = as.data.frame(maf(vcf[samples=s]))
  #make into table and read snpeff table
  form_freq_table = cbind(rownames(mfreq), mfreq$Frequency, sfreq$Frequency)
  x = fread(paste0(hapvcf, '.snpeff'), fill = TRUE)
  #bind snpeff table to freq table
  maftab = cbind(form_freq_table,x[,c(1,3,11,12)])
  maftab = cbind(maftab, vcf@fix[,4:5])
  colnames(maftab) = c('allele_id','Maf_Ancol','Maf_Angam', 'Pos', 'Variant_Effect', 'NtChange', 'AAChange', 'ALT', 'REF')
  #add colnames and clean again
  maftab$Maf_Ancol <- as.numeric(maftab$Maf_Ancol)
  maftab$Maf_Angam <- as.numeric(maftab$Maf_Angam)
  newmaftab = maftab %>% dplyr::filter(Variant_Effect == 'missense_variant' | Variant_Effect == '3_prime_UTR_variant' | Variant_Effect == '5_prime_UTR_variant')
  #extract 'functional' variants - missenses and UTR variants
  newmaftab$varid = paste0(newmaftab$Variant_Effect, '|', newmaftab$NtChange, '|', newmaftab$AAChange)
  newmaftab$shortid = sub("^\\D+", "", newmaftab$NtChange)

  return(newmaftab) #return final variant table
}


acetab = make_variant_table(ace_vcf)
acetab = acetab[Variant_Effect == 'missense_variant']
cyptab = make_variant_table(cyp_vcf)
cyptab = pivot_longer(cyptab, cols=2:3)
acetab = pivot_longer(acetab, cols=2:3)

#plot table ####to-do - FIX Y LABS
colcyt = ggplot(cyptab[cyptab$name == 'Maf_Ancol',], aes(x=name, y=varid, fill=value))+
  geom_tile(color = "grey99",lwd = 1,linetype = 1)+
  scale_fill_gradient(low = "#f7fcb9", high="#31a354")+
  geom_text(aes(label = round(value, digits = 2)), color = "#1a1a1a", size = 4) +
  theme_classic()+
  annotate("segment", x = 2.55, xend = 2.55, y = 0.5, yend = 7.5)+
  annotate("segment", x = 2.55, xend = 2.55, y = 7.6, yend = 8.5)+
  annotate("segment", x = 2.55, xend = 2.55, y = 8.6, yend = 10.5)+
  annotate(geom="text", x=3, y=3, label='3\' UTR ' )+
  annotate(geom="text", x=3, y=8, label='5\' UTR ' )+
  annotate(geom="text", x=3, y=9.5, label='Missense ' )+
  scale_y_discrete(labels=c("187T>G","163A>C","118G>A","91T>C","32T>A","28T>C","6C>A","973G>C","671A>T","45T>C"))+
  labs(x='', y='')+ 
  theme(legend.position="none",
  axis.text.x = element_text(color = "grey20", face = "plain"),
  axis.text.y = element_text(color = "grey20", angle = 20, hjust = 0.9, vjust = 0.8, face = "plain"),
  text = element_text(size = 15))+
  coord_cartesian(clip = 'off')+
  theme(plot.margin = unit(c(1,5,1,1), "lines"))


gamcyt = ggplot(cyptab[cyptab$name == 'Maf_Angam',], aes(x=name, y=varid, fill=value))+
  geom_tile(color = "grey99",lwd = 1,linetype = 1)+
  scale_fill_gradient(low = "#fff7bc", high="#d95f0e")+
  geom_text(aes(label = round(value, digits = 2)), color = "#1a1a1a", size = 4) +
  theme_classic()+
  annotate("segment", x = 2.55, xend = 2.55, y = 0.5, yend = 7.5)+
  annotate("segment", x = 2.55, xend = 2.55, y = 7.6, yend = 8.5)+
  annotate("segment", x = 2.55, xend = 2.55, y = 8.6, yend = 10.5)+
  annotate(geom="text", x=3, y=3, label='3\' UTR ' )+
  annotate(geom="text", x=3, y=8, label='5\' UTR ' )+
  annotate(geom="text", x=3, y=9.5, label='Missense ' )+
  scale_y_discrete(labels=c("187T>G","163A>C","118G>A","91T>C","32T>A","28T>C","6C>A","973G>C","671A>T","45T>C"))+
  labs(x='', y='')+ 
  theme(legend.position="none",
        axis.text.x = element_text(color = "grey20", face = "plain"),
        axis.text.y = element_text(color = "grey20", angle = 20, hjust = 0.9, vjust = 0.8, face = "plain"),
        text = element_text(size = 15))+
  coord_cartesian(clip = 'off')+
  theme(plot.margin = unit(c(1,5,1,1), "lines"))


cyptab$shortid

#plot table ####to-do - FIX Y LABS
act = ggplot(acetab, aes(x=name, y=varid, fill=value))+
  geom_tile(color = "grey99",lwd = 1,linetype = 1)+
  scale_fill_gradient(low = "#fff7bc", high="#d95f0e")+
  geom_text(aes(label = round(value, digits = 2)), color = "#1a1a1a", size = 4) +
  scale_y_discrete(labels=c("104T>C","193G>T","838G>A","1431G>T"))+
  theme_classic()+
  annotate("segment", x = 2.55, xend = 2.55, y = 0.5, yend = 4.5)+
  annotate(geom="text", x=3, y=2.5, label='Missense ' )+
  labs(x='', y='')+ 
  coord_cartesian(clip = 'off')+
  theme(legend.position="none",
        axis.text.x = element_text(color = "grey20", size = 10, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 20, hjust = 0.9, vjust = 0.8, face = "plain"),
      text = element_text(size = 15))+
  theme(plot.margin = unit(c(1,5,1,1), "lines"))

act

cowplot::plot_grid(act, cyt, rel_heights  = c(2.5, 1), align = 'h', nrow = 1)

#----------------------
#Get GT table, convert to nucletodide haplotypes
#----------------------


#another horrible function to make a cyp9 gt table (has gt lookip)
#make_cyp_gt_table <- function(vcf_file) {
#readvcf
vcf = read.vcfR(paste0(cyp_vcf, '.vcf.gz')) #load vcf
gt_table = as.data.frame(extract.gt(vcf, return.alleles = TRUE)) #make GT df)
gt_table$allel_id = rownames(gt_table) #use rownames to make chr/pos ids
gt_table = separate(gt_table, col = allel_id, into = c('chr', 'chrx', 'pos')) 
gt_table = cbind(gt_table, cypeff[,c(1,2,3)]) #bind to eff info
gt_table_sub = gt_table %>% dplyr::filter(V3 == 'missense_variant' | V3 == '3_prime_UTR_variant' | V3 == '5_prime_UTR_variant') #filter out garbage
gt_tb_long = gt_table_sub %>% pivot_longer(cols=1:314)
septab = separate_rows(gt_tb_long, value) #sep each hap into a row
hapseptab =  cbind(septab,rep(c('h1', 'h2'), nrow(septab))) #add hap ids
hapseptab$hapid = paste0(hapseptab$name, '_', hapseptab$`rep(c("h1", "h2"), nrow(septab))`)

gt_tb = left_join(hapseptab, hapmeta, by = c('hapid'= 'indid'))
haptab = gt_tb %>% select(name, pos, value, V3, Site, hap) %>% distinct()
#return(gt_table_nucs)
b = haptab %>% filter(hap == 'LXXV')#}
c = haptab %>% filter(hap == 'VIII')#}

m = haptab %>% select(hap, value, pos) %>% distinct() %>% pivot_wider(names_from = hap)

################################################

#Individual analysis of haplotype memberships, with haps
################################################
hapmembership = hapmembership[hapmembership$Freq > 0,]

head(hapmembership)
hapmembership = dplyr::left_join(hapmembership, hapmeta, by=c('individuals' = 'indid'))

dplyr::group_by(hapmembership, hap.x)
x = hapmembership %>% group_by(Form, Site, hap.x, Lat, Long) %>% count()

hapmembership



scatpwide = pivot_wider(x,id_cols = c('Form', 'Site', 'Lat', 'Long'), names_from = hap.x, values_from = n)
scatpwide[is.na(scatpwide)] = 0
scatpwide = scatpwide[scatpwide$Lat > 2,]
scatpwide$rad = apply(scatpwide[5:45], 1, sum)
scatpwide$rad
mscatpwide = scatpwide#[scatpwide$Form == 'S',]


#############
#draw map
#############
#which individuals contain missense mutation hap or transcription hap

adafoah <- data.frame(x1 = 0.7, y1 = 5.4, x2 = 0.66410, y2 = 5.80350, site = 'adafoah')
anloga <- data.frame(x1 = 1, y1 = 5.5, x2 = 0.86514, y2 = 5.78630, site = 'anloga')
ashaiman <- data.frame(x1 = 0.3, y1 = 5.4, x2 = -0.05237, y2 = 5.69530, site = 'ashaiman')
dzorwulu <-  data.frame(x1 = -0.05, y1 = 5.33, x2 = -0.20592, y2 = 5.60499, site = 'dzorwulu')
mamprobi <-  data.frame(x1 = -0.35, y1 = 5.2, x2 = -0.23770, y2 = 5.54100, site = 'mamprobi')

#create new long lat columns
mscatpwide$newlong <- mscatpwide$Long
mscatpwide$newlat <- mscatpwide$Lat
#manually adjust long
mscatpwide$newlong[mscatpwide$Site == "Ada Foah"] <- 0.7
mscatpwide$newlong[mscatpwide$Site == "Anloga"] <- 1
mscatpwide$newlong[mscatpwide$Site == "Ashaiman"] <- 0.3 
mscatpwide$newlong[mscatpwide$Site == "Dzorwulu"] <- -0.05
mscatpwide$newlong[mscatpwide$Site == "Mamprobi"]  <- -0.35

#manually adjust lat
mscatpwide$newlat[mscatpwide$Site == "Ada Foah"] <- 5.4
mscatpwide$newlat[mscatpwide$Site == "Anloga"] <- 5.5
mscatpwide$newlat[mscatpwide$Site == "Ashaiman"] <- 5.4 
mscatpwide$newlat[mscatpwide$Site == "Dzorwulu"] <- 5.33
mscatpwide$newlat[mscatpwide$Site == "Mamprobi"]  <- 5.2

pt = ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-3.6, 1.5), ylim = c(4, 8), expand = FALSE)+
 # geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = adafoah)+
  #geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = anloga)+
  #geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = ashaiman)+
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = dzorwulu)+
  #geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = mamprobi)+
  geom_scatterpie(aes(x=newlong, y=newlat, group=Site, r=0.12), data=mscatpwide, position_dodge(width=2),
                  cols=colnames(mscatpwide[5:132]))+theme_classic()+
  labs(x='Longitude', y='Latitude', fill = 'Haplotype')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = 'none')+
        facet_wrap(~Form)
pt
ggsave('~/Desktop/geohapplot.svg', pt, device='svg', units = 'px')

#----------------------
#Analysis of haplotype variation over geography, subset to specific variants
#----------------------

#snpeff tables
aceeff = fread('~/Projects/MOVE/td_je_angam_2022/data/ace1_alignment/ace1_samtools.b5phased.finfix.eff.snpeff', fill=TRUE)
cypeff = fread('~/Projects/MOVE/td_je_angam_2022/data/cyp9k_hap_alignment/from_samtools/cyp9k_samtools.minmaf0.01.b5phased.finfix.eff.snpeff', fill=TRUE)

#function  to loop over forms and sites, extract maf from vcf for each
get_sitewise_maf_table <- function(efftab, vcf) {
vcf = read.vcfR(paste0(vcf, '.vcf.gz')) #load vcf
freqlist = list()
for (form in c('M', 'S')) {
  for (site in unique(metadata$Site)) {
    t =  filter(metadata, species == form & Site == site)
    if(nrow(t) >0) {
      print("found samples yay!")
      freq = as.data.frame(maf(vcf[samples=t$seq_id]))
      fnam = paste0(form, '_', site, '_', 'Frequency')
      colnames(freq) <- c('nAllele', 'gub', 'Count', "Freq")
      freq$name = fnam
      freq$site = site
      freq$form = form
      freq$allele = rownames(freq)
      freq = cbind(freq, efftab[,c(1,3,11,12)])
      freqlist[[fnam]] <- freq
    }
    else{
      print("No samples for this one")
      
    }
  }
}
sitefreqtab = do.call("rbind", freqlist)
return(sitefreqtab)
}

cypfreqtab =get_sitewise_maf_table(cypeff, cyp_vcf)
acefreqtab =get_sitewise_maf_table(aceeff, ace_vcf)

 #subset variants. in the interest of plottability we want just the functional CYP9K variants, for ACE1 we want the missense variants from ACE1
newcyptab =  dplyr::filter(cypfreqtab, V3 == 'missense_variant' | V3 == '3_prime_UTR_variant' | V3 == '5_prime_UTR_variant' & Freq > 0.02)
newacetab =  dplyr::filter(cypfreqtab, V3 == 'missense_variant' & Freq > 0.05)

newcyptab$varid = paste0(newcyptab$V3, '|', newcyptab$V11, '|', newcyptab$V12)
newacetab$varid = paste0(newacetab$V3, '|', newacetab$V11, '|', newacetab$V12)



#tidy by site af data
#join to lat long
submaftab = newcyptab %>% dplyr::select(site, form, varid, Freq)
sitedat = metadata %>% dplyr::select(Site, Lat, Long) %>% unique()
submaftab = left_join(submaftab, sitedat, by =c('site' = 'Site'))
mafmembership = pivot_wider(submaftab, id_cols = c('form', 'site', 'Lat', 'Long'), names_from = varid, values_from = Freq)
mafmembership$rad = apply(mafmembership[5:14], 1, sum)


asnile = submaftab %>% filter(varid == 'missense_variant|c.671A>T|p.Asn224Ile')
asnile$`671T` = asnile$Freq
asnile$`671A` = 1 - asnile$Freq
asnmembership = pivot_wider(asnile, id_cols = c('form', 'site', 'Lat', 'Long'), names_from = varid, values_from = c('671T', '671A'))
ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-3.6, 1.5), ylim = c(4, 8), expand = FALSE)+
  geom_scatterpie(aes(x=Long, y=Lat, group=site, r=0.1), data=asnmembership, position_dodge(width=2),
                  cols=colnames(asnmembership[5:6]))+theme_classic()+facet_wrap(~form)


asnile = submaftab %>% filter(varid == 'missense_variant|c.973G>C|p.Val325Leu')
asnile$`973C` = asnile$Freq
asnile$`973G` = 1 - asnile$Freq
asnmembership = pivot_wider(asnile, id_cols = c('form', 'site', 'Lat', 'Long'), names_from = varid, values_from = c('973C', '973G'))
ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-3.6, 1.5), ylim = c(4, 8), expand = FALSE)+
  geom_scatterpie(aes(x=Long, y=Lat, group=site, r=0.1), data=asnmembership, position_dodge(width=2),
                  cols=colnames(asnmembership[5:6]))+theme_classic()+facet_wrap(~form)


asnile = submaftab %>% filter(varid == '3_prime_UTR_variant|c.*91T>C|')
asnile$`91T` = asnile$Freq
asnile$`91C` = 1 - asnile$Freq
asnmembership = pivot_wider(asnile, id_cols = c('form', 'site', 'Lat', 'Long'), names_from = varid, values_from = c('91T', '91C'))


ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-3.6, 1.5), ylim = c(4, 8), expand = FALSE)+
  geom_scatterpie(aes(x=Long, y=Lat, group=site, r=0.3), data=asnmembership, position_dodge(width=2),
                  cols=colnames(asnmembership[5:6]))+theme_classic()+facet_wrap(~form)

asnile = submaftab %>% filter(varid == '3_prime_UTR_variant|c.*28T>C|')
asnile$`28T` = asnile$Freq
asnile$`28C` = 1 - asnile$Freq
asnmembership = pivot_wider(asnile, id_cols = c('form', 'site', 'Lat', 'Long'), names_from = varid, values_from = c('28T', '28C'))
ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-3.6, 1.5), ylim = c(4, 8), expand = FALSE)+
  geom_scatterpie(aes(x=Long, y=Lat, group=site, r=0.1), data=asnmembership, position_dodge(width=2),
                  cols=colnames(asnmembership[5:6]))+theme_classic()+facet_wrap(~form)

#tidy by site af data
#join to lat long
submaftab = newcyptab %>% dplyr::select(site, form, varid, Freq)
sitedat = metadata %>% select(Site, Lat, Long) %>% unique()
submaftab = left_join(submaftab, sitedat, by =c('site' = 'Site'))
mafmembership = pivot_wider(submaftab, id_cols = c('form', 'site', 'Lat', 'Long'), names_from = varid, values_from = Freq)
mafmembership$rad = apply(mafmembership[5:14], 1, sum)


asnile = submaftab %>% filter(varid == 'missense_variant|c.671A>T|p.Asn224Ile')
asnile$`671T` = asnile$Freq
asnile$`671A` = 1 - asnile$Freq
asnmembership = pivot_wider(asnile, id_cols = c('form', 'site', 'Lat', 'Long'), names_from = varid, values_from = c('671T', '671A'))
ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-3.6, 1.5), ylim = c(4, 8), expand = FALSE)+
  geom_scatterpie(aes(x=Long, y=Lat, group=site, r=0.1), data=asnmembership, position_dodge(width=2),
                  cols=colnames(asnmembership[5:6]))+theme_classic()+facet_wrap(~form)

b = data.frame(ind.hap)


#----------------------
#AHapstats
#----------------------
m_hap_file = '~/Projects/MOVE/td_je_angam_2022/data/genotype_likelihood_data/phasing_b4/mform_X_b4phasedgls.vcf.gz'
s_hap_file = '~/Projects/MOVE/td_je_angam_2022/data/genotype_likelihood_data/phasing_b4/sform_X_b4phasedgls.vcf.gz'
m_hap_file = '~/Projects/MOVE/td_je_angam_2022/data/genotype_likelihood_data/phasing_b4/mform_X_b4phasedgls.vcf.gz'


mh <- data2haplohh(hap_file = m_hap_file,
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

sh <- data2haplohh(hap_file = s_hap_file,
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")




#m_hh_scan = scan_hh(mh, polarized = FALSE)
#s_hh_scan = scan_hh(sh, polarized = FALSE)
s_hh_scan=fread('~/Projects/MOVE/td_je_angam_2022/data/genome_scans/s_hh.txt')
m_hh_scan=fread('~/Projects/MOVE/td_je_angam_2022/data/genome_scans/m_hh.txt')



#which plaotrype is hap of interest (middle of cyp9k)
mpos = which(mh@positions == 15241588, arr.ind=TRUE)

mres <- calc_ehh(mh, 
                mrk = mpos, 
                include_nhaplo = TRUE, polarized = FALSE)
plot(mres)
spos = which(sh@positions == 15241461, arr.ind=TRUE)
sres <- calc_ehh(sh, 
                 mrk = spos, 
                 include_nhaplo = TRUE, polarized = FALSE)
plot(sres)

pal=c('#4daf4a', '#984ea3')
s= pivot_longer(sres$ehh, cols=2:3) %>% ggplot(aes(x=POSITION, y=value, color=name))+
  geom_line(size=1)+
  theme_classic()+
  scale_colour_manual(values =c('#cc4c02', '#fec44f'), labels=c('Major', 'Minor'))+
  labs(x='Position', y='EHH')
s  
me = pivot_longer(mres$ehh, cols=2:3) %>% ggplot(aes(x=POSITION, y=value, color=name))+
  geom_line()+
  geom_line(size=1)+
  theme_classic()+
  scale_colour_manual(values =c('#238443', '#addd8e'), labels=c('Major', 'Minor'))+
  theme(legend.title = element_text('Allele'))+
  labs(x='Position', y='EHH')
me
cowplot::plot_grid(me, s, align='v', ncol=1)
#
