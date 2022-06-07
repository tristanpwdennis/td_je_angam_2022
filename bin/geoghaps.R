#ls pkg
pkg = c("ape", "BiocManager", 'phangorn', 'phytools', 'pegas' ,'vcfR')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

metadata = read.csv('~/Projects/MOVE/td_je_angam_2022/metadata/metadata_with_insecticide_info.csv')
meta_h1 = metadata
meta_h2 = metadata
meta_h1$hapid = 'h1'
meta_h2$hapid = 'h2'
hapmeta = rbind(meta_h1, meta_h2)
hapmeta$indid = paste0('ind', "", hapmeta$bamorder, '_', hapmeta$hapid)
treedata = dplyr::select(hapmeta, indid, hapid, habitat, Form, Site)
treedata$indid




phy = ape::read.tree('~/Projects/MOVE/td_je_angam_2022/data/cyp9k_hap_alignment/cyp9k_allhaps.fas.treefile')
phy_p = midpoint.root(phy)

#create ordered treedata
treedata = treedata[match(phy$tip.label, treedata$indid), ] 

phy_p$species = treedata$Form
phy_p$site  =treedata$Site
phy_p$habitat = treedata$habitat
phy_p$tip.colors = c('orangered3', 'dodgerblue3')
phy_p$edge.length[phy_p$edge.length >  0.005] = 0.005
phy_p$edge.length[phy_p$edge.length == 0]    = 5e-5

phy_p$tip.label
as.factor(phy_p$species)

spec <- factor(phy_p$habitat)
spec
(mycol <- c("orangered3", "dodgerblue3", "turquoise2", "deeppink4")[spec])
(myshape <- c())
plot.phylo(phy_p, type = 'unr',
           use.edge.length = T, 
           show.tip.label = F,
           lab4ut = "axial",
           edge.color = "slategray3", font=0.5, edge.width = 0.5, node.depth = 1)
tiplabels(pch=16,col=mycol)



myfas = read.FASTA('~/Projects/MOVE/td_je_angam_2022/data/cyp9k_hap_alignment/from_samtools/af_001.grouped.rn.fasta')
names = gsub(pattern = "af_001/" , '',names(myfas))
names = gsub(pattern = "geneonly" , '',names)
names(myfas) = names

cyphaps = haplotype(myfas)
net = haploNet(cyphaps)
#pn = plot(net, attr(net, "freq"), fast = FALSE, show.mutation = 2, scale.ratio = 20, cex = 0.8, labels=FALSE, threshold = 0)

#make network metadata for pies
treedata$hapnetnames = paste0(hapmeta$seq_id, '_',hapmeta$hapid)
networkdata = treedata[match(names(myfas), treedata$hapnetnames), ]
#change fasta names to form id
names(myfas) <- networkdata$Form
#make table of now many individuals from each form belong to each cluster
ind.hap<-with(
     stack(setNames(attr(cyphaps, "index"), rownames(cyphaps))),
     table(hap=ind, individuals=names(myfas)[values]))

?plot.haploNet

pegas::diffHaplo(cyphaps)

net = haploNet(haplotype(myfas))


plot(net, attr(net, "freq"), fast = TRUE, legend = c(-700, 30), show.mutation = 2, scale.ratio = 20, cex = 0.8, labels=FALSE, threshold = 0,  pie=ind.hap, col = c('orangered3', 'dodgerblue3'))

labels(cyphaps)

vcf = vcfR::read.vcfR('~/Projects/MOVE/td_je_angam_2022/data/cyp9k_hap_alignment/from_freebayes/cyp9k_geneonly_q30_snps_mafall.vcf.gz')
maf(vcf)



m = metadata[metadata$Form == 'M',]$seq_id
s = metadata[metadata$Form == 'S',]$seq_id

mfreq = as.data.frame(maf(vcf[samples=m]))
sfreq = as.data.frame(maf(vcf[samples=s]))
form_freq_table = cbind(rownames(mfreq), mfreq$Frequency, sfreq$Frequency)
x = fread('~/Projects/MOVE/td_je_angam_2022/data/cyp9k_hap_alignment/from_freebayes/cyp9k_geneonly_q30_snps_mafall.snpeff', fill=T)
maftab = cbind(form_freq_table,x[,c(1,3,11,12)])
colnames(maftab) = c('allele_id','Maf_Ancol','Maf_Angam', 'Pos', 'Variant_Effect', 'NtChange', 'AAChange')
maftab$Maf_Ancol <- as.numeric(maftab$Maf_Ancol)
maftab$Maf_Angam <- as.numeric(maftab$Maf_Angam)

mafmat = as.matrix(maftab[,2:3], rownames = maftab$Variant_Effect)
heatmap(mafmat)

