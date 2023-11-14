import glob
import itertools as it
import os


for chrom in ['AgamP4_3L', 'AgamP4_3R', 'AgamP4_2L', 'AgamP4_2R', 'AgamP4_X']:
    ecos = ['savannah', 'rainforest', 'decidforest','mangrove']
    for speca, specb in it.combinations_with_replacement(['mform','sform'], 2):
        for a, b in it.combinations(ecos, 2):
            prefa = f'{speca}-{a}{chrom}'
            prefb = f'{specb}-{b}{chrom}'

            safa = f'{prefa}.saf.idx'
            safb = f'{prefb}.saf.idx'

            print(safa, safb)
            outsfs = chrom + '_' + prefa + '_' + prefb + '.2dsfs'
            fstfile = chrom + '_' + prefa + '_' + prefb
           	    	 #call realsfs
            os.system('~/packages/angsd/misc/realSFS ' + safa + ' ' + safb + ' -fold 1 -P 20 '+ ' > fst/' + outsfs)
            os.system('~/packages/angsd/misc/realSFS ' + 'fst index -whichFst 1 ' + safa + ' ' + safb  + ' -sfs fst/' + outsfs + ' -fstout fst/' + fstfile + ' -P 20')
            os.system('~/packages/angsd/misc/realSFS ' + 'fst stats2 fst/' + fstfile + '.fst.idx ' + '-win 10000 -step 5000 > fst/' + fstfile + '.window10step0p2.txt' + ' -P 20')
            os.system('~/packages/angsd/misc/realSFS ' + 'fst stats fst/' + fstfile + '.fst.idx ' + '> fst/' + fstfile + '.globalstat.txt')

            #do dxy
            #os.system('Rscript do_dxy.R -p ' + safa + ' -q ' + safb + ' -t 252604662 -o ' + fstfile)
