#PBS -N angsd_wg_ecozone_nodpfilt
#PBS -S /bin/bash
#PBS -m abe
#PBS -M tristan.dennis@glasgow.ac.uk
#PBS -l nodes=1:ppn=20:centos7,walltime=60:00:00,cput=900:00:00

source ~/.bashrc
conda activate angsd

#per form for wg, incorporating monomorphic sites for dxy and thetas calculation
for LIST in `ls /export/projects/III-data/lamberton/tdd3v/anopheles/by_ecotype/m*pop | sed s/.pop//g`
do
                MININD=$(expr `cat $LIST  | wc -l` / 4)
                REF=/export/projects/III-data/lamberton/tdd3v/anopheles/ref/vb_agamp3.fna
                TNUM=20
                ~/angsd/angsd  -b $LIST'.pop' -gl 1 -anc $REF -domajorminor 1  -doCounts 1 -trim 0 -doPost 1 -domaf 1 -nThreads $TNUM -doSaf 1 -uniqueonly 1 -minQ 30 -minInd $MININD -out $LIST'_nomaf_nodp_allsites' -ref $REF

done


for LIST in `ls /export/projects/III-data/lamberton/tdd3v/anopheles/by_ecotype/s*pop | sed s/.pop//g`
do
                MININD=$(expr `cat $LIST  | wc -l` / 4)
                REF=/export/projects/III-data/lamberton/tdd3v/anopheles/ref/vb_agamp3.fna
                TNUM=20
                ~/angsd/angsd  -b $LIST'.pop' -gl 1 -anc $REF -domajorminor 1  -doCounts 1 -trim 0 -doPost 1 -domaf 1 -nThreads $TNUM -doSaf 1 -uniqueonly 1 -minQ 30 -minInd $MININD  -out $LIST'_nomaf_nodp_allsites' -ref $REF

done

