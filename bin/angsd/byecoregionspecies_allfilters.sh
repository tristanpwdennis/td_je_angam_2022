#PBS -N angsd_wg_ecoregion
#PBS -S /bin/bash
#PBS -m abe
#PBS -M tristan.dennis@glasgow.ac.uk
#PBS -l nodes=1:ppn=20:centos7,walltime=50:00:00,cput=600:00:00

source ~/.bashrc
export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH

for LIST in `ls /export/projects/III-data/tdd3v/anopheles/by_ecotype/m*pop | sed s/.pop//g`
do
		MININD=$(expr `cat $LIST  | wc -l` / 4)
		REF=/export/projects/III-data/lamberton/tdd3v/anopheles/ref/vb_agamp3.fna
		TNUM=20
		for CHROM in AgamP4_2L AgamP4_2R AgamP4_3L AgamP4_3R AgamP4_3R AgamP4_X
			do
				~/angsd/angsd  -b $LIST'.formpop' -gl 1 -anc $REF -domajorminor 4 -setMaxDepth 6000  -doCounts 1 -trim 0 -doPost 1 -domaf 1 -nThreads $TNUM -doSaf 1 -uniqueonly 1 -doHWE 1 -doCheck 0 -minmaf 0.05 -minQ 30 -minInd $MININD  -r $CHROM -snp_pval 0.05 -out $$LIST.$CHROM-ref $REF
			done
done

for LIST in `ls /export/projects/III-data/tdd3v/anopheles/by_ecotype/s*pop | sed s/.pop//g`
do
		MININD=$(expr `cat $LIST  | wc -l` / 4)
		REF=/export/projects/III-data/lamberton/tdd3v/anopheles/ref/vb_agamp3.fna
		TNUM=20
		for CHROM in AgamP4_2L AgamP4_2R AgamP4_3L AgamP4_3R AgamP4_3R AgamP4_X
			do
				~/angsd/angsd  -b $LIST'.formpop' -gl 1 -anc $REF -domajorminor 4 -setMaxDepth 6000  -doCounts 1 -trim 0 -doPost 1 -domaf 1 -nThreads $TNUM -doSaf 1 -uniqueonly 1 -doHWE 1 -doCheck 0 -minmaf 0.05 -minQ 30 -minInd $MININD  -r $CHROM -snp_pval 0.05 -out $LIST.$CHROM -ref $REF
			done
done
