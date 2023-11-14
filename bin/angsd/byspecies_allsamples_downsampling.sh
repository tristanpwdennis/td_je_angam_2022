#PBS -N bychr_allsamples_downsampled
#PBS -S /bin/bash
#PBS -m abe
#PBS -M tristan.dennis@glasgow.ac.uk
#PBS -l nodes=1:ppn=20:centos7,walltime=20:00:00,cput=200:00:00


#allsamples together
for CHROM in AgamP4_2L AgamP4_2R AgamP4_3L AgamP4_3R AgamP4_X
	do
	for COV in 0.5 1 2.5 5 7.5 12 full
		do
		/export/home2/tdd3v/angsd/angsd -b /export/projects/III-data/$COV/tdd3v/anopheles/bam/lists/$COVbamlist_new.txt -gl 1 -domajorminor 4 -setMaxDepth 6000 -doCounts 1 -trim 0 -doPost 1 -doDepth 1 -doQsDist -C 50 -dumpCounts 2 -domaf 1 -nThreads 16 -doGlf 2 -uniqueonly 1 -minmaf 0.05 -doHWE 1 -doCheck 0 -minQ 30 -dosnpstat 1 -r $CHROM -snp_pval 0.05 -out /export/projects/III-data/tdd3v/anopheles/scripts/$SPEC.$CHROM_glf2_minorfilts_angsd_samtools_downsampled_core_$COV -ref /export/projects/III-data/tdd3v/anopheles/ref/vb_agamp3.fna

		done
	done


#allsamples together
for SPEC in m s
	do
	for CHROM in AgamP4_2L AgamP4_2R AgamP4_3L AgamP4_3R AgamP4_X
		do
		for COV in 0.5 1 2.5 5 7.5 12 full
			do
			/export/home2/tdd3v/angsd/angsd -b /export/projects/III-data/$COV/tdd3v/anopheles/bam/lists/$SPEC.$COVbamlist_new.txt -gl 1 -domajorminor 4 -setMaxDepth 6000 -doCounts 1 -trim 0 -doPost 1 -doDepth 1 -doQsDist -C 50 -dumpCounts 2 -domaf 1 -nThreads 16 -doGlf 2 -uniqueonly 1 -minmaf 0.05 -doHWE 1 -doCheck 0 -mininds 40 -minQ 30 -dosnpstat 1 -r $CHROM -snp_pval 0.05 -out /export/projects/III-data/tdd3v/anopheles/scripts/$SPEC.$CHROM_glf2_minorfilts_angsd_samtools_downsampled_core_$COV -ref /export/projects/III-data/tdd3v/anopheles/ref/vb_agamp3.fna

			done
		done
	done
