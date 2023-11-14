#PBS -N by_species_forngsrelate_noinversions
#PBS -S /bin/bash
#PBS -m abe
#PBS -M tristan.dennis@glasgow.ac.uk
#PBS -l nodes=1:ppn=20:centos7,walltime=20:00:00,cput=200:00:00

#gambiae
/export/home2/tdd3v/angsd/angsd -b /export/projects/III-data/tdd3v/anopheles/per_pop_downsampled/poplists/s_form.fullset.pop -gl 1 -anc /export/projects/III-data/tdd3v/anopheles/ref/vb_agamp3.fna -domajorminor 1 -trim 0 -doPost 1 -C 50 -domaf 1 -nThreads 20 -doGlf 3 -doSaf 1 -uniqueonly 1 -minmaf 0.05 -doHWE 1 -doCheck 0 -minQ 30 -minInd 38 -rf /export/home2/tdd3v/scripts/wg_no_inversions.txt -dogeno 1 -snp_pval 0.05 -out /export/projects/III-data/tdd3v/anopheles/per_pop_downsampled/poplists/s_form.fullset.pop.wg_noinversions -ref /export/projects/III-data/tdd3v/anopheles/ref/vb_agamp3.fna

#coluzzi
/export/home2/tdd3v/angsd/angsd -b /export/projects/III-data/tdd3v/anopheles/per_pop_downsampled/poplists/m_form.fullset.pop -gl 1 -anc /export/projects/III-data/tdd3v/anopheles/ref/vb_agamp3.fna -domajorminor 1 -trim 0 -doPost 1 -C 50 -domaf 1 -nThreads 20 -doGlf 3 -doSaf 1 -uniqueonly 1 -minmaf 0.05 -doHWE 1 -doCheck 0 -minQ 30 -minInd 40 -rf /export/home2/tdd3v/scripts/wg_no_inversions.txt -dogeno 1 -snp_pval 0.05 -out /export/projects/III-data/tdd3v/anopheles/per_pop_downsampled/poplists/s_form.fullset.pop.wg_noinversions -ref /export/projects/III-data/tdd3v/anopheles/ref/vb_agamp3.fna

