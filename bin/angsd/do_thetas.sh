for f in `ls *.saf.idx | sed s/_nomaf_nodp_allsites.saf.idx//g`
do
~/packages/angsd/misc/realSFS $f"_nomaf_nodp_allsites.saf.idx" -fold 1 -cores 10 > $f"_nomaf_nodp_allsites.sfs"
~/packages/angsd/misc/realSFS saf2theta $f"_nomaf_nodp_allsites.saf.idx" -outname $f"_nomaf_nodp_allsites" -sfs $f"_nomaf_nodp_allsites.sfs" -fold 1
~/packages/angsd/misc/thetaStat do_stat $f"_nomaf_nodp_allsites.thetas.idx"
~/packages/angsd/misc/thetaStat do_stat $f"_nomaf_nodp_allsites.thetas.idx" -win 10000 -step 5000  -outnames $f"_nomaffilter.thetasWindow.gz"
done
