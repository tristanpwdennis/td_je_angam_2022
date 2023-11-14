import subprocess
import itertools as it
import numpy as np
import os

#coluzzi
r = np.arange(0,158,1)

basedir = "/move_bioinfo/anopheles_03_21/angsd_output/by_species_wg_noinversions/"


freqfile = basedir + "m_form.freq"
glffile = basedir + "m_form.fullset.pop.wg_noinversions.beagle.gz"
outdir= basedir + 'm_form_outres'
ninds = 159
tnum=20


for a,b in it.permutations(r,2):
    outfile = outdir + str(a) + '_' + str(b) + '.res'
    subprocess.run(["~/packages/NgsRelate/ngsRelate", "-f", freqfile, "-G", glffile, "-n",str(ninds), "-a", str(a), "-b", str(b), "-O", outfile, "-p", str(tnum)], shell=True)

#gambiae
r = np.arange(0,155,1)

basedir = "/move_bioinfo/anopheles_03_21/angsd_output/by_species_wg_noinversions/"


freqfile = basedir + "s_form.freq"
glffile = basedir + "s_form.fullset.pop.wg_noinversions.beagle.gz"
outdir= basedir + 's_form_outres'
ninds = 155
tnum=20


for a,b in it.permutations(r,2):
    outfile = outdir + str(a) + '_' + str(b) + '.res'
    subprocess.run(["~/packages/NgsRelate/ngsRelate", "-f", freqfile, "-G", glffile, "-n",str(ninds), "-a", str(a), "-b", str(b), "-O", outfile, "-p", str(tnum)], shell=True)


