#!/usr/bin/env python
import glob
import os

infilepath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/FreeLigands/'
outfilepath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/FreeLigands/Renamed/'
renamedfilepath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/FreeLigands/Renamed/'

atomfile = open('renamingscript.sh', 'w')
print(atomfile)
atomfile.write('#!/bin/sh\n')
ligandlist = glob.glob('{}*.mol2'.format(infilepath))
for ligandfile in ligandlist:
    ligandname = ligandfile.replace(infilepath, '').replace('.mol2', '')
    InputFileName = '{}{}.mol2'.format(infilepath, ligandname)
    OutputFileName = '{}{}_rn.ac'.format(outfilepath, ligandname)
    # subprocess.call(['home/mdmannin/amber16/bin/atomtype', '-i ' + red_filename, '-f mol2',  '-o ' + ante_filename, '-p gaff'])
    atomfile.write('atomtype -i ' + InputFileName + ' -f mol2 -o ' + OutputFileName + ' -p gaff2 \n')
atomfile.close()
os.chmod('renamingscript.sh', 0744)
# rename = subprocess.Popen(['./renamingscript.sh'])
# rename.wait()
