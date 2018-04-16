#!/usr/bin/env python
import os
import itertools
# #execfile('/home/mdmannin/amber16/amber.sh')
# #subprocess.call(['source', 'home/mdmannin/amber16/amber.sh'])
batchfilename = 'batch_input_disc.txt'
batchfile = open(batchfilename, 'r')
filepath = '/home/mdmannin/Desktop/Nanoparticles/Nanodisc/'
thicklist, radiuslist, curvelist, densitylist, ligandlist,ratiolist = [], [], [], [],[], []
for line in batchfile:
    line = line.strip()
    if line.find("THICKNESS=") != -1:
        thickbite = line[line.index('=')+1:].split(',')
        for i in thickbite: thicklist.append(int(i))
    elif line.find("RADIUS=") != -1:
        radiusbite = line[line.index('=')+1:].split(',')
        for i in radiusbite:radiuslist.append(int(i))
    elif line.find("CURVE=") != -1:
        curvebite = line[line.index('=')+1:].split(',')
        for i in curvebite:curvelist.append(int(i))
    elif line.find("DENSITY=") != -1:
        densitybite = line[line.index('=')+1:].split(',')
        for i in densitybite:densitylist.append(int(i))
    elif line.find("LIGAND1=") != -1:
        ligandbite = line[line.index('=')+1:].split(',')
        for i in ligandbite:ligandlist.append(i)
    elif line.find("RATIO=") != -1:
        ratiobite = line[line.index('=')+1:].split(',')
        for i in ratiobite:ratiolist.append(int(100*float(i)))
atomfile = open('discrenamingscript.sh', 'w')
atomfile.write('#!/bin/sh\n')
for thick, radius, curve, density, ligand, ratio in itertools.product(thicklist, radiuslist, curvelist, densitylist, ligandlist, ratiolist):
    goldname = ('{}x{}-{}curve-{}ratio{}dense{}'.format(thick, radius, curve, ratio, density, ligand))
    InputFileName = '{}discmol2/{}.mol2'.format(filepath,goldname)
    OutputFileName = '{}Atomtype/{}_rn.ac'.format(filepath,goldname)
    #subprocess.call(['home/mdmannin/amber16/bin/atomtype', '-i ' + InputFileName, '-f mol2',  '-o ' + OutputFileName, '-p gaff'])
    atomfile.write('atomtype -i ' + InputFileName + ' -f mol2 -o ' + OutputFileName + ' -p gaff \n')
atomfile.close()
os.chmod('discrenamingscript.sh', 0744)
#rename = subprocess.Popen(['./renamingscript.sh'])
#rename.wait()

    
    
