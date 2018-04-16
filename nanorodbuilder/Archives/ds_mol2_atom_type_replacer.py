#!/usr/bin/env python
import os
import itertools
# #execfile('/home/mdmannin/amber16/amber.sh')
# #subprocess.call(['source', 'home/mdmannin/amber16/amber.sh'])
batch_file_name = 'batch_input.txt'
batch_file = open(batch_file_name, 'r')
filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/'
lengthlist, radiuslist, densitylist, ligandlist,ratiolist = [], [], [], [],[]
for line in batch_file:
    line = line.strip()
    if line.find("LENGTH=") != -1:
        lengthbite = line[line.index('=')+1:].split(',')
        for i in lengthbite:lengthlist.append(int(i))
    elif line.find("RADIUS=") != -1:
        radiusbite = line[line.index('=')+1:].split(',')
        for i in radiusbite:radiuslist.append(int(i))
    elif line.find("DENSITY=") != -1:
        densitybite = line[line.index('=')+1:].split(',')
        for i in densitybite:densitylist.append(int(i))
    elif line.find("LIGAND1=") != -1:
        ligandbite = line[line.index('=')+1:].split(',')
        for i in ligandbite:ligandlist.append(i)
    elif line.find("RATIO=") != -1:
        ratiobite = line[line.index('=')+1:].split(',')
        for i in ratiobite:ratiolist.append(int(100*float(i)))
atomfile = open('renamingscript.sh', 'w+')
atomfile.write('#!/bin/sh\n')
for length, radius, density, ligand, ratio in itertools.product(lengthlist, radiuslist, densitylist, ligandlist, ratiolist):
    mol2_name = ('{}x{}-{}ratio{}dense{}'.format(length, radius, ratio, density, ligand))
    InputFileName = '{}dsmol2/{}.mol2'.format(filepath,mol2_name)
    OutputFileName = '{}Atomtype/{}_rn.ac'.format(filepath,mol2_name)
    #subprocess.call(['home/mdmannin/amber16/bin/atomtype', '-i ' + InputFileName, '-f mol2',  '-o ' + OutputFileName, '-p gaff'])
    atomfile.write('atomtype -i ' + InputFileName + ' -f mol2 -o ' + OutputFileName + ' -p gaff \n')
atomfile.close()
os.chmod('renamingscript.sh', 0744)
#rename = subprocess.Popen(['./renamingscript.sh'])
#rename.wait()

    
    