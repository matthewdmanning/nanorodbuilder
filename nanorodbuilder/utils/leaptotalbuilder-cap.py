#!/usr/bin/env python
import itertools
import math
import random

import numpy as np
from numpy import int

from nanorodbuilder.utils import nanorodbatchrefactor

bondcutoff = float(3.0)
solventbuffer = int(10)
overlapdist = 1.5
throw = 1.5
### Outline
# 1. Fetch dimensions and files (gold mol2 files (surf and core) and ligand mol2).
# 2. Check for overlap and delete non-ligand gold.
# 3. Open shell file and import parameters and FFs. (Load frcmod files from RESP output)
# 3. Throw ligands, bind to gold, and 'impose bond length'.
# 4. Bind gold atoms.
#    Save mol2? or prmtop? file and use antechamber to rename/combine gold residues?
# 5. Solvate/salt and build topology

# Merging nanorodbatch.py
# Leap Commands after confirming part of nanorod.
# system = createUnit system?
# AU = createResidue
# createAtom name ligand_type(Au) charge
# ??? set ATOM element Au
# set ATOM position { x y z }

filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/'
parmpath = '{}ParmFiles/'.format(filepath)
print
'\n\n INPUT FILES COMING FROM RENAMED MOL2. MAKE SURE THEY HAVE BEEN CREATED.\n\n'


# ===========================================================================
# Check ligands for overlaps and add to list of deleted atoms.
def overlapcheck(atoms, goldid, protected, goldsulflig, residuelist, namelist, residuedelete):
    atomsdelete = []
    ausurvivors = []
    atomsize = atoms.shape[0]
    for counter, atom in enumerate(atoms):
        if counter in atomsdelete or residuelist[int(atom[0])] in residuedelete: continue
        for i in range(counter + 1, atomsize - 1):
            checkatom = atoms[i]
            if i in atomsdelete or residuelist[int(checkatom[0])] in residuedelete: continue
            if residuelist[int(atom[0])] == residuelist[int(checkatom[0])] and int(atom[0]) not in goldid: continue
            if math.sqrt((atom[1] - checkatom[1]) ** 2 + (atom[2] - checkatom[2]) ** 2 + (
                    atom[3] - checkatom[3]) ** 2) - 1.0 < overlapdist or (
                    atom[1] == checkatom[1] and atom[2] == checkatom[2] and atom[3] == checkatom[3]):
                if int(atom[0]) in goldid and int(checkatom[0]) not in goldid:
                    # print 'Overlapping ligand containing Atom #{}.'.format(atom[0])
                    residuedelete.append(residuelist[int(checkatom[0])])
                    print
                    'Overlapping ligand found. Deleting ligand residue: nanorod.{}.\n'.format(
                        residuelist[int(checkatom[0])])
                elif int(atom[0]) not in goldid and int(checkatom[0]) in goldid:
                    residuedelete.append(residuelist[int(atom[0])])
                    print
                    'Overlapping ligand found. Deleting ligand residue: nanorod.{}.\n'.format(residuelist[int(atom[0])])
                elif int(atom[0]) in goldid and int(checkatom[0]) in goldid:
                    ###Delete overlapping atom if not bound to ligand.
                    if int(checkatom[0]) not in protected:
                        atomsdelete.append(i)
                        ausurvivors.append(counter)
                        # print atom
                        # print checkatom
                        # print 'Gold overlap deleted ({})'.format(namelist[i])
                    ###Atom being check is overlapping with a ligand-bound atom.
                    elif int(atom[0]) not in protected:
                        atomsdelete.append(counter)
                        ausurvivors.append(i)
                        # print 'Gold overlap deleted ({})'.format(namelist[counter])
                        break
                    ### Both atoms are bound to ligands. -> Not good. Choose one gold/ligand and delete.
                    else:
                        row = np.ndarray.tolist(goldsulflig[:, 0]).index(int(checkatom[0]))
                        print
                        row
                        print
                        goldsulflig[row]
                        print
                        'Overlapping ligand-bound Au atoms #{} and #{}'.format(namelist[counter],
                                                                               namelist[int(checkatom[0])])
                        print
                        'Deleting residue nanorod.{} and Au nanorod.{}.{}'.format(goldsulflig[row, 2],
                                                                                  residuelist[int(checkatom[0])],
                                                                                  namelist[int(checkatom[0])])
                        atomsdelete.append(i)
                        ausurvivors.append(counter)
                        residuedelete.append(goldsulflig[row, 2])
    print
    residuedelete
    # print 'Counter = {}. protected = {}'.format(counter, protected)
    return 1, atomsdelete, residuedelete, ausurvivors;


# Calculate Cartesian distances
def distancecheck(coords1, coords2):
    x1, y1, z1 = coords1
    x2, y2, z2 = coords2
    distance = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
    return distance


def main():
    # Fetch parameters from input text files.
    def parameters():
        lengthlist, radiuslist, densitylist, ligand1list, ligand2list, capligandlist, sulfindex1list, sulfindex2list, capsulfindexlist, ratiolist = [], [], [], [], [], [], [], [], [], []
        batchfilename = 'batch_input.txt'
        batchfile = open(batchfilename, 'r')
        for line in batchfile:
            if line.find("LENGTH=") != -1:
                lengthbite = line[line.index('=') + 1:].split(',')
                for i in lengthbite: lengthlist.append(int(i.rstrip()))
            elif line.find("RADIUS=") != -1:
                radiusbite = line[line.index('=') + 1:].split(',')
                print
                radiusbite
                for i in radiusbite: radiuslist.append(int(i.rstrip()))
            elif line.find("DENSITY=") != -1:
                densitybite = line[line.index('=') + 1:].split(',')
                for i in densitybite: densitylist.append(int(i.rstrip()))
            elif line.find("LIGAND1=") != -1:
                ligand1bite = line[line.index('=') + 1:].split(',')
                for i in ligand1bite: ligand1list.append(i.rstrip())
            elif line.find("LIGAND2=") != -1:
                ligand2bite = line[line.index('=') + 1:].split(',')
                for i in ligand2bite: ligand2list.append(i.rstrip())
            elif line.find("CAPLIGAND=") != -1:
                capligandbite = line[line.index('=') + 1:].split(',')
                for i in capligandbite: capligandlist.append(i.rstrip())
            elif line.find("SULFINDEX1=") != -1:
                sulf1bite = line[line.index('=') + 1:].split(',')
                for i in sulf1bite: sulfindex1list.append(i.rstrip())
            elif line.find("SULFINDEX2=") != -1:
                sulf2bite = line[line.index('=') + 1:].split(',')
                for i in sulf2bite: sulfindex2list.append(i.rstrip())
            elif line.find("CAPSULFINDEX=") != -1:
                capsulfbite = line[line.index('=') + 1:].split(',')
                for i in capsulfbite: capsulfindexlist.append(i.rstrip())
            elif line.find("RATIO=") != -1:
                ratiobite = line[line.index('=') + 1:].split(',')
                for i in ratiobite: ratiolist.append(100 * float(i.rstrip()))
        batchfile.close()
        # for i in itertools.chain(lengthlist, radiuslist, densitylist, ligand1list, ligand2list, sulfindex1list, sulfindex2list, ratiolist):
        #    i = i.rstrip()
        #    print i
        return lengthlist, radiuslist, densitylist, ligand1list, ligand2list, capligandlist, sulfindex1list, sulfindex2list, capsulfindexlist, ratiolist

    # Create a master shell file that will execute each tLeap script generated.
    def MasterFile():
        ### Create master file for tLeap execution.
        filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/'
        parmpath = '{}ParmFiles/'.format(filepath)
        shellpath = '{}LeapScripts/'.format(filepath)
        mastershellfile = open('{}masterscript.sh'.format(shellpath), 'w+')
        for length, radius, density, ligand, ratio in itertools.product(lengthlist, radiuslist, densitylist,
                                                                        ligand1list, ratiolist):
            iterfile = '{}goldbinder{}x{}-{}ratio{}dense{}'.format(shellpath, length, radius, int(ratio), density,
                                                                   ligand)
            mastershellfile.write('tleap -s -f {}.sh > {}.log'.format(iterfile, iterfile))
            mastershellfile.write('sleep 5')
        mastershellfile.close()
        print
        'Master shell script created: {}.\n'.format(shellpath)

    # Search for nearest neighbors and add to list of bonded atoms (bondpairs)
    # Keeps track of number of bonds and avoids having more than 8.
    def goldbonder(shellarray, corearray, anchorarray):
        combinedarray = shellarray + corearray
        anchorbonds, atombonds = [], []
        for i in range(len(anchorarray)): anchorbonds.append(0)
        for i in range(len(combinedarray)):  atombonds.append(0)
        bondpairs = []
        neighborlist = []
        for ai, atom in enumerate(anchorarray):
            for ni, neigh in enumerate(shellarray):
                if neigh in anchorarray: continue
                if atombonds[ni] > 7: continue
                dist = distancecheck(atom, neigh)
                ### if dist - 1.0 < 0.001: #Check for overlapping golds.
                if bondcutoff - abs(dist) > 0.001:
                    neighborlist.append(ni)
            # print 'Neighbor list of Atom {} is: {}.\n'.format(counter+1,neighborlist)
            random.shuffle(neighborlist)
            # anchorbonds[ai] = neighborlist[0:7]
            for bondnum in range(7):
                if len(neighborlist) < bondnum + 1: break
                bondpairs.append([combinedarray.index(atom), neighborlist[bondnum]])
                atombonds[neighborlist[bondnum]] = + 1
            neighborlist[:] = []

        for i, atom in enumerate(combinedarray):
            if atom in anchorarray: continue
            if atombonds[i] > 7: continue
            numbonds = atombonds[i]
            valence = 8 - numbonds
            for n, neigh in enumerate(combinedarray):
                if n <= i: continue
                if neigh in anchorarray: continue
                if atombonds[n] > 7: continue
                dist = distancecheck(atom, neigh)
                # if dist - 1.0 < 0.0001:
                if bondcutoff - abs(dist) > 0.001:
                    neighborlist.append(n)
            # print 'Neighbor list of Atom {} is: {}.\n'.format(counter+1,neighborlist)
            random.shuffle(neighborlist)
            atombonds[i] = 8
            for bondnum in range(0, valence):
                if len(neighborlist) < bondnum + 1: break
                bondpairs.append([i, neighborlist[bondnum]])
                atombonds[neighborlist[bondnum]] = + 1
            neighborlist[:] = []
        for pair in bondpairs:
            if [pair[1], pair[0]] in bondpairs: print
            'Duplicate\n'
        return bondpairs, combinedarray

    # Write tLeap input file.
    def shellwriter(shellarray, corearray, caparray, boundarray, boundcaparray, ligand1, ligand2, capligand, sulfindex1,
                    sulfindex2, capsulfindex, ratio, filepath, goldname):

        ligandpath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/'
        allboundarray = boundarray + boundcaparray
        outerarray = shellarray + caparray
        bondpairs, combinedarray = goldbonder(outerarray, corearray, allboundarray)
        ### Load parameters and forcefields
        leapshell = open(leapshellname, 'w+')
        leapshell.write('logfile {}.log\n'.format(goldname))
        print
        goldname
        # leapshell.write('#source /home/mdmannin/amber16/dat/leap/cmd/oldff/leaprc.ff14SB\n')
        leapshell.write('#source /home/mdmannin/amber16/dat/leap/cmd/oldff/leaprc.ff99bsc0\n')
        # leapshell.write('#source /home/mdmannin/amber16/dat/leap/cmd/oldff/leaprc.ff99chiOL3\n')
        leapshell.write('loadamberparams /home/mdmannin/amber16/dat/leap/parm/AuNP.frcmod\n')
        leapshell.write('loadamberparams /home/mdmannin/amber16/dat/leap/parm/frcmod.ionsjc_tip3p\n')
        leapshell.write('loadamberparams gaff.dat\n')
        ### Create gold residue and atoms and move to coordinates in list.
        leapshell.write('nanorod = createUnit rod\n')
        leapshell.write('aures = createResidue goldres\n')
        leapshell.write('add nanorod aures\n')
        liglist, caplist = [], []
        for n, atom in enumerate(combinedarray):
            leapshell.write('Au{} = createAtom g{} Au 0\n'.format(n, n))
            leapshell.write('set Au{} element Au\n'.format(n))
            leapshell.write('set Au{} position {{ {} {} {} }}\n'.format(n, atom[0], atom[1], atom[2]))
            leapshell.write('add nanorod.1 Au{} \n'.format(n))
            if atom in boundarray:
                liglist.append(n)
            elif atom in boundcaparray:
                caplist.append(n)
        # Find gold atom pairs and write commands for bonding.
        # leapshell.write('bondByDistance nanorod {}\n'.format(bondcutoff))
        # leapshell.write('saveMol2 nanorod bondedgold.mol2 1\n')

        ### Place ligands around nanorod. Choose between charged/uncharged by random num. Bond to gold.
        print
        ligand1, ligand2
        leapshell.write('philicligand = loadmol2 {}{}.mol2\n'.format(ligandpath, ligand1))
        leapshell.write('set philicligand.1 name "{}"\n'.format(ligand1))
        leapshell.write('phobicligand = loadmol2 {}{}.mol2\n'.format(ligandpath, ligand2))
        leapshell.write('set phobicligand.1 name "{}"\n'.format(ligand2))

        cappath = '/home/mdmannin/Desktop/Photoswitching/'
        leapshell.write('capligand = loadmol2 {}{}.mol2\n'.format(cappath, capligand))
        leapshell.write('set capligand.1 name "{}"\n'.format(capligand))

        sulfname1 = 'S8'
        sulfname2 = 'S8'
        sulfdist = '-2.0'
        throw = 1.75
        ligand1num, ligand2num, capnum = 0, 0, 0

        print
        liglist
        for n, ligau in enumerate(liglist):
            randindex = random.random()
            x, y, z = combinedarray[ligau]
            throwx = x * throw
            throwy = y * throw
            throwz = z * throw
            # Copy and move ligand depending on which ligand ligand_type was chosen.
            if randindex - ratio / 100 < 0.00001:
                # leapshell, ligand1num = ligandPlacer(leapshell, combinedarray, 'philicligand', ligand1num, ligand1, n, sulfindex1, ligau)
                leapshell.write('ligand{} = copy {}\n'.format(n, 'philicligand'))
                if z < 0.0001: leapshell.write(
                    'transform ligand{}.1 {{ {{ 1 0 0 }} {{ 0 1 0 }} {{ 0 0 -1 }} }}\n'.format(n))
                leapshell.write('translate ligand{}.1 {{ {} {} {} }}\n'.format(n, throwx, throwy, 1.5 * throwz))
                leapshell.write('{}{} = ligand{}.1\n'.format(ligand1, ligand1num, n))
                leapshell.write('remove ligand{} {}{}\n'.format(n, ligand1, ligand1num))
                leapshell.write('add nanorod {}{}\n'.format(ligand1, ligand1num))
                leapshell.write('bond nanorod.{}.{} nanorod.goldres.g{}\n'.format(n + 2, sulfindex1, ligau))
                ligand1num = ligand1num + 1
                # leapshell.write('impose nanorod {{ 1 {} }} {{ {{ "g{}" "{}" {} }} }}\n'.format(n+2, ligau, sulfname1, sulfdist))
            else:
                # leapshell, ligand2num = ligandPlacer(leapshell, combinedarray, 'phobicligand', ligand2num, ligand2, n, sulfindex2, ligau)
                leapshell.write('ligand{} = copy {}'.format(n, 'phobicligand\n'))
                if z < 0.0001: leapshell.write(
                    'transform ligand{}.1 {{ {{ 1 0 0 }} {{ 0 1 0 }} {{ 0 0 -1 }} }}\n'.format(n))
                leapshell.write('translate ligand{}.1 {{ {} {} {} }}\n'.format(n, throwx, throwy, 1.5 * throwz))
                leapshell.write('{}{} = ligand{}.1\n'.format(ligand2, ligand2num, n))
                leapshell.write('remove ligand{} {}{}\n'.format(n, ligand2, ligand2num))
                leapshell.write('add nanorod {}{}\n'.format(ligand2, ligand2num))
                leapshell.write('bond nanorod.{}.{} nanorod.goldres.g{}\n'.format(n + 2, sulfindex2, ligau))
                ligand2num = ligand2num + 1
                # leapshell.write('impose nanorod {{ 1 {} }} {{ {{ "g{}" "{}" {} }} }}\n'.format(n+2, ligau, sulfname2, sulfdist))
        print
        caplist
        for n, ligau in enumerate(caplist):
            n = n + np.size(liglist)
            x, y, z = combinedarray[ligau]
            throwx = x * throw
            throwy = y * throw
            throwz = z * throw
            # leapshell, capnum = ligandPlacer(leapshell, combinedarray, 'capligand', capnum, capligand, n, capsulfindex, ligau)
            leapshell.write('ligand{} = copy {}'.format(n, 'capligand\n'))
            if z < 0.0001: leapshell.write(
                'transform ligand{}.1 {{ {{ 1 0 0 }} {{ 0 1 0 }} {{ 0 0 -1 }} }}\n'.format(n))
            leapshell.write('translate ligand{}.1 {{ {} {} {} }}\n'.format(n, throwx, throwy, 1.5 * throwz))
            leapshell.write('{}{} = ligand{}.1\n'.format(capligand, capnum, n))
            leapshell.write('remove ligand{} {}{}\n'.format(n, capligand, capnum))
            leapshell.write('add nanorod {}{}\n'.format(capligand, capnum))
            leapshell.write('bond nanorod.{}.{} nanorod.goldres.g{}\n'.format(n + 2, capsulfindex, ligau))
            capnum = capnum + 1

        print
        ligand1num, ligand2num, capnum
        # leapshell.write('relax nanorod') Can't get this to work.
        for pair in bondpairs: leapshell.write('bond nanorod.1.g{} nanorod.1.g{}\n'.format(pair[0], pair[1]))

        ligandlist = ''
        for lig in range(len(liglist) - 1): ligandlist = ligandlist + 'ligand{} '.format(lig)
        # for n, ligau in enumerate(liglist):

        # for residue in residuedelete: leapshell.write('remove gold gold.{}\n'.format(residue))
        # aliasdict = {'counter':'alias'}
        ### Alias all gold atoms for reference after move to new residue.
        # Bond gold atoms to neighbors.
        # for bonding in bondlist:
        #    leapshell.write('bond gold.{}.{} gold.{}.{}\n'.format(residuelist[bonding[0]], namelist[bonding[0]],residuelist[bonding[1]], namelist[bonding[1]]))
        leapshell.write('saveMol2 nanorod {}BondedMol2/{}_bond.mol2 1\n'.format(filepath, goldname))
        leapshell.write('saveMol2 nanorod {}dsmol2/{}_bond.mol2 0\n'.format(filepath, goldname))
        leapshell.write(
            'saveamberparm nanorod {}{}_vac.prmtop {}{}_vac.rst7\n'.format(parmpath, goldname, parmpath, goldname))
        # leapshell.write('solvateBox nanorod TIP3PBOX {}\n'.format(int(solventbuffer)))
        # leapshell.write('saveamberparm nanorod {}{}_vac.prmtop {}{}_vac.rst7\n'.format(parmpath,goldname, parmpath, goldname))
        # leapshell.write('check nanorod\n')
        # leapshell.write('quit\n')
        leapshell.close()
        print
        'Tleap input file ready!.\n'
        return

    def ligandPlacer(leapshell, combinedarray, ligandname, ligandnum, ligandoriginal, n, sulfindex, goldindex):
        x, y, z = combinedarray[goldindex]
        throwx = x * throw
        throwy = y * throw
        throwz = z * throw
        leapshell.write('ligand{} = copy {}\n'.format(n, ligandoriginal))
        if z < 0.0001: leapshell.write('transform ligand{}.1 {{ {{ 1 0 0 }} {{ 0 1 0 }} {{ 0 0 -1 }} }}\n'.format(n))
        leapshell.write('translate ligand{}.1 {{ {} {} {} }}\n'.format(n, throwx, throwy, 1.5 * throwz))
        leapshell.write('{}{} = ligand{}.1\n'.format(ligandname, ligandnum, n))
        leapshell.write('remove ligand{} {}{}\n'.format(n, ligandname, ligandnum))
        leapshell.write('add nanorod {}{}\n'.format(ligandname, ligandnum))
        leapshell.write('bond nanorod.{}.{} nanorod.goldres.g{}\n'.format(n + 2, sulfindex, goldindex))
        ligandnum = ligandnum + 1
        return leapshell, ligandnum

    ### Main Code
    lengthlist, radiuslist, densitylist, ligand1list, ligand2list, capligandlist, sulfindex1list, sulfindex2list, capsulfindexlist, ratiolist = parameters()
    # print lengthlist, radiuslist, densitylist, ligandlist, ratiolist
    print
    list(itertools.product(lengthlist, radiuslist, densitylist, enumerate(ligand1list), ratiolist))

    # Nested loops iterate over all combinations of geometries and ligand choices. Nesting minimizes number of times xyzbuilder and ligandpicker are called.
    for length, radius in itertools.product(lengthlist, radiuslist):
        parmpath = '{}ParmFiles/'.format(filepath)
        shellarray, corearray, normarray, diagarray, caparray = nanorodbatchrefactor.xyzbuilder(length, radius)
        # shellarray = shellarray.append(caparray)
        for density in densitylist:
            ### Choose ligand sites.
            boundarray = nanorodbatchrefactor.ligandpicker(length, radius, density, shellarray, filepath)
            boundcaparray = nanorodbatchrefactor.ligandpicker(length, radius, density, caparray, filepath)
            for ratio, captup, ligandtup in itertools.product(ratiolist, enumerate(capligandlist),
                                                              enumerate(ligand1list)):
                ligandindex = ligandtup[0]
                capindex = captup[0]
                ligand1 = ligand1list[ligandindex]
                ligand2 = ligand2list[ligandindex]
                capligand = capligandlist[capindex]
                sulfindex1 = sulfindex1list[ligandindex]
                sulfindex2 = sulfindex2list[ligandindex]
                capsulfindex = capsulfindexlist[capindex]
                goldname = (
                    '{}x{}-{}ratio{}dense{}-cap{}'.format(length, radius, int(ratio), density, ligand1, capligand))
                leapshellname = '{}LeapScripts/goldbinder{}.in'.format(filepath, goldname)
                print
                'Saving script to {}.'.format(leapshellname)
                shellwriter(shellarray, corearray, caparray, boundarray, boundcaparray, ligand1, ligand2, capligand,
                            sulfindex1, sulfindex2, capsulfindex, ratio, filepath, goldname)  # MasterFile()


if __name__ == "__main__":
    main()
