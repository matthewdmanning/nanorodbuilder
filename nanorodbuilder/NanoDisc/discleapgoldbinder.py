#!/usr/bin/env python
#import subprocess, os
import itertools

from Archives import leapgoldbinder as leap
from numpy import int

bondcutoff = float(3)

batchfilename = 'batch_input_disc.txt'
batchfile = open(batchfilename, 'r')

filepath = '/home/mdmannin/Desktop/Nanoparticles/Nanodisc/'
parmpath = '{}ParmFiles/'.format(filepath)
print '\n\n INPUT FILES COMING FROM RENAMED MOL2. MAKE SURE THEY HAVE BEEN CREATED.\n\n'

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

def main():
    for thick, radius, curve, density, ligand, ratio in itertools.product(thicklist, radiuslist, curvelist, densitylist, ligandlist, ratiolist):
        #try:
        goldname = ('{}x{}-{}curve-{}ratio{}dense{}'.format(thick, radius, curve, ratio, density, ligand))
        goldfilename = ('{}Renamed/{}_rn'.format(filepath, goldname))
        leapshellname = '{}LeapScripts/goldbinder{}.sh'.format(filepath, goldname)
        parmpath = '{}ParmFiles/'.format(filepath)
        print 'Input file: {}.mol2.'.format(goldfilename)
        print 'Saving script to {}.'.format(leapshellname)
        leap.leapscript(filepath, goldname, goldfilename, leapshellname)
        
if __name__ == "__main__":
    main()
#===============================================================================
# 
#     #===========================================================================
#     def overlapcheck(atoms, goldid, protected, goldsulflig):
#         atomsdelete = []
#         ausurvivors = []
#         atomsize = atoms.shape[0]
#         for counter, atom in enumerate(atoms):
#             if counter in atomsdelete or residuelist[int(atom[0])] in residuedelete: continue
#             for i in range(counter+1, atomsize - 1):
#                 checkatom = atoms[i]
#                 if i in atomsdelete or residuelist[int(checkatom[0])] in residuedelete: continue
#                 if residuelist[int(atom[0])] == residuelist[int(checkatom[0])] and int(atom[0]) not in goldid: continue
#                 if math.sqrt((atom[1] - checkatom[1])**2  + (atom[2] - checkatom[2])**2 + (atom[3] - checkatom[3])**2) - 1.0 < 0.001 or (atom[1] == checkatom[1] and atom[2] == checkatom[2] and atom[3] == checkatom[3]):
#                     if int(atom[0]) in goldid and int(checkatom[0]) not in goldid:
#                         #print 'Overlapping ligand containing Atom #{}.'.format(atom[0])
#                         residuedelete.append(residuelist[int(checkatom[0])])
#                         print 'Overlapping ligand found. Deleting ligand residue: gold.{}.\n'.format(residuelist[int(checkatom[0])])
#                     elif int(atom[0]) not in goldid and int(checkatom[0]) in goldid:
#                         residuedelete.append(residuelist[int(atom[0])])
#                         print 'Overlapping ligand found. Deleting ligand residue: gold.{}.\n'.format(residuelist[int(atom[0])])
#                     elif int(atom[0]) in goldid and int(checkatom[0]) in goldid:
#                         ###Delete overlapping atom if not bound to ligand.
#                         if int(checkatom[0]) not in protected:
#                             atomsdelete.append(i)
#                             ausurvivors.append(counter)
#                             #print atom
#                             #print checkatom
#                             #print 'Gold overlap deleted ({})'.format(namelist[i])
#                         ###Atom being check is overlapping with a ligand-bound atom.
#                         elif int(atom[0]) not in protected:
#                             atomsdelete.append(counter)
#                             ausurvivors.append(counter)
#                             #print 'Gold overlap deleted ({})'.format(namelist[counter])
#                             break
#                         ### Both atoms are bound to ligands. -> Not good. Choose one gold/ligand and delete.
#                         else:
#                             row = np.ndarray.tolist(goldsulflig[:,0]).index(checkatom[0])
#                             print 'Overlapping ligand-bound Au atoms #{} and #{}'.format(namelist[counter], namelist[i])
#                             print 'Deleting residue gold.{} and Au gold.{}.Au{}'.format(goldsulflig[row,2], residuelist[int(checkatom[0])], residuelist[int(checkatom[0])])
#                             atomsdelete.append(i)
#                             ausurvivors.append(counter)
#                             residuedelete.append(goldsulflig[row,2])
# 
#         print ausurvivors
#         #print 'Counter = {}. protected = {}'.format(counter, protected)
#         return 1,atomsdelete, residuedelete, ausurvivors;
# 
# def distancecheck(coords1, coords2, atomid1, atomid2):
#     x1, y1, z1 = coords1
#     x2, y2, z2 = coords2
#     distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
#     #if distance < bondcutoff: print 'Coords: {},{},{}; {},{},{}.  '.format(x1,y1,z1,x2,y2,z2)
#     #print distance
#     #print distance<bondcutoff
#     #if distance < 1.0:
#     #    print 'Missed deletion: {} and {}.'.format(namelist[int(atomid1)], namelist[int(atomid2)])
#     return distance
# 
# # ##Replace atom names in mol2 file.
# for thick, radius, curve, density, ligand, ratio in itertools.product(thicklist, radiuslist, curvelist, densitylist, ligandlist, ratiolist):
#     #try:
#     goldname = ('{}x{}-{}curve-{}ratio{}dense{}'.format(thick, radius, curve, ratio, density, ligand))
#     goldfilename = ('{}Renamed/{}_rn'.format(filepath, goldname))
#     leapshellname = '{}LeapScripts/goldbinder{}.sh'.format(filepath, goldname)
#     parmpath = '{}ParmFiles/'.format(filepath)
#     print 'Input file: {}.mol2.'.format(goldfilename)
#     print 'Saving script to {}.'.format(leapshellname)
#     leap.leapscript(filepath, goldname, goldfilename, leapshellname)
#     goldfile = open('{}.mol2'.format(goldfilename), 'r')
#     trigger = 0
#     bondtrigger = 0
#     bondnum = 0
#     atomslist = []
#     goldidlist = []
#     sulfidlist = []
#     bondlist = []
#     protected = []
#     #List of atoms with eight bonds (AMBER max).
#     fullbond = []
#     neighborlist = []
#     ausurvivors = []
#     #List of gold names in original mol2
#     namelist = []
#     #Contains residues corresponding to atoms. Do not modify.
#     residuelist = []
#     goldres = []
#     #List of deleted residues.
#     atomsdelete = []
#     residuedelete = []
#     #Numpy arrays holding [atom #, residue #]
#     goldsulfliglist = []
#     #List for atoms names of all gold in file.
#     namelist = []
#     sulfname = {'resnum':'sulfname'}
#         # ## Array format: Atom ID #, XYZ coords, # of bonds, ID#'s of bonded atoms
#     for line in goldfile:
#         spline = line.split()
#         if '@<TRIPOS>ATOM' in line:
#             trigger = 1
#             continue
#         if '@<TRIPOS>BOND' in line:
#             #print atoms
#             bondtrigger = 1
#             trigger = 0
#             goldid = np.vstack(goldidlist)
#             sulfid = np.vstack(sulfidlist)
#             atoms = np.vstack(atomslist)
#             goldsulflig = np.vstack(goldsulfliglist)
#             goldresidnum = np.unique(np.vstack(goldres))
#             goldresstr = []
#             for x in goldresidnum: goldresstr.append("gold.{}".format(x))
#             continue
#         if trigger == 1:
#             atomslist.append([int(spline[0])-1, float(spline[2]), float(spline[3]), float(spline[4]),0, 0, 0, 0, 0, 0, 0, 0, 0])
#             residuelist.append(int(spline[6]))
#             namelist.append(spline[1])
#             if 'Au' in line:
#                 goldidlist.append(int(spline[0])-1)
#                 goldres.append(int(spline[6]))
#             elif 'S' in spline[1]:
#                 sulfidlist.append(int(spline[0])-1)
#                 goldsulfliglist.append([0, int(int(spline[0])-1), int(spline[6])])
#                 sulfname.update({int(spline[6]):spline[1]})
#         elif bondtrigger == 1:
#             # ## Assemble list of bonds and update atom lists                                
#             if '@' in line:
#                 break
#             bondnum = int(spline[0])
#             # ## Place Au in protected list if it is bound to S
#             atom1 = int(int(spline[1]) - 1)
#             atom2 = int(int(spline[2]) - 1)
#             if atom1 in goldid and atom2 in sulfid:
#                 protected.append(atom1)
#                 row = np.ndarray.tolist(goldsulflig[:,1]).index(atom2)
#                 goldsulflig[row,0] = atom1    
#             if atom2 in goldid and atom1 in sulfid:
#                 protected.append(atom2)
#                 row = np.ndarray.tolist(goldsulflig[:,1]).index(atom1)
#                 goldsulflig[row,0] = atom2    
#             atom1bonds = int(atoms[atom1, 4])
#             atom2bonds = int(atoms[atom2, 4])
#             atoms[atom1, 5 + atom1bonds] = atom2
#             atoms[atom2, 5 + atom2bonds] = atom1
#             atoms[atom1, 4] = int(atoms[atom1, 4]) + 1
#             atoms[atom2, 4] = int(atoms[atom2, 4]) + 1
#     success, atomsdelete, residuedelete, ausurvivors = overlapcheck(atoms, goldid, protected, goldsulflig)
#      
#     #if success == 0: continue
#     atomsize = atoms.shape[0]
#     ###Bonding of Gold. No need to check for deleted residues, only atoms.
#     for counter in range(atomsize - 1):
#         if counter not in goldid: continue
#         atom = atoms[counter]
#         neighborlist[:] = []
#         if int(atom[0]) in atomsdelete or int(atom[0]) in fullbond: continue
#         for i in range(counter+1, atomsize - 1):
#             if i not in goldid: continue
#             checkatom = atoms[i]
#             if  int(checkatom[0]) in atomsdelete or int(checkatom[0]) in fullbond: continue
#             elif int(checkatom[4]) >= 8:
#                 fullbond.append(int(checkatom[0]))
#                 continue
#             dist = distancecheck(atom[1:4], checkatom[1:4], int(atom[0]), int(checkatom[0]))
#             if dist - 1.0 < 0.001:
#                 atomsdelete.append(i)
#                 ausurvivors.append(counter)
#                 #leapshell.write('remove gold.{} gold.{}.{}\n'.format(residuelist[int(checkatom[0])], residuelist[int(checkatom[0])],namelist[int(checkatom[0])]))
#                 #print'Deleting atom: {}.\n'.format(namelist[int(checkatom[0])])
#             elif bondcutoff - dist > 0.001:
#                 neighborlist.append(int(checkatom[0]))
#         #print 'Neighbor list of Atom {} is: {}.\n'.format(counter+1,neighborlist)
#         random.shuffle(neighborlist)
#         for i in range(len(neighborlist)-1):
#             bondid = int(neighborlist[i])
#             if atoms[counter, 4] >= 8: 
#                 fullbond.append(int(atom[0]))
#                 break
#             bondnum += 1
#             #leapshell.write('bond gold.{}.{} gold.{}.{}\n'.format(residuelist[int(atom[0])], namelist[int(atom[0])], residuelist[bondid],  namelist[bondid]))
#             #bondlist.append(['  {} {} {} 1\n'.format(bondnum,int(atom[0]+1),int(bondid+1))])
#             bondlist.append([int(atom[0]), bondid])
#             atoms[counter, int(atoms[counter,4]+5)] = bondid
#             atoms[bondid, int(atoms[bondid, 4]+5)] = atom[0]
#             atoms[bondid, 4] = atoms[bondid, 4] + 1
#             atoms[counter, 4] = atoms[counter, 4] + 1
#     namesdelete = []
#     for x in atomsdelete: namesdelete.append(namelist[x])
#     #print 'Atoms deleted: {}'.format(", ".join(namesdelete))
#     #print 'Residues deleted: {}'.format(", ".join(residuedelete))
#     ### Write Leap Input File
#     leapshell = open(leapshellname,'w+')
#     leapshell.write('loadamberparams /home/mdmannin/amber16/dat/leap/parm/AuNP.frcmod\n')
#     leapshell.write('loadamberparams gaff.dat\n')
#     leapshell.write('gold = loadMol2 {}.mol2\n'.format(goldfilename))
#     for residue in residuedelete: leapshell.write('remove gold gold.{}\n'.format(residue))
#     #leapshell.write('goldatoms = createResidue goldatoms\n')
#     #leapshell.write('add gold goldatoms\n')
#     goldnum = 0
#     #aliasdict = {'counter':'alias'}
#     for counter, atom in enumerate(atomslist):
#         if counter not in goldid: continue
#         #Deletes gold atoms on deletelist.
#         if counter in atomsdelete or residuelist[int(atom[0])] in residuedelete:
#             leapshell.write('remove gold.{} gold.{}.{}\n'.format(residuelist[counter], residuelist[counter], namelist[int(atom[0])]))
#             continue
#         ### Alias all gold atoms for reference after move to new residue.
#         #alias = 'gold{}'.format(goldnum)
#         #aliasdict.update({counter:alias})
#         #leapshell.write('{} = gold.{}.{}\n'.format(alias, residuelist[counter], namelist[int(atom[0])]))
#         #leapshell.write('set {} name "g{}"\n'.format(alias, goldnum))
#         #leapshell.write('remove gold.{} gold.{}.g{}\n'.format(residuelist[counter], residuelist[counter],goldnum))
#         #leapshell.write('add goldatoms {}\n'.format(alias))
#         #if counter in protected:
#         #    row = np.ndarray.tolist(goldsulflig[:,0]).index(int(atom[0]))
#         #    sulfres = goldsulflig[row,2]
#         #    leapshell.write('bond {} gold.{}.{}\n'.format(alias, sulfres, sulfname[sulfres]))
#         #goldnum += 1
#     #for x in goldresidnum: leapshell.write('remove gold gold.{}\n'.format(x))
#     # Bond gold atoms to neighbors.
#     for bonding in bondlist:
#         leapshell.write('bond gold.{}.{} gold.{}.{}\n'.format(residuelist[bonding[0]], namelist[bonding[0]],residuelist[bonding[1]], namelist[bonding[1]]))
#     leapshell.write('saveMol2 gold {}BondedMol2/{}_bond.mol2 1\n'.format(filepath,goldname))
#     leapshell.write('saveamberparm gold {}{}_vac.prmtop {}{}_vac.rst7\n'.format(parmpath,goldname, parmpath, goldname))
#     leapshell.write('check gold\n')
#     leapshell.write('quit\n')
#     leapshell.close()
#     print 'Tleap input file ready!.\n'
#     goldfile.close()    
#     print 'Number of ligands: {}\n'.format(sulfid.shape[0])
#     print 'Number of deleted atoms: {}\n\n.'.format(len(atomsdelete))
#                 
# ### Create master file for tLeap execution.
# shellpath = '{}LeapScripts/'.format(filepath)
# mastershellfile = open('{}masterscript.sh'.format(shellpath), 'w')
# #mastershellfile = open('/home/mdmannin/Dropbox/NanorodScripts/masterleapscript.sh','w')
# for thick, radius, curve, density, ligand, ratio in itertools.product(thicklist, radiuslist, curvelist, densitylist, ligandlist, ratiolist):
#     iterfile = '{}goldbinder{}x{}-{}curve-{}ratio{}dense{}'.format(shellpath, thick, radius, curve, ratio, density, ligand)
#     mastershellfile.write('tleap -f {}.sh > {}.log'.format(iterfile, iterfile))
#     mastershellfile.write('sleep 5')
# mastershellfile.close()
# #os.chmod('/~/Dropbox/NanorodScripts/masterleapscript.sh', 0744)
# print 'Master shell script created: {}.\n'.format(shellpath)
#===============================================================================
