import glob
import itertools
import math
import os
import random
import sys

import numpy as np

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

RNA = True
# Avogadro's Number (Reduced to avoid memory overflow. Compensated for in volume conversion)
NA = 6.023
cellbuff = 10
bp = 30
n = 60
rnarad = 8

rnafile = '/home/mdmannin/Desktop/1aoi/RNA30.pdb'

# rnafile = '/home/mdmannin/Desktop/1aoi/RNA{}.pdb'.format(bp)
dnafile = '/home/mdmannin/Desktop/1aoi/1aoi_{}_heat.pdb'.format(bp)

if RNA == True:
    rise = 2.6
    nafile = rnafile
    na = 'rna'
if RNA == False:
    rise = 3.4
    nafile = dnafile
    na = 'dna'
if RNA == True:
    # leappath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/rna_test_ligands/Thiols/'
    # input_file_path = '/home/mdmannin/Desktop/Nanoparticles/Ligands/rna_test_ligands/Thiols/'
    leappath = '/home/mdmannin/Desktop/rna_paper/Ligands/'
    input_file_path = '/home/mdmannin/Desktop/rna_paper/Ligands/'
if RNA == False:
    leappath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/dna_test_ligands/'
    input_file_path = '/home/mdmannin/Desktop/Nanoparticles/Ligands/dna_test_ligands/'


def tleapcon(filelist):
    print('Starting tleap scripts.')
    for filename in filelist:
        print(filename)
        subprocess.check_call('tleap -f {}'.format(filename), shell=True)
    return


### Need to build separate functions for: Ligands+RNA and Ligands Only.
def ligand_mover(move, leap_script, ligand_alias, mol2_file_name):
    m = randRotate()
    translate = np.random.uniform(low=-1, high=1, size=3) * move
    # print translate
    translate = np.add(translate, np.array([rnarad * np.sign(translate[0]), rnarad * np.sign(translate[1]), 0]))
    translate_string = np.array2string(translate, precision=3, suppress_small=True, separator=' ',
                                       formatter={'float_kind': lambda x: "%.3f" % x})
    translate_string = translate_string.replace('[', '')
    translate_string = translate_string.replace(']', '')
    # copyname = 'system{}'.format(str(lig_num))
    leap_script.write('{} = loadmol2 {}\n'.format(ligand_alias, mol2_file_name))
    leap_script.write(
        'transform {} {{ {{ {} {} {} }} {{ {} {} {} }} {{ {} {} {} }} }}\n'.format(ligand_alias, m[0, 0], m[0, 1],
                                                                                   m[0, 2], m[1, 0], m[1, 1], m[1, 2],
                                                                                   m[2, 0], m[2, 1], m[2, 2]))
    leap_script.write('translate {} {{  {}  }}\n'.format(ligand_alias, translate_string))
    # leap_script.write('ligand{} = {}.1 \n'.format(lig_num, copyname))
    # leap_script.write('desc ligand{}\n'.format(i))
    # leap_script.write('remove {} ligand{}\n'.format(copyname, lig_num))
    # leap_script.write('add ligunit ligand{}\n'.format(lig_num))
    # leap_script.write('desc ligunit\n')
    return leap_script


def sysSalter(leapfile, salt, n, bp, system, cellbuff):
    leapfile.write('solvateBox {} TIP3PBOX {}\n'.format(system, cellbuff))
    leapfile.write('addions {} Cl- 0\n'.format(system))
    leapfile.write('addions {} Na+ 0\n'.format(system))
    leapfile.write('addions {} Na+ {}\n'.format(system, salt))
    leapfile.write('addions {} Cl- 0\n'.format(system))
    return


def ligandRotate(coords):
    M = randRotate()
    return np.matmul(M, coords)


def randRotate():
    angle = (random.random() - 0.5) * (2 * math.pi)
    direction = np.random.random(3) - 0.5

    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = unit_vector(direction[:3])
    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array([[0.0, -direction[2], direction[1]],
                   [direction[2], 0.0, -direction[0]],
                   [-direction[1], direction[0], 0.0]])
    M = np.identity(4)
    M[:3, :3] = R
    return np.array(M)


def unit_vector(data, axis=None, out=None):
    if out is None:
        data = np.array(data, dtype=np.float64, copy=True)
        if data.ndim == 1:
            data /= math.sqrt(np.dot(data, data))
            return data
    else:
        if out is not data:
            out[:] = np.array(data, copy=False)
        data = out
    length = np.atleast_1d(np.sum(data * data, axis))
    np.sqrt(length, length)
    if axis is not None:
        length = np.expand_dims(length, axis)
    data /= length
    if out is None:
        return data


def builder(ligandfile, leapfile, ligandname, saltconc, volume, n, bp, na):
    ### File names
    parmpath = '{}ParmFiles/'.format(leappath)
    savefile = '{}{}-{}{}bp-{}M'.format(parmpath, ligandname, na, bp, str(saltconc).replace('.', ''))

    ### Calculations
    LigandVolume = volume * n
    CellHeight = bp * rise
    CellSideL = np.sqrt(LigandVolume / (CellHeight)) + 2 * rnarad
    cellvolume = (CellSideL + 2 * cellbuff) ** 2 * (CellHeight + 2 * cellbuff)
    move = np.array([CellSideL / 2, CellSideL / 2, bp * rise / 2])
    salt = saltconc * NA * (cellvolume) * 10 ** (-4)
    print('Molarity: ', saltconc)
    print('Number of ions: ', salt)
    print('Expected Volume :', cellvolume)
    ### Load ligand file.
    # leapfile = open(leapfilename,'w')
    leapfile.write('ligand = loadMol2 {}\n'.format(ligandfile))
    leapfile.write('ligunit = createUnit RNA-Ligs\n')

    ### Single Ligand Sim
    if volume == 0:
        sysSalter(leapfile, salt, n, bp, 'ligand', cellbuff)
        leapfile.write('saveamberparm ligand {}.prmtop {}.rst7\n'.format(savefile, savefile))
        leapfile.close()
        return

    ### Copy Ligands and Randomly Move
    # print 'Number of ligands being placed: {}.\n'.format(n)
    for i in range(n):
        ligand_mover(move, leapfile)

    ### Add RNA to System.
    leapfile.write(
        'rna_seq = { U U C U A C C A A A A G U G U A U U U G G A A A C U G C U C G A G C A G U U U C C A A A U A C A C U U U U G G U A G A A } \n')
    leapfile.write('rna = loadPdbUsingSeq {} rna_seq \n'.format(rnafile))
    # leapfile.write('rna = loadpdb {}\n'.format(nafile))
    leapfile.write('alignaxes rna\n')
    leapfile.write('system = combine { ligunit rna }\n')
    leapfile.write('check system\n')
    # leapfile.write('desc system\n')
    ### Solvate and Salt.
    # leapfile.write('ExpectedVolume = {}\nlist ExpectedVolume\n'.format(cellvolume))
    sysSalter(leapfile, salt, n, bp, 'system', cellbuff)
    ### Save topology file.
    leapfile.write('saveamberparm system {}.prmtop {}.rst7 \n'.format(savefile, savefile))
    leapfile.write('quit')
    # leapfile.write('clearVariables {ligand* combined} \n')
    leapfile.close()
    return


def main():
    # Spacer for RNA/DNA molecule
    # '''
    # def rotation(theta, R = np.zeros((3,3))):
    # cx,cy,cz = np.cos(theta)
    # sx,sy,sz = np.sin(theta)
    # R.flat = (cx*cz - sx*cy*sz, cx*sz + sx*cy*cz, sx*sy,
    # -sx*cz - cx*cy*sz, -sx*sz + cx*cy*cz,
    # cx*sy, sy*sz, -sy*cz, cy)
    # return R
    # '''
    ### Set parameters for run.
    vollist = [1500]  # Cubic angstrom per ligand
    saltlist = [0.1]
    # Define path and get list of mol2 files in folder.

    ligandfilelist = glob.glob('{}*.mol2'.format(input_file_path))
    filelist = []
    for ligandfile, volume, salt in itertools.product(ligandfilelist, vollist, saltlist):
        ligandname = ligandfile.replace(input_file_path, '').replace('.mol2', '')
        leapfilename = '{}{}.in'.format(leappath, ligandname)
        print(leapfilename)
        filelist.append(leapfilename)
        leapfile = open('{}'.format(leapfilename), 'w+')
        leapfile.write('source /home/mdmannin/amber16/dat/leap/cmd/oldff/leaprc.ff14SB \n')
        leapfile.write('loadamberparams gaff.dat\n')
        leapfile.write('loadamberparams AuNP.frcmod\n')
        if RNA == True:
            # leapfile.write('source leaprc.RNA.OL3 \n')
            leapfilename = '{}{}-{}M-RNA.in'.format(leappath, ligandname, salt)
        if RNA == False:
            leapfilename = '{}{}-{}M-DNA.in'.format(leappath, ligandname, salt)
        # leapfile.write('source leaprc.DNA.OL15 \n')
        leapfile.write('loadamberparams /home/mdmannin/amber16/dat/leap/parm/frcmod.ionsjc_tip3p \n')
        # leapfile.write('loadamberparams {}frcmod.known\n'.format(input_file_path))
        builder(ligandfile, leapfile, ligandname, salt, volume, n, bp, na)
        tleapcon(filelist)


if __name__ == "__main__":
    main()
