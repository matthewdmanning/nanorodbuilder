#!/usr/bin/env python
import itertools
# #execfile('/home/mdmannin/amber16/amber.sh')
# #subprocess.call(['source', 'home/mdmannin/amber16/amber.sh'])
batchfilename = 'batch_input_disc.txt'
batchfile = open(batchfilename, 'r')
infilepath = '/home/mdmannin/Desktop/Nanoparticles/Nanodisc/discmol2/'
outfilepath = '/home/mdmannin/Desktop/Nanoparticles/Nanodisc/Atomtype/'
renamedfilepath = '/home/mdmannin/Desktop/Nanoparticles/Nanodisc/Renamed/'

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

    #===========================================================================
    # elif line.find("LIGANDS2=") == 1:
    #     ligand2list = line[line.index('=')+1,:].split(delimiter=',')
    # elif line.find("SULF1INDEX=") == 1:
    #     sulf1indexlist = line[line.index('=')+1,:].split(delimiter=',')
    # elif line.find("SULF2INDEX=") == 1:
    #     sulf2indexlist = line[line.index('=')+1,:].split(delimiter=',')  
    #===========================================================================
#===============================================================================
# lengthlist = [100, 150]
# radiuslist = [20, 40]
# ligandlist = ['hexylamine', 'undecylamine']
# densitylist = [60]
#===============================================================================
###Replace atom names in mol2 file.
for thick, radius, curve, density, ligand, ratio in itertools.product(thicklist, radiuslist, curvelist, densitylist, ligandlist, ratiolist):
    try:
        goldname = ('{}x{}-{}curve-{}ratio{}dense{}'.format(thick, radius, curve, ratio, density, ligand))
        InputFileName = '{}{}.mol2'.format(infilepath,goldname)
        OutputFileName = '{}{}_rn.ac'.format(outfilepath,goldname)
        InputFile = open(InputFileName, 'r')
        OutputFile = open(OutputFileName, 'r')
        RenamedFileName = '{}{}_rn.mol2'.format(renamedfilepath,goldname)
        RenamedFile = open(RenamedFileName, 'w+')
        newnamelist = []
        for line in OutputFile:
            #print line[:4] == 'ATOM'
            if line[:4] == 'ATOM':
                newname = line.split()
                newnamelist.append(newname[-1])
                #print line[-3:]
        readindex = 0
        InputFile.seek(0)
        trigger = 0
        DUnum = 0
        print 'The size of the renamed list is {}.\n'.format(len(newnamelist))
        for line in InputFile:
            if line[:13] == '@<TRIPOS>ATOM':
                trigger = 1
                RenamedFile.write(line)
                continue
            if line[:13] == '@<TRIPOS>BOND':
                trigger = 0
            if trigger == 1:
                if newnamelist[readindex].find('DU') == 1:
                    newline = line[:55]+'hc'+line[58:]
                    DUnum = DUnum + 1
                else:
                    if newnamelist[readindex] == 's':
                        newline = line[:55]+'ss'+line[58:]
                    else:
                        newline = line[:55]+newnamelist[readindex]+line[58:]
                RenamedFile.write(newline)
                readindex = readindex + 1
            else:
                RenamedFile.write(line)
            #Leave off all info related to substructure,formal charges, etc.
            #if line[:21] == '@<TRIPOS>SUBSTRUCTURE': break
        print '{} atoms reassigned as from DU to hc.'.format(DUnum)
        RenamedFile.close()
        OutputFile.close()
        InputFile.close()
    except IOError:
        print 'File not found.'
        continue
        
