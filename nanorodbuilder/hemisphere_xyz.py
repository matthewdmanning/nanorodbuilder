'''
Created on Jun 12, 2017

@author: mdmannin
'''
import itertools

import numpy as np

frcmod_file = "h4_8a"
if __name__ == '__main__':
    radius = 29.
    box_size = 70. - 5.
    # Set dimensions of nanoparticle
    caplength = radius / 2.
    # ## Set dimensions of FCC unit cell
    unit_vector = 4.
    bond_length = 4 * np.sqrt(2)
    unit_array = np.array(
        [[0., 0., 0.], [unit_vector / 2., unit_vector / 2., 0.], [unit_vector / 2., 0., unit_vector / 2.],
         [0., unit_vector / 2., unit_vector / 2.]])
    cellsize = 2 * int(radius / unit_vector) + 2
    print(cellsize)
    #cells = range(-cellsize, cellsize)
    cells = range(0,cellsize)
    sphere_coords_array, solvent_coords_array = [], []
    for x_extend, y_extend, z_extend in itertools.product(cells, cells, cells):
        translate = np.array([x_extend * unit_vector, y_extend * unit_vector, z_extend * unit_vector])
        for coords in unit_array:
            new_coords = coords + translate
            dist = np.sqrt(np.sum(np.multiply(new_coords, new_coords)))
            if dist  - radius <  0.0001 and radius - dist < 5.000: #and radius - dist - 10. < 0.001:
                sphere_coords_array.append(new_coords.tolist())
    #print(sphere_coords_array)
    filename = '/home/mdmannin/Desktop/H4_Histone/h4_tail_leap.in'
    leapshell = open(filename, 'w')
    #leapshell.write('logfile {}.log\n'.format('Au_sphere'))
    leapshell.write('verbosity 0\n')
    leapshell.write('source /home/mdmannin/amber16/dat/leap/cmd/oldff/leaprc.ff10\n')
    leapshell.write('source leaprc.water.opc\n')
    #leapshell.write('source leaprc.protein.ff14SB\n')
    leapshell.write('set default PBradii bondi\n')
    leapshell.write('WAT=OPC\nHOH=OPC\n')    
    leapshell.write('mods = loadamberparams frcmod.phosaa10\n')
    leapshell.write('mods = loadamberparams {}.frcmod\n'.format(frcmod_file))
    # ## Create gold residue and atoms and move to coordinates in list.
    #leapshell.write('set default nocenter on\n')
    leapshell.write('sphere = createUnit sph\n')
    leapshell.write('sp_res = createResidue sphres\n')
    leapshell.write('add sphere sp_res\n')
    for atom_index, atom_coords in enumerate(sphere_coords_array):
        #print row
        leapshell.write('z{0} = createAtom l{0} ZZ 0\n'.format(atom_index))
        #leapshell.write('set z{} element ZZ\n'.format(atom_index))
        leapshell.write('set z{0} position {{ {1[0]} {1[1]} {1[2]} }}\n'.format(atom_index, atom_coords))
        leapshell.write('add sphere.1 z{}\n'.format(atom_index))
        #if abs(atom_coords[0]) - 1.0 < 0.0001 and abs(atom_coords[1]) - 1.0 < 0.0001 and radius - abs(atom_coords[2]) > 1.0:
        if radius - np.sqrt(atom_coords[0]**2 + atom_coords[1]**2  + atom_coords[2]**2) < 1.0 and abs(atom_coords[0] - atom_coords[1]) - 2.0 < 0.0001 and abs(atom_coords[1] - atom_coords[2]) - 2.0 < 0.0001:
                terminal_num = atom_index + 1
                terminal_coords = atom_coords[:]
                print('Linker atom coords: {}'.format(terminal_coords))
    #leapshell.write('tail = loadpdb original-last-linker.pdb\n')
    leapshell.write('tail = loadpdbusingseq original_original-last-new.pdb {NSER GLY ARG GLY LYS GLY GLY LYS GLY LEU GLY LYS GLY GLY ALA LYS ARG HIP ARG LYS VAL LEU ARG ASP ASN ILE}\n')
    #leapshell.write('transform tail { { 0 0 1 } { 0 1 0 } { 1 0 0 }}\n')
    #leapshell.write('alignaxes tail\n')
    #leapshell.write('translate tail { 0 0 50 }\n')
    leapshell.write('link = createunit l\n')
    leapshell.write('linkres = createresidue lr\n')
    leapshell.write('add link linkres\n')
    leapshell.write('clink = createAtom CL1 CT 0\n')
    leapshell.write('set clink element C\n')
    leapshell.write('set clink position { 0 0 0}\n')
    leapshell.write('add link.1 clink\n')
    leapshell.write('h1 = createatom h1 HC 0\n')
    leapshell.write('translate h1 {1 0 0}\n')
    leapshell.write('h2 = createatom h2 HC 0\n')
    leapshell.write('translate h2 { -1 0 0}\n')
    leapshell.write('add link.1 h1\n')
    leapshell.write('add link.1 h2\n')
    leapshell.write('bond link.1.CL1 link.1.h1\n')
    leapshell.write('bond link.1.CL1 link.1.h2\n')
    #leapshell.write('set clink position { 0 0 30 }\n')
    #leapshell.write('translate tail { 0 0 50 }\n')
    carbon_coords = radius / np.sqrt(3) + 1.5
    leapshell.write('translate link.1 {{ {0} {0} {0} }}\n'.format(carbon_coords))
    leapshell.write('combined = combine { sphere link tail }\n')
    leapshell.write('bond combined.2.CL1 combined.3.N\n')
    leapshell.write('bond combined.1.N combined.2.CL1\n')
    leapshell.write('bondbydistance combined.1 {}\n'.format(bond_length))
    leapshell.write('check combined\n')
    leapshell.write('savemol2 combined au_sphere_h4_8a.mol2 0\n')
    # corner_list = [[box_size, 0, 0], [0, box_size, 0], [0, 0, box_size], [box_size, box_size, 0], [box_size, 0, box_size], [box_size, box_size, box_size]]
    # leapshell.write('gridmol = createresidue grid\n')
    # leapshell.write('add combined gridmol\n')
    # for index, coords in enumerate(corner_list):
    #    leapshell.write('g{0} = createatom ga{0} AA 0\n'.format(index))
    #    leapshell.write('set g{0} position {{ {1[0]} {1[1]} {1[2]} }}\n'.format(index, coords))
    #    leapshell.write('add gridmol g{0}\n'.format(index))
    ### Create unit with dummy atoms for create water box of specific size.
    #===========================================================================
    # leapshell.write('gridunit = createunit grid_unit\n ')
    # leapshell.write('gridmol = createresidue grid_mol\n')
    # leapshell.write('add gridunit gridmol\n')
    # unit_list = unit_array.tolist()
    # grid_size = int(box_size / unit_vector)
    # last_index = 0
    # for index, x, y, z, unit_coords in enumerate(itertools.product(range(grid_size), range(grid_size), range(grid_size), unit_list)):
    #     leapshell.write('g{0} = createatom D D 0\n'.format(index))
    #     atom_array = np.array([x * unit_vector, y * unit_vector, z * unit_vector])  + np.array(unit_coords)
    #     atom_coords = atom_array.tolist()
    #     leapshell.write('set g{0} positions {{ {1[0]} {1[1]} {1[2]} }}\n'.format(index, atom_coords))
    #     last_index = index
    #===========================================================================
    #leapshell.write('grid_combined = combine { combined gridunit }\n')
    #leapshell.write('solvateBox grid_combined OPCBOX ')
    #for index in range(last_index): leapshell.write('remove grid_combined g{0}\n'.format(index))
    leapshell.close()