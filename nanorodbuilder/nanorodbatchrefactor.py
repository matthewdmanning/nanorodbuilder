# python nanoscript filter
import collections
import itertools
import math as m
import random

import numpy as np

import geometry_tools
import nanorod_builder
from geometry_tools import nanorod_xyz_builder
from script_writing_tools import xyz_writer

__surfthick__ = float(2 * 1000)
capangle = m.pi / 4

diag_add = float(-1 * 1000)


def get_iterable(x):
    if isinstance(x, collections.Iterable):
        return x
    else:
        return (x,)


def initializer():
    filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/NanorodXYZFiles/'
    print('Files saving to: {}.\n'.format(filepath))
    batch_file_name = 'batch_input.txt'
    batch_file = open(batch_file_name, 'r')
    lengthlist, radiuslist, densitylist = [], [], []

    for line in batch_file:
        line = line.strip()
        if line.find("LENGTH=") != -1:
            lengthbite = line[line.index('=') + 1:].split(',')
            for i in lengthbite: lengthlist.append(int(i))
        elif line.find("RADIUS=") != -1:
            radiusbite = line[line.index('=') + 1:].split(',')
            for i in radiusbite: radiuslist.append(int(i))
        elif line.find("DENSITY=") != -1:
            densitybite = line[line.index('=') + 1:].split(',')
            for i in densitybite: densitylist.append(int(i))
    batch_file.close()
    return filepath, lengthlist, radiuslist, densitylist




# Define length and radius in Angstroms


def ligand_picker(in_length, in_radius, in_density, surface_coords_list, filepath, in_ratio_percent=100,
                  rotatedcoords=[]):
    # Selects ligand bound atom from a list of surface atoms.
    # Ratio allows for more even distribution of charged ligands. This was a problem for smaller ratios.
    scaled_length = float(in_length * 1000)
    # print('Length: {}'.format(scaled_length))
    scaled_radius = float(in_radius * 1000)
    # Input density is in ligand per square nanometer.
    # Needs to be converted to ligands per square (Angstroms / 1000) or (Nanometers / 10000)
    ligand_surface_density = float(in_density)  # / (10**8)
    print('Target ligand surface density: {}\n'.format(ligand_surface_density))
    charged_ratio_decimal = float(in_ratio_percent) / 100.
    charged_num, uncharged_num = nanorod_builder.calc_ligand_num(in_length, in_radius, ligand_surface_density,
                                                                 charged_ratio_decimal)
    total_ligand_number = charged_num + uncharged_num
    print('Number of ligands from calc_ligand_num: {}.\n'.format(total_ligand_number))
    ligand_separation_correction_factor = 1.03
    min_ligand_site_separation = (ligand_separation_correction_factor / (
        ligand_surface_density * charged_ratio_decimal * np.pi)) ** 0.5
    if min_ligand_site_separation - 3500. < 0.0001: min_ligand_site_separation = 3.500 * 1000
    print('Smallest cutoff distance: {} Angstroms.'.format(min_ligand_site_separation / 1000))
    ligand, dist = 0, 0

    # for n in arrayrange: surface_atom_coords_array[n] = surface_coords_list[n]
    bound_gold_coords_list = []
    long_ligand_separation = min_ligand_site_separation * 3
    short_ligand_separation = min_ligand_site_separation
    for min_ligand_site_separation in [long_ligand_separation, short_ligand_separation]:
        surface_atom_coords_array = np.empty([len(surface_coords_list), 3], dtype='int')
        np.random.shuffle(surface_atom_coords_array)
        for i, tripel in enumerate(surface_coords_list):
            surface_atom_coords_array[i, :] = np.array(tripel, dtype='int')
        while surface_atom_coords_array.shape[0] > 0 and len(
                bound_gold_coords_list) < total_ligand_number:  # ligand < facetligandnum-1 and
            separated = True
            randindex = random.randint(0, surface_atom_coords_array.shape[0] - 1)
            # Check that randomly selected gold atom is not within the min_ligand_site_separation distance
            if len(bound_gold_coords_list) == 0 and len(rotatedcoords) == 0:
                # print('First ligand chosen. Both arrays are empty.')
                bound_gold_coords_list.append(surface_atom_coords_array.tolist()[randindex])
                # print(surface_atom_coords_array.tolist()[randindex])
                # print(bound_gold_coords_list[0])
                surface_atom_coords_array = np.delete(surface_atom_coords_array, (randindex), axis=0)
                continue
            # Begin looping through previously selected ligand sites and check distance with candidate site.
            distcount = 0
            while (distcount < len(bound_gold_coords_list) - 1) and separated == True:
                dist = geometry_tools.calculate_distance(surface_atom_coords_array[randindex, :],
                                                         bound_gold_coords_list[distcount])
                # print(dist)
                if dist - min_ligand_site_separation < 0.00001:
                    surface_atom_coords_array = np.delete(surface_atom_coords_array, (randindex), axis=0)
                    separated = False
                else:
                    distcount += 1
            r = 0
            while r < len(rotatedcoords) - 1 and len(rotatedcoords) > 0 and separated == True:
                # dist = ((surface_atom_coords_array[randindex, 0] - rotatedcoords[r][0]) ** 2 + (surface_atom_coords_array[randindex, 1] - rotatedcoords[r][1]) ** 2 + (surface_atom_coords_array[randindex, 2] - rotatedcoords[r][2]) ** 2) ** 0.5
                dist = geometry_tools.calculate_distance(surface_atom_coords_array[randindex, :], rotatedcoords[r])
                if dist - min_ligand_site_separation < 0.00001:
                    surface_atom_coords_array = np.delete(surface_atom_coords_array, (randindex), axis=0)
                    separated = False
                else:
                    r += 1
            if separated == True:
                bound_gold_coords_list.append(surface_atom_coords_array.tolist()[randindex])
                surface_atom_coords_array = np.delete(surface_atom_coords_array, (randindex), axis=0)
                ligand += 1
    surface_area = 12 * scaled_radius * m.tan(m.pi / 8) * scaled_length + 2 * m.pi * scaled_radius ** 2
    actual_ligand_density = len(bound_gold_coords_list) * 10 ** 8 / surface_area
    print('Actual ligand density: {} per square nanometer.\n\n'.format(actual_ligand_density))
    return bound_gold_coords_list


def ligand_rotater(length, radius, densitylist, file_path, normalarray, diagonalarray):
    scaled_radius = radius * 1000
    scaled_length = length * 1000
    rotatedcoords = []
    boundarray = []
    surface_area = 12 * scaled_radius * m.tan(m.pi / 8) * scaled_length + 2 * m.pi * scaled_radius ** 2
    side_length = 2 * scaled_radius * m.tan(m.pi / 8)
    cap_length = scaled_radius / 2
    for ligand_density in get_iterable(densitylist):
        normalboundfile = "{}{}x{}normalboundgold".format(file_path, scaled_length, scaled_radius)
        rotatedfile = "{}{}x{}-{}dense-rotated.xyz".format(file_path, scaled_length, scaled_radius, ligand_density)
        diagonalboundfile = "{}{}x{}diagonalboundgold".format(file_path, scaled_length,
                                                              scaled_radius)  # ## Divided by constant, because theoretical density wasn't giving accurate density.
        cutoff = (100 / ligand_density) ** 0.5 / 1.3
        if cutoff - 4075. < 0.0001: cutoff = 4075.
        print('Minimum ligand separation distance is {} Angstrom.'.format(cutoff))
        rotatedcoords[:] = []
        # Loops through to generate four facets for each family or orientations
        for orient in range(0, 8):
            # if orient == 0: print '\n Number of normal atoms selected for binding: '
            # if orient == 1: print '\n Number of diagonal atoms selected for binding: '
            # print 'Coords: {}. '.format(len(rotatedcoords))
            # print len(rotatedcoords)
            if orient > 3:
                filename = '{}{}-{}.xyz'.format(normalboundfile, ligand_density, int(orient - 4))
                arraycopy = np.copy(normalarray)
                facetligandnum = int((side_length * scaled_length + 2 * cap_length * side_length / (
                    3 * m.sin(capangle))) * ligand_density / 100)
                print('Expect {} atoms. Actual: '.format(facetligandnum))
            if orient < 4:
                filename = '{}{}-{}.xyz'.format(diagonalboundfile, ligand_density, int(orient))
                arraycopy = np.copy(diagonalarray)
                facetligandnum = int(side_length * scaled_length * ligand_density / 100)
                print('Expect {} atoms. Actual: '.format(facetligandnum))

            boundarray = ligand_picker(scaled_length, scaled_radius, ligand_density, arraycopy, file_path, 1,
                                       rotatedcoords)
            # Rotated coords to correct face.
            transformmatrix = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            rotationmatrix = np.matmul([[1, 0, 0], [0, 1, 0], [0, 0, 1]], transformmatrix)
            for i in range(int(orient)): rotationmatrix = np.matmul(rotationmatrix, transformmatrix)
            xyzfile = open(filename, 'w+')
            xyzfile.write('{}'.format(len(boundarray)))
            xyzfile.write(' \n \n')
            print('{}.'.format(len(boundarray)))
            for atom in range(len(boundarray)):
                newcoords = np.ndarray.tolist(np.matmul(np.hstack(boundarray[atom]), rotationmatrix))
                rotatedcoords.append([newcoords[0], newcoords[1], newcoords[2]])
                xyzfile.write(
                    'Au {:.4f} {:.4f} {:.4f} \n'.format(boundarray[atom][0], boundarray[atom][1], boundarray[atom][2]))
            xyzfile.close()
        ligandnum = ligand_density * surface_area / 100
        random.shuffle(rotatedcoords)
        while len(rotatedcoords) > ligandnum:
            rotatedcoords.pop()
        rotateout = open(rotatedfile, 'w+')
        rotateout.write('{}'.format(len(rotatedcoords)))
        rotateout.write(' \n \n')
        for atom in rotatedcoords:
            print(atom)
            rotateout.write('Au {:.4f} {:.4f} {:.4f} \n'.format(atom[0], atom[1], atom[2]))
        rotateout.close()
        print('Density = {}. {} ligands on {} A^2'.format(len(rotatedcoords) * 100 / surface_area, len(rotatedcoords),
                                                          surface_area))


def main():
    filepath, lengthlist, radiuslist, densitylist = initializer()
    for length, radius in itertools.product(lengthlist, radiuslist):
        # Establish minimum distance between ligands.
        shell_out_name = "{}{}x{}nanorodshell.xyz".format(filepath, length, radius)
        core_out_name = "{}{}x{}nanorodcore.xyz".format(filepath, length, radius)
        cap_out_name = "{}{}x{}nanorodcap.xyz".format(filepath, length, radius)

        shellarray, corearray, normalarray, diagonalarray, caparray = nanorod_xyz_builder(length, radius)
        xyz_writer(shell_out_name, shellarray)
        xyz_writer(core_out_name, corearray)
        xyz_writer(cap_out_name, caparray)
        # for density in densitylist:
        # ligand_rotater(length, radius, density, filepath, normalarray, diagonalarray)
        # bound_array = ligand_picker(length, radius, density, shellarray, filepath, ratio=1, rotatedcoords=[])


if __name__ == "__main__":
    main()
