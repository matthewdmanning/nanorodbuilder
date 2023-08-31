'''
Created on Mar 22, 2017

@author: mdmannin
'''

### Diffusion within a binary substitutional solution.
### This file provides functions for simulating diffusion in a binary simple cubic system.
### The jump frequency between AA, AB, and BB atoms have been weighted.
###
### Two modes are available:
### diffuser() will cause every atom to switch places during one time step (provided that all other atoms around it have not already moved in the same time step.
### This is simple, but does not take into account the differences in jump frequencies between the two atom types.
### weighteddiffuser() randomly selected a percentage of atoms in the crystal to jump and randomly chooses the type of jump (AA, AB, or BB) based on the jump frequency.
### These two functions return positions for all atoms at all time steps in the simulations.
### The stats() function then returns the mean-square-displacements of type A and type B atoms.

import itertools
import random

import matplotlib.pylab as plt
import numpy as np

np.set_printoptions(precision=6, threshold=100, edgeitems=6, suppress=True)


def main():
    # omega = 1.2*10**4
    # T = 900
    # Set conditions to save data and check that filepaths are valid before beginning simulation. Will throw error if path does not exist.
    savedata = True
    plotdata = True
    weighted = True
    dpath = '/home/mdmannin/Dropbox/'
    msdfile = open('{}msd.txt'.format(dpath), 'w+')
    msdfile.close()
    posfile = open('{}pos.txt'.format(dpath), 'w+')
    posfile.close()

    # Set timesteps, relative jump frequency, dimensions of supercell, number of atoms, and ratio of A:B,
    tlength = 1000
    freq = np.array([2, 3, 5])
    dimen = np.array([6, 6, 6])
    numatom = np.prod(dimen)
    ratio = np.array([0.5, 0.5])
    numAB = numatom * ratio

    # Parameters for weighted diffusion by atom type.
    moveratio = 0.25
    movenum = int(numatom * moveratio)
    weights = [5, 7]
    weighting = [weights[0] / (weights[0] + weights[1]), weights[1] / (weights[0] + weights[1])]

    # Iniitialize coordinates in supercell as list.
    coords = list(itertools.product(range(dimen[0]), range(dimen[1]), range(dimen[2])))

    # Initialize array for staying positions of each atom in supercell at each time point. Indices: [atom number, time step, ijk]
    loca = np.empty([numatom, tlength, 3], dtype='int')  # (atomid, timestep, ijk)
    # Randomize atom location and store as first time point.
    current = coords[:]
    random.shuffle(current)
    loca[:, 0, :] = np.array(current)

    # print loca[:,0,:]
    ### Data Structures
    # coords: 2-D array listing all points in lattice and their coordinates
    # current: 1-D array containing tuples of ijk coords at atoms with index = atom number
    # loca: 3-D array with axis=0 -> atom number, axis=1 -> timestep, axis=3 -> [i,j,k]
    # movelist: list of atom positions which have been moved

    def boundary(co):
        tracker = [0, 0, 0]
        new = list(co)
        for i in range(3):
            if co[i] == dimen[i]:
                new[i] = 0
                tracker[i] = 1
            #		print 'Dimension {} is {}.'.format(i,co[i])
            elif co[i] < 0:
                new[i] = dimen[i] - 1
                tracker[i] = -1
        #		print 'Dimension {} is less than 0. Resetting to {}.'.format(i, dimen[i]-1)
        # if np.any(tracker == 1): print 'New coords: {}.'.format(tracker)
        # print new
        return tracker, tuple(new)

    def jump(center, neighbors):
        # 1)Construct probabilities for jump to each.
        # 2) Random jump. Add atom coords to movelist and update positions in loca and current
        # Inputs as atom types not coordinates.
        # Output is index(of neighbors) of atom chosen for replacement
        nprob = np.zeros(len(neighbors))
        ncum = []
        for i, neighbor in enumerate(neighbors):
            if center == 0 and neighbor == 0:
                nprob[i] = nprob[i - 1] + freq[0]
            elif center == 1 and neighbor == 1:
                nprob[i] = nprob[i - 1] + freq[2]
            else:
                nprob[i] = nprob[i - 1] + freq[1]
            ncum.append(np.sum(nprob))
        draw = np.random.random() * np.sum(nprob)
        for i in range(len(neighbors)):
            if draw - ncum[i] < 0.001:
                return i

    def stats(positions):
        # Input matrix of atom positions at each time point.
        # Calculates and returns the mean-sum-of-squares for both atom types, separately.
        inipos = positions[:, 0, :]
        dev = np.empty([numAB[0], 2])
        ss = np.empty(tlength)
        msd = np.empty([tlength, 2])
        for t in range(tlength):
            finalpos = positions[:, t, :]
            move = np.subtract(finalpos, inipos)
            deviate = np.sum(np.square(move), axis=1)
            dev[:, 0] = deviate[:numAB[0]]
            dev[:, 1] = deviate[numAB[0]:]
            for atype in [0, 1]:
                ss[atype] = np.sum(dev[:, atype], axis=0)
                msd[t, atype] = ss[atype] / numAB[atype]
        # print msd.shape
        print
        msd[:, 0]
        print
        msd[:, 1]

        return msd

    def plotter(msd, weighted=False):
        # Plot msd data and save figure.
        weigh = ''
        if weighted == True: weigh = '-weighted'
        fig = plt.figure()
        plt.plot(np.arange(tlength), msd[:, 0], 'x', label='Type A')
        plt.plot(np.arange(tlength), msd[:, 1], '*', label='Type B')
        plt.title('RMS of A and B Atoms in a Binary Crystal')
        plt.savefig('{}msd{}'.format(dpath, weigh))
        return

    def saver(positions, msd, weighted=False):
        # Write data to files.
        weigh = ''
        if weighted == True: weigh = '-weighted'
        msdfile = open('{}msd-typea{}.txt'.format(dpath, weigh), 'w+')
        np.savetxt(msdfile, msd[:, 0], fmt=('%1.5f'))
        msdfile.close()
        msdfile = open('{}msd-typeb{}.txt'.format(dpath, weigh), 'w+')
        np.savetxt(msdfile, msd[:, 1], fmt=('%1.5f'))
        msdfile.close()
        # posfile = open('{}pos.txt'.format(dpath),'w+')
        # print positions[1]
        # with posfile:
        # for slice_2d in positions:
        #	np.savetxt(posfile,positions.astype(np.int16))
        # posfile.close()
        return

    def diffuser():
        # Function conductions diffusion simulation using global variables and returns to positions of atoms at each time point.
        ### Workflow
        # 0) Copy current into last. Clear current and movelist.
        # 1) Loop through all atom positions. Find nearest neighbors, lookup atom numbers and atom types.
        # 2) Call Jump function and obtain move for pair of atoms.
        # 3) Store updated locations of atoms in loca matrix.
        # 4) Return loca.
        nearest = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        ### Loop through all timesteps.
        for t in range(tlength):
            # Reset matrices that store list of previously moved atoms, lists of neighbor coords and neighbor atom types
            movelist = []
            neighbors = [0, 0, 0, 0, 0, 0]
            neightype = [0, 0, 0, 0, 0, 0]
            # Loop through all atoms in supercell, skipping over atoms which have already moved during the current timestep.
            for original in coords:
                # Find atom number of selected atom and skip if it has already been moved.
                oindex = current.index(original)
                if original in movelist: continue
                # Find atom type of center atom.
                otype = 0
                if oindex > numAB[0]: otype = 1
                movelist.append(original)
                # ID coordinates of nearest neighbors and correct for periodic boundary conditions.
                tracker = [0] * 6
                # print tracker
                # print boundary(np.add(original,nearest[0]).tolist())
                for n, near in enumerate(nearest): tracker[n], neighbors[n] = boundary(np.add(original, near).tolist())
                # Get atom types of nearest neighbors
                for n, nei in enumerate(neighbors):
                    nei = tuple(nei)
                    # print boundary(tuple(nei))[1]
                    nindex = current.index(boundary(tuple(nei))[1])
                    if nindex > numAB[0]: neightype[n] = 1
                # Choose which neighbor will be selected for jump.
                chosen = jump(otype, neightype)
                # Find coordinates and atom number of chosen atom.
                jumpcoord = tuple(neighbors[chosen])
                # print jumpcoord
                jumpind = current.index(jumpcoord)
                # Update movelist, current positions of atoms, position time series, and boundary tracker.
                movelist.append(jumpind)
                current[jumpind] = original
                current[oindex] = jumpcoord
                loca[oindex, t, :] = np.add(loca[oindex, t - 1, :], np.array(nearest[chosen]))
                loca[jumpind, t, :] = np.subtract(loca[jumpind, t - 1, :], np.array(nearest[chosen]))
        return loca

    def jumpweight(type1, neightype, type2):
        allowed = []
        for i, nei in enumerate(neightype):
            if nei == type2:
                allowed.append(i)
        if len(allowed) < 1: return False, np.NAN
        jumper = int(allowed[np.random.randint(len(allowed))])
        return True, jumper

    def weighteddiffuser():
        # Function conductions diffusion simulation using global variables and returns to positions of atoms at each time point.
        ### Workflow
        # 0) Copy current into last. Clear current and movelist.
        # 1) Loop through all atom positions. Find nearest neighbors, lookup atom numbers and atom types.
        # 2) Call Jump function and obtain move for pair of atoms.
        # 3) Store updated locations of atoms in loca matrix.
        # 4) Return loca.
        nearest = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        ### Loop through all timesteps.
        for t in range(tlength):
            print
            't = {}'.format(t)
            # Reset matrices that store list of previously moved atoms, lists of neighbor coords and neighbor atom types
            movelist = []
            neighbors = [0, 0, 0, 0, 0, 0]
            neightype = [0, 0, 0, 0, 0, 0]
            loca[:, t, :] = np.array(current)
            if t == 0: continue
            # Loop through all atoms in supercell, skipping over atoms which have already moved during the current timestep.
            for i in range(movenum):
                validjump = False
                while validjump == False:
                    rnum = np.random.random()
                    if rnum < 0.2:
                        index1 = int(np.random.random() * numAB[0])
                        type1 = 0
                        type2 = 1
                    if rnum < 0.5:
                        index1 = int(np.random.random() * numAB[0])
                        type1 = 0
                        type2 = 1
                    else:
                        index1 = int(np.random.random() * numAB[0] + numAB[0])
                        type1 = 1
                        type2 = 1
                    if index1 in movelist: continue
                    tracker = [0] * 6
                    original = current[index1]
                    for n, near in enumerate(nearest): tracker[n], neighbors[n] = boundary(
                        np.add(original, near).tolist())
                    # Get atom types of nearest neighbors
                    for n, nei in enumerate(neighbors):
                        nei = tuple(nei)
                        # print boundary(tuple(nei))[1]
                        nindex = current.index(boundary(tuple(nei))[1])
                        if nindex > numAB[0]: neightype[n] = 1
                    # Choose which neighbor will be selected for jump.
                    validjump, chosen = jumpweight(type1, neightype, type2)
                    # Find coordinates and atom number of chosen atom.
                    if type(chosen) != int: continue
                    jumpcoord = tuple(neighbors[chosen])
                    # print jumpcoord
                    index2 = current.index(jumpcoord)
                    if current[index2] in movelist: validjump = False
                movelist.append(index1)
                movelist.append(index2)
                # Update movelist, current positions of atoms, position time series, and boundary tracker.
                movelist.append(jumpcoord)
                current[index2] = original
                current[index1] = jumpcoord
                loca[index1, t, :] = np.add(loca[index1, t - 1, :], np.array(nearest[chosen]))
                loca[index2, t, :] = np.subtract(loca[index2, t - 1, :], np.array(nearest[chosen]))
        return loca

    #### Outline.
    if weighted == False:
        positions = diffuser()
        msd = stats(positions)
        # print msd
        # print positions[0]
        if savedata == True: saver(positions, msd)
        if plotdata == True: plotter(msd)
    else:
        positions = weighteddiffuser()
        print
        'Weighted Diffusion'
        msd = stats(positions)
        plt.figure()
        plt.plot(range(tlength), loca[5, :, 0])
        plt.plot(range(tlength), loca[5, :, 1])
        plt.plot(range(tlength), loca[5, :, 2])
        plt.show()
        if savedata == True: saver(positions, msd, True)


# if plotdata == True: plotter(msd, True)


if __name__ == '__main__':
    main()
