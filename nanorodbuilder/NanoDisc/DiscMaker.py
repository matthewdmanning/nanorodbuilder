#!/usr/bin/env python
import numpy as np
import itertools

surfthick = 2

'''
Plan
	1) Check diameter in xy plane.
		-- x**2 + y**2 < r**2
	2) Check thickness in z(x,y) or constant thickness.
		-- Constant thickness = t
		-- Lens shaped disc: Effective radius = radius of curvature - t/2 + x.
	3) Select surface atoms and separate into flat or rim.
	4) Randomly choose ligand atoms for surfaces.
'''

def initializer():
	filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/DiscXYZFiles/'
	print 'Files saving to: {}.\n'.format(filepath)
	batchfilename = 'batch_input_disc.txt'
	batchfile = open(batchfilename, 'r')
	thicklist, radiuslist, densitylist, curvelist = [], [], [], []
	for line in batchfile:
		line = line.strip()
		#print line
		if line.find("THICKNESS=") != -1:
			thickbite = line[line.index('=')+1:].split(',')
			for i in thickbite:thicklist.append(int(i))
		elif line.find("RADIUS=") != -1:
			radiusbite = line[line.index('=')+1:].split(',')
			for i in radiusbite:radiuslist.append(int(i))
		elif line.find("CURVE=") != -1:
			curvebite = line[line.index('=')+1:].split(',')
			for i in curvebite:curvelist.append(int(i))
		elif line.find("DENSITY=") != -1:
			densitybite = line[line.index('=')+1:].split(',')
			for i in densitybite:densitylist.append(int(i))
	batchfile.close()
	return filepath, thicklist, radiuslist, curvelist, densitylist

def radiuscheck (x,y, thick, radius, curve):
	return x**2 + y**2 - radius**2 <= 0.001
def effectivexyz (x,y,z, thick, curve):
	return x, y, curve - thick/2 + z
	#return curve - thick/2 + x, curve - thick/2 + y, curve - thick/2 + z
def lenscheck (x,y,z, thick, radius, curve):
	xeff, yeff, zeff = effectivexyz(x,y,z,thick,curve)
	return xeff**2 + yeff**2 + zeff**2 - curve**2 <= 0.001
def rimsurfcheck (x,y,z, thick, radius, curve):
	return x**2 + y**2 - (radius - surfthick)**2 >= 0.001
def lenssurfcheck (x,y,z, thick, radius, curve):
	xeff, yeff, zeff = effectivexyz(x,y,z,thick,curve)
	return xeff**2 + yeff**2 + zeff**2 - (curve - surfthick)**2	>= 0.001

def filewriter(filename, array):
	outfile = open(filename,'w')
	outfile.write('{} \n \n'.format(len(array)))
	for row in array: outfile.write('Au {:.4f} {:.4f} {:.4f} \n'.format(row[0], row[1], row[2]))
	outfile.close()
	
def distcheck(point1, point2, cutoff):
	return ((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 + (point1[2] - point2[2])**2)**0.5 - cutoff <= 0.001
'''
def coordconverter(array,a,b,c,system='cubic'):
	if system=='cubic':
'''	

#Construct Nanodisc
def main():
	filepath, thicklist, radiuslist, curvelist, densitylist = initializer()
	for thick, radius, curve in itertools.product(thicklist, radiuslist,curvelist):
		curve = curvelist[0]
		parameters = thick, radius, curve, densitylist, filepath
		#File names for XYZ.
		shelloutfile = "{}{}x{}nanodiscshell.xyz".format(filepath, thick, radius)
		coreoutfile = "{}{}x{}nanodisccore.xyz".format(filepath, thick, radius)
		### Set dimensions of FCC unit cell
		vertx, verty, vertz = 4.07, 4.07, 4.07
		#unitlist = np.array([[0,0,0],[2.035, 2.035, 0], [2.035, 0 , 2.035], [0, 2.035, 2.035]])
		#Set how many times to extend lattice in the x,y,z direction;
		xextender = int(radius/vertx)+2
		yextender = int(radius/verty)+2
		zextender = int(thick/vertz)+2
		#Initialize arrays
		corearray,flatarray,rimarray = [],[],[]
		#Construct coords of atoms in xy-plane
		xlist = [vertx * x for x in range(-xextender, xextender + 1)]
		ylist = [verty * y for y in range(-yextender, yextender + 1)]
		xyzlist = [(x,y,0) for x,y in itertools.product(xlist,ylist)]
		for xy in xyzlist: xyzlist.append([xy[0] + vertx/2, xy[1] + verty/2, 0])
		for xyz in xyzlist: xyzlist.append([xyz[0] + vertx/2,0, xyz[2] + vertz/2])
		zlist = [vertz * range(zextender+1)]
		for x,y,z in xyzlist:
			if radiuscheck(x,y, thick, radius, curve) == 0: continue
			for zadd in zlist:
				z = z + zadd
				if lenscheck(x,y,z, thick, radius, curve) == 0: continue
				if rimsurfcheck(x, y, z, thick, radius, curve):
					rimarray.append([x,y,z])
					rimarray.append([x,y,-z])
				elif lenssurfcheck(x,y,z, thick, radius, curve):
					flatarray.append([x,y,z])
					flatarray.append([x,y,-z])
				else:
					corearray.append([x,y,z])
					corearray.append([x,y,-z])				

		print 'Number of core atoms: {}'.format(len(corearray))
		print 'Number of flat atoms: {}'.format(len(flatarray))
		print 'Number of rim atoms: {}'.format(len(rimarray))
		filewriter(coreoutfile,np.vstack(corearray))
		filewriter(shelloutfile,np.vstack([flatarray,rimarray]))
		randomligands([rimarray,flatarray],parameters)

#Loop through all atoms in array and randomly place ligands that meet cutoff requirements.
#Write selected atoms to output file and return append listed of selected atoms.
def arrayloop(array, filename, boundarray, cutoff):
	np.random.shuffle(array)
	newbonds = []
	for atom in array:
		boundarray, array, newbonds = overlapcheck(atom,boundarray,array,cutoff,newbonds)
	filewriter(filename,newbonds)
	return boundarray

#Check atom for cutoff requirements against list of previously selected atoms.
def overlapcheck(atom,boundarray,array,cutoff,newbonds):
	for bonded in boundarray:
		if distcheck(atom,bonded,cutoff) == 0: continue
		else:
			return boundarray, array, newbonds
	boundarray.append(atom)
	newbonds.append(atom)
	return boundarray, array, newbonds

#Randomly select atoms for ligand binding from a list of arrays and write to XYZ file.
def randomligands(arraylist,parameters):
	thick, radius, curve, densitylist, filepath = parameters
	boundarray = []
	for density in densitylist:
		cutoff = (100/density)**0.5
		rimboundfile = "{}{}x{}-{}dense-rimboundgold.xyz".format(filepath, thick, radius, density)
		flatboundfile = "{}{}x{}-{}dense-flatboundgold.xyz".format(filepath, thick, radius, density)
		filelist = [rimboundfile,flatboundfile]
		####Ligand Atom Selector
		boundarray[:] = []
		for array, filename in zip(arraylist, filelist):
			boundarray = arrayloop(array, filename, boundarray, cutoff)
		
if __name__ == "__main__": main()
