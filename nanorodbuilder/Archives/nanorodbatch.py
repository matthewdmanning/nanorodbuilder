#python nanoscript filter
import numpy as np
import math as m
import random
import itertools
from scipy.spatial import distance

surfthick = 2.5

def initializer():
	filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/NanorodXYZFiles/'
	print 'Files saving to: {}.\n'.format(filepath)
	batchfilename = 'batch_input.txt'
	batchfile = open(batchfilename, 'r')
	lengthlist, radiuslist, densitylist = [], [], []
	
	for line in batchfile:
		line = line.strip()
		print line
		if line.find("LENGTH=") != -1:
			lengthbite = line[line.index('=')+1:].split(',')
			for i in lengthbite:lengthlist.append(int(i))
		elif line.find("RADIUS=") != -1:
			radiusbite = line[line.index('=')+1:].split(',')
			for i in radiusbite:radiuslist.append(int(i))
		elif line.find("DENSITY=") != -1:
			densitybite = line[line.index('=')+1:].split(',')
			for i in densitybite:densitylist.append(int(i))
	batchfile.close()
	return filepath, lengthlist, radiuslist, densitylist

#===============================================================================
# def bodycheck (x,y,z):
# 	return (abs(y) + abs(z)-(2**0.5)*apothem < 0.00001) and (abs(y) - apothem < 0.00001) and (abs(z) - apothem < 0.00001) 
# def surfbondcheck (x,y,z):
# 	return ((abs(y) + abs(z)-(2**0.5)*apothem + surfthick > 0.00001) or (abs(y) - apothem + surfthick > 0.00001) or (abs(z) - apothem + surfthick > 0.0001)) and abs(x) - bodylength/2 < 0.0001
# def corecheck (x,y,z):
# 	return ((abs(y) + abs(z)-(2**0.5)*apothem + surfthick < 0.00001) or (abs(y) - apothem + surfthick < 0.00001) or (abs(z) - apothem + surfthick < 0.0001)) and abs(x) - bodylength/2 < 0.0001
# def surfbondnormalcheck(x,y,z):
# 	return (z - apothem + 1.5 > 0.0001)
# def surfbonddiagonalcheck(x,y,z):
# 	return (y + z - (2**0.5)*apothem + 2 > 0.00001)
# def capcheck(x,y,z):
# 	slopeheight =  apothem-(x-bodylength/2)/np.tan(capangle)
# 	return (x - bodylength/2 > 0.0001) and (x - bodylength/2 - caplength < 0.0001) and (abs(y)-slopeheight < 0.0001 and abs(z)-slopeheight < 0.0001)
# def capbondcheck(x,y,z):
# 	slopeheight =  apothem-(x-bodylength/2)/np.tan(capangle)
# 	return (z - slopeheight > 0.0001) or (x - bodylength/2 + 2 - caplength > 0.0001 and y + z > 0.0001 and z -y > 0.0001)
# def arrayerror(array):
# 	if len(array) == 0:
# 		return 0
# 	else:
# 		return 1
#===============================================================================

def bodycheck (x,y,z, bodylength, apothem, caplength, capangle):
	return (abs(y) + abs(z)-(1+m.tan(m.pi/8))*apothem < 0.0001) and (abs(y) - apothem < 1) and (abs(z) - apothem < 1) 
def corecheck (x,y,z, bodylength, apothem, caplength, capangle):
	return ((abs(y) + abs(z)-(1+m.tan(m.pi/8))*apothem + surfthick < 0.00001) and (abs(y) - apothem + surfthick < 0.00001) and (abs(z) - apothem + surfthick < 0.0001)) and abs(x) - (bodylength)/2 + caplength < 0.0001
def surfbondcheck (x,y,z, bodylength, apothem, caplength, capangle):
	return ((abs(y) + abs(z)-(1+m.tan(m.pi/8))*apothem + surfthick > 0.0001) or (abs(y) - apothem > -surfthick) or (abs(z) - apothem > -surfthick)) and abs(x) - (bodylength)/2 + caplength - surfthick < 0.0001
def surfbondnormalcheck(x,y,z, bodylength, apothem, caplength, capangle):
	return (z - apothem + surfthick > 0.0001)
def surfbonddiagonalcheck(x,y,z, bodylength, apothem, caplength, capangle):
	return (y + z - (1+m.tan(m.pi/8))*apothem + surfthick > 0.0001)
def capcheck(x,y,z, bodylength, apothem, caplength, capangle):
	slopeheight =  apothem - np.tan(capangle)*(x-bodylength/2+caplength)
	return (x - bodylength/2 + caplength > 0.0001) and (abs(y)-slopeheight < 0.0001 and abs(z)-slopeheight < 0.0001) and ((abs(y) - abs(z)-(1+m.tan(m.pi/8))*apothem < 0.00001) or (abs(y) + abs(z)-(1+m.tan(m.pi/8))*apothem < 0.00001) and (abs(y) - apothem < 0.00001) and (abs(z) - apothem  < 0.0001))
def capbondcheck(x,y,z, bodylength, apothem, caplength, capangle):
	slopeheight =  apothem - np.tan(capangle)*(x-bodylength/2+caplength)
	capthick = surfthick*1.5
	return (z - slopeheight + capthick > 0.0001) or (y + z - (1+m.tan(m.pi/8))*apothem + capthick > 0.0001 and z - y - (1+m.tan(m.pi/8))*apothem + capthick > 0.0001)
	#Define length and radius in Angstroms

def xyzbuilder(length, radius, density):
	#Set dimensions of nanoparticle
	apothem = radius
	caplength = apothem/2
	capangle = m.pi/4
	### Set dimensions of FCC unit cell
	maxy, maxx, maxz = 4.07, 4.07, 4.07
	unitlist = np.array([[0,0,0],[2.035, 2.035, 0], [2.035, 0 , 2.035], [0, 2.035, 2.035]])
	#Set how many times to extend lattice in the x,y,z direction;
	xextender = int((length)/maxx)+2
	yextender = 2*int(2*apothem/maxy)+2
	zextender = 2*int(2*apothem/maxz)+2
	#Calculate facet dimensions and densities.
	sidelength = 2*apothem*m.tan(m.pi/8)
	#Initialize arrays
	corearray,shellarray,normalarray,diagonalarray = [],[],[],[]
	if (xextender % 2 == 1):
		xextender += 1
	#Initialize counters.
	normalnum = 0
	diagonalnum = 0
	corenum = 0
	shellnum = 0
	capnum = 0	
	for x,y,z in unitlist:
		for yprime in range (0, yextender):
			newy = y+(yprime-yextender/2)*maxy
			for zprime in range (0, zextender):
					newz = z+(zprime-zextender/2)*maxz
					if bodycheck(x,newy,newz, length, apothem, caplength, capangle) == 1:
						for i in range (0, xextender):
							newx = x+(i-xextender/2)*maxx		
							if surfbondcheck(newx,newy,newz, length, apothem, caplength, capangle):
								if (newx, newy, newz) not in shellarray:
									shellarray.append([newy, newx, newz])
									shellnum += 1
								if surfbondnormalcheck(newx, newy, newz, length, apothem, caplength, capangle) and (newx, newy, newz) not in normalarray and (newx, newy, newz) not in diagonalarray:
									normalarray.append([newy, newx, newz])
									normalnum += 1
								elif surfbonddiagonalcheck(newx, newy, newz, length, apothem, caplength, capangle) and (newx, newy, newz) not in diagonalarray and (newx, newy, newz) not in normalarray:
									diagonalarray.append([newy, newx, newz])
									diagonalnum += 1
							elif corecheck(newx,newy,newz, length, apothem, caplength, capangle) and (newx, newy, newz) not in corearray:
								corearray.append([newy, newx, newz])
								corenum += 1
							elif capcheck(newx, newy,newz, length, apothem, caplength, capangle):
								if capbondcheck(newx,newy,newz, length, apothem, caplength, capangle) and (newx, newy, newz) not in normalarray:
									normalarray.append([newy, newx, newz])
									normalarray.append([newy, -newx, newz])
									normalnum += 2
								elif (newx, newy, newz) not in corearray and (newx, newy, newz) not in shellarray:
									corearray.append([newy, newx, newz])
									corearray.append([newy, -newx, newz])
									corenum += 2

	print 'Number of shell atoms: {}'.format(len(shellarray))
	print 'Number of core atoms: {}'.format(len(corearray))
	print 'Number of normal atoms: {}'.format(len(normalarray))
	print 'Number of diagonal atoms: {}'.format(len(diagonalarray))
	print 'Number of cap atoms: {}'.format(capnum)
	return shellarray, corearray

def xyzwriter():
	shellout.write('{}'.format(len(shellarray)))
	shellout.write(' \n \n')
	for row in shellarray:
		shellout.write('Au {:.4f} {:.4f} {:.4f} \n'.format(row[0], row[1], row[2]))
	coreout.write('{}'.format(len(corearray)))
	coreout.write(' \n \n')
	for row in corearray:
		coreout.write('Au {:.4f} {:.4f} {:.4f} \n'.format(row[0], row[1], row[2]))
	return

def ligandpicker():
	rotatedcoords = []
	boundarray = []
	
	surfacearea = 12*apothem*m.tan(m.pi/8)*length + 2 * m.pi * apothem**2
	for liganddensity in densitylist:
		### Divided by constant, because theoretical density wasn't giving accurate density.
		cutoff = (100/density)**0.5/1.3
		print 'Minimum ligand separation distance is {} Angstrom.'.format(cutoff)
		rotatedcoords[:] = []
		####Ligand Atom Selector
		for orient in range (0,8):
		#if orient == 0: print '\n Number of normal atoms selected for binding: '
		#if orient == 1: print '\n Number of diagonal atoms selected for binding: ' 
		#Loop through to generate four facets for each family or orientations
			#print 'Coords: {}. '.format(len(rotatedcoords))
			#print len(rotatedcoords)
			if orient > 3:
				filename = '{}{}-{}.xyz'.format(normalboundfile,liganddensity,int(orient-4))
				arraycopy = np.copy(normalarray)
				facetligandnum = int((sidelength*bodylength+2*caplength*sidelength/(3*m.sin(capangle)))*liganddensity/100)
				print 'Expect {} atoms. Actual: '.format(facetligandnum)
			if orient < 4:
				filename = '{}{}-{}.xyz'.format(diagonalboundfile,liganddensity,int(orient))
				arraycopy = np.copy(diagonalarray)
				facetligandnum = int(sidelength*bodylength*liganddensity/100)
				print 'Expect {} atoms. Actual: '.format(facetligandnum)
			#for p in arraycopy: print
			#Loop until all ligand sites are selected
			ligand=0
			np.random.shuffle(arraycopy)
			dist = 0
			boundarray[:] = []
			while arraycopy.shape[0] > 0: #ligand < facetligandnum-1 and
				separated = True
				randindex = random.randint(0,arraycopy.shape[0]-1)
				#Check that randomly selected gold atom is not within the cutoff distance
				distcount = 0
				if len(boundarray) == 0:
					boundarray.append(np.ndarray.tolist(arraycopy[randindex]))						
					arraycopy = np.delete(arraycopy, (randindex), axis=0)
					continue
				while (distcount < len(boundarray) - 1) and separated == True:
					dist = ((arraycopy[randindex,0]-boundarray[distcount][0])**2 + (arraycopy[randindex,1]-boundarray[distcount][1])**2 + (arraycopy[randindex,2]-boundarray[distcount][2])**2)**0.5
					if dist - cutoff < 0.00001:
						arraycopy = np.delete(arraycopy, (randindex), axis=0)
						separated = False
					else :
						distcount += 1
				r = 0
				while r < len(rotatedcoords) - 1 and len(rotatedcoords) > 0 and separated == True:
					dist = ((arraycopy[randindex,0]-rotatedcoords[r][0])**2 + (arraycopy[randindex,1]-rotatedcoords[r][1])**2 + (arraycopy[randindex,2]-rotatedcoords[r][2])**2)**0.5
					if dist - cutoff < 0.00001:
						arraycopy = np.delete(arraycopy, (randindex), axis=0)
						separated = False
					else :
						r += 1
				if separated == True:
					boundarray.append(np.ndarray.tolist(arraycopy[randindex]))						
					arraycopy = np.delete(arraycopy, (randindex), axis=0)
					ligand += 1
			transformmatrix = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
			rotationmatrix = np.matmul([[1,0,0],[0,1,0],[0,0,1]],transformmatrix)
			for i in range(int(orient)): rotationmatrix = np.matmul(rotationmatrix,transformmatrix)
			#print rotationmatrix
			#print rotationmatrix
			xyzfile = open(filename, 'w+')
			xyzfile.write('{}'.format(len(boundarray)))
			#print len(boundarray)
			xyzfile.write(' \n \n')
			print '{}.'.format(len(boundarray))
			for atom in range(len(boundarray)):
				#print np.hstack(boundarray[atom])
				#print np.matmul(np.hstack(boundarray[atom]),rotationmatrix)
				newcoords = np.ndarray.tolist(np.matmul(np.hstack(boundarray[atom]),rotationmatrix))
				#print newcoords
				#print newcoords[0]
				rotatedcoords.append([newcoords[0],newcoords[1],newcoords[2]])
				xyzfile.write('Au {:.4f} {:.4f} {:.4f} \n'.format(boundarray[atom][0],boundarray[atom][1],boundarray[atom][2]))
			xyzfile.close()
		ligandnum = density*surfacearea/100
		random.shuffle(rotatedcoords)
		while len(rotatedcoords) > ligandnum:
			rotatedcoords.pop()
		rotateout = open(rotatedfile,'w+')
		rotateout.write('{}'.format(len(rotatedcoords)))
		rotateout.write(' \n \n')
		for atom in rotatedcoords:
			rotateout.write('Au {:.4f} {:.4f} {:.4f} \n'.format(atom[0],atom[1],atom[2]))
		rotateout.close()
		print('Density = {}. {} ligands on {} ^2'.format(len(rotatedcoords)*100/surfacearea, len(rotatedcoords), surfacearea))
		#pdf = []
		#for atom1, atom2 in itertools.combinations(rotatedcoords, 2):
		#	pdf.append(distance.euclidean(atom1,atom2))

def main():
	filepath, lengthlist, radiuslist, densitylist = initializer()
	for length, radius, density in itertools.product(lengthlist, radiuslist, densitylist):

		#Establish minimum distance between ligands.
		shelloutfile = "{}{}x{}nanorodshell.xyz".format(filepath, length, radius)
		coreoutfile = "{}{}x{}nanorodcore.xyz".format(filepath,length, radius)
		normalboundfile = "{}{}x{}normalboundgold".format(filepath,length, radius)
		diagonalboundfile = "{}{}x{}diagonalboundgold".format(filepath,length, radius)
		rotatedfile = "{}{}x{}-{}dense-rotated.xyz".format(filepath,length, radius,density)
		shellout = open(shelloutfile,'w+')
		coreout = open(coreoutfile,'w+')



	return shellarray, corearray
if __name__ == "__main__":
	main()
