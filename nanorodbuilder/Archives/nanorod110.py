#python nanoscript filter
import numpy as np
import math as m
import random
import itertools

__surfthick__ = 2

def initializer():
    filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/NanorodXYZFiles/110/'
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
# def body_check (x,y,z):
#     return (abs(y) + abs(z)-(2**0.5)*apothem < 0.00001) and (abs(y) - apothem < 0.00001) and (abs(z) - apothem < 0.00001) 
# def surf_bond_check (x,y,z):
#     return ((abs(y) + abs(z)-(2**0.5)*apothem + __surfthick__ > 0.00001) or (abs(y) - apothem + __surfthick__ > 0.00001) or (abs(z) - apothem + __surfthick__ > 0.0001)) and abs(x) - bodylength/2 < 0.0001
# def core_check (x,y,z):
#     return ((abs(y) + abs(z)-(2**0.5)*apothem + __surfthick__ < 0.00001) or (abs(y) - apothem + __surfthick__ < 0.00001) or (abs(z) - apothem + __surfthick__ < 0.0001)) and abs(x) - bodylength/2 < 0.0001
# def surf_bond_normal_check(x,y,z):
#     return (z - apothem + 1.5 > 0.0001)
# def surf_bond_diagonal_check(x,y,z):
#     return (y + z - (2**0.5)*apothem + 2 > 0.00001)
# def cap_check(x,y,z):
#     slopeheight =  apothem-(x-bodylength/2)/np.tan(capangle)
#     return (x - bodylength/2 > 0.0001) and (x - bodylength/2 - caplength < 0.0001) and (abs(y)-slopeheight < 0.0001 and abs(z)-slopeheight < 0.0001)
# def cap_bond_check(x,y,z):
#     slopeheight =  apothem-(x-bodylength/2)/np.tan(capangle)
#     return (z - slopeheight > 0.0001) or (x - bodylength/2 + 2 - caplength > 0.0001 and y + z > 0.0001 and z -y > 0.0001)
# def arrayerror(array):
#     if len(array) == 0:
#         return 0
#     else:
#         return 1
#===============================================================================

def body_check (x,y,z, bodylength, apothem, caplength, capangle):
    return (abs(y) + abs(x)-(1+m.tan(m.pi/8))*apothem < 0.0001) and (abs(y) - apothem < 0.001) and (abs(x) - apothem < 0.0001) 
def core_check (x,y,z, bodylength, apothem, caplength, capangle):
    return ((abs(y) + abs(x)-(1+m.tan(m.pi/8))*apothem + __surfthick__*1.5 < 0.00001) and (abs(y) - apothem + __surfthick__ < 0.00001) and (abs(x) - apothem + __surfthick__ < 0.0001)) and abs(z)-bodylength/2+caplength < 0.0001
def surf_bond_check (x,y,z, bodylength, apothem, caplength, capangle):
    return ((abs(y) + abs(x)-(1+m.tan(m.pi/8))*apothem + __surfthick__*1.5 > 0.0001) or (abs(y) - apothem > -__surfthick__) or (abs(x) - apothem > -__surfthick__)) and abs(z)-bodylength/2+caplength < 0.0001
def surf_bond_normal_check(x,y,z, bodylength, apothem, caplength, capangle):
    return (y - apothem + __surfthick__ > 0.0001)
def surf_bond_diagonal_check(x,y,z, bodylength, apothem, caplength, capangle):
    return (y + x - (1+m.tan(m.pi/8))*apothem + __surfthick__*1.5 > 0.0001)
def cap_check(x,y,z, bodylength, apothem, caplength, capangle):
    slopeheight =  apothem - np.tan(capangle)*(abs(z)-bodylength/2+caplength)
    print (z-(bodylength-2*caplength)/2)
    print slopeheight
    if slopeheight < 2.0: return False
    #print (abs(z) - bodylength/2 + caplength > 0.0001) and (abs(y)-slopeheight < 0.0001 and abs(x)-slopeheight < 0.0001) and (abs(z) - bodylength/2 < 0.0001)#and (abs(y) - abs(x)-(1+m.tan(m.pi/8))*apothem < 0.00001) #or (abs(y) + abs(x)-(1+m.tan(m.pi/8))*apothem < 0.00001)# and (abs(y) - apothem < 0.00001) and (abs(x) - apothem  < 0.0001))
    return (abs(y)-slopeheight < 0.0001 and abs(x)-slopeheight < 0.0001) and (abs(z) - bodylength/2 < 0.0001)#and (abs(y) - abs(x)-(1+m.tan(m.pi/8))*apothem < 0.00001) #or (abs(y) + abs(x)-(1+m.tan(m.pi/8))*apothem < 0.00001)# and (abs(y) - apothem < 0.00001) and (abs(x) - apothem  < 0.0001))
def cap_bond_check(x,y,z, bodylength, apothem, caplength, capangle):
    slopeheight =  apothem - np.tan(capangle)*(abs(z)-bodylength/2+caplength)
    if slopeheight < 0.1: return False
    capthick = __surfthick__*1.5
    return (y - slopeheight + capthick > 0.0001) or (x + y - (1+m.tan(m.pi/8))*apothem + capthick > 0.0001 and x - y - (1+m.tan(m.pi/8))*apothem + capthick > 0.0001)
    #Define length and radius in Angstroms
def main():
    filepath, lengthlist, radiuslist, densitylist = initializer()
    for length, radius, density in itertools.product(lengthlist, radiuslist, densitylist):
        #Set dimensions of nanoparticle
        apothem = radius
        caplength = apothem/1.5
        bodylength = length
        capangle = m.pi/3.5
        #Establish minimum distance between ligands.
        cutoff = (100/density)**0.5
        print cutoff
        shelloutfile = "{}{}x{}nanorodshell.xyz".format(filepath, length, radius)
        coreoutfile = "{}{}x{}nanorodcore.xyz".format(filepath,length, radius)
        normalboundfile = "{}{}x{}normalboundgold".format(filepath,length, radius)
        diagonalboundfile = "{}{}x{}diagonalboundgold".format(filepath,length, radius)
        rotatedfile = "{}{}x{}rotated.xyz".format(filepath,length, radius)
        shellout = open(shelloutfile,'w+')
        coreout = open(coreoutfile,'w+')
        ### Set dimensions of FCC unit cell
        maxx, maxy, maxz = 5.756, 4.07, 2.88
        unitlist = np.array([[0,0,0],[2.88, 0, 0], [1.439, 2.035, 1.439], [4.317, 2.035, 1.439]])
        unitlist = unitlist + [2.88, 2.035, 0]
        #unitlist = np.array([[0,0,0],[2.88, 0, 0], [1.439, 2.035, 1.439], [4.317, 2.035, 1.439]])
        #Set how many times to extend lattice in the x,y,z direction;
        xextender = int((bodylength)/maxx)+2
        yextender = 4*int(2*apothem/maxy)+2
        zextender = 4*int(2*apothem/maxz)+2
        #Calculate facet dimensions and densities.
        sidelength = 2*apothem*m.tan(m.pi/8)
        #Initialize arrays
        corearray,shellarray,normalarray,diagonalarray,caparray,capbondarray = [],[],[],[],[],[]
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
                for xprime in range (0, xextender):
                        newx = x+(xprime-xextender/2)*maxx   
                        #print newx     
                        if body_check(newx,newy,z, length, apothem, caplength, capangle) == 1:
                            for zprime in range (0, zextender):
                                newz = z+(zprime-zextender/2)*maxz
                                if abs(newz) > (bodylength/2): continue      
                                if surf_bond_check(newx,newy,newz, length, apothem, caplength, capangle):
                                    if (newx, newy, newz) not in shellarray:
                                        shellarray.append([newx, newy, newz])
                                        shellnum += 1
                                    if surf_bond_normal_check(newx, newy, newz, length, apothem, caplength, capangle) and (newx, newy, newz) not in normalarray and (newx, newy, newz) not in diagonalarray:
                                        normalarray.append([newx, newy, newz])
                                        normalnum += 1
                                    elif surf_bond_diagonal_check(newx, newy, newz, length, apothem, caplength, capangle) and (newx, newy, newz) not in diagonalarray and (newx, newy, newz) not in normalarray:
                                        diagonalarray.append([newx, newy, newz])
                                        diagonalnum += 1
                                elif core_check(newx,newy,newz, length, apothem, caplength, capangle) and (newx, newy, newz) not in corearray:
                                    corearray.append([newx, newy, newz])
                                    corenum += 1
                                elif cap_check(newx, newy,newz, length, apothem, caplength, capangle):
                                    if cap_bond_check(newx,newy,newz, length, apothem, caplength, capangle) and (newx, newy, newz) not in normalarray:
                                        normalarray.append([newx, newy, newz])
                                        #normalarray.append([newy, -newx, newz])
                                        normalnum += 1
                                    elif (newx, newy, newz) not in corearray and (newx, newy, newz) not in shellarray:
                                        corearray.append([newx, newy, newz])
                                        #corearray.append([newy, -newx, newz])
                                        corenum += 1
    
        print 'Number of shell atoms: {}'.format(len(shellarray))
        print 'Number of core atoms: {}'.format(len(corearray))
        print 'Number of normal atoms: {}'.format(len(normalarray))
        print 'Number of diagonal atoms: {}'.format(len(diagonalarray))
        print 'Number of cap atoms: {}'.format(capnum)
        shellout.write('{}'.format(len(shellarray)))
        shellout.write(' \n \n')
        for row in shellarray:
            shellout.write('Au {:.4f} {:.4f} {:.4f} \n'.format(row[0], row[1], row[2]))
        coreout.write('{}'.format(len(corearray)))
        coreout.write(' \n \n')
        for row in corearray:
            coreout.write('Au {:.4f} {:.4f} {:.4f} \n'.format(row[0], row[1], row[2]))
        rotatedcoords = []
        boundarray = []
        for liganddensity in densitylist:
            rotatedcoords[:] = []
            ####Ligand Atom Selector
            for orient in range (0,8):
            #if orient == 0: print '\n Number of normal atoms selected for binding: '
            #if orient == 1: print '\n Number of diagonal atoms selected for binding: ' 
            #Loop through to generate four facets for each family or orientations
                #print 'Coords: {}. '.format(len(rotatedcoords))
                #print len(rotatedcoords)
                if orient > 3:
                    filename = '{}{}-{}.xyz'.format(normalboundfile,liganddensity,int(orient/2))
                    arraycopy = np.copy(normalarray)
                    facetligandnum = int((sidelength*bodylength+2*caplength*sidelength/(3*m.sin(capangle)))*liganddensity/100)
                    print 'Expect {} atoms. Actual: '.format(facetligandnum)
                if orient < 4:
                    filename = '{}{}-{}.xyz'.format(diagonalboundfile,liganddensity,int((orient-1)/2))
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
                rotationmatrix = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
                transformmatrix = np.array([[0,1,0],[-1,0,0],[0,0,1]])
                for i in range(orient): rotationmatrix = np.matmul(rotationmatrix,transformmatrix)
                print rotationmatrix
                #print rotationmatrix
                xyzfile = open(filename, 'w+')
                xyzfile.write('{}'.format(len(boundarray)))
                #print len(boundarray)
                xyzfile.write(' \n \n')
                #print '{}.'.format(len(boundarray))
                for atom in range(len(boundarray)):
                    #print np.hstack(boundarray[atom])
                    #print np.matmul(np.hstack(boundarray[atom]),rotationmatrix)
                    newcoords = np.ndarray.tolist(np.matmul(np.hstack(boundarray[atom]),rotationmatrix))
                    #print newcoords
                    #print newcoords[0]
                    rotatedcoords.append([newcoords[0],newcoords[1],newcoords[2]])
                    xyzfile.write('Au {:.4f} {:.4f} {:.4f} \n'.format(boundarray[atom][0],boundarray[atom][1],boundarray[atom][2]))
                xyzfile.close()
            rotateout = open(rotatedfile,'w+')
            rotateout.write('{}'.format(len(rotatedcoords)))
            rotateout.write(' \n \n')
            for atom in rotatedcoords:
                rotateout.write('Au {:.4f} {:.4f} {:.4f} \n'.format(atom[0],atom[1],atom[2]))
            rotateout.close()
            
if __name__ == "__main__":
    main()
