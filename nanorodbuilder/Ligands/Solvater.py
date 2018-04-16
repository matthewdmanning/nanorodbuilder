import glob
import itertools


def main(densitylist, leappath):

	def builder(ligandfile, leapfile, ligandname, density):
		#leapfile = open(leapfilename,'w')
		leapfile.write('ligand = loadMol2 {}\n'.format(ligandfile))
		if density == 0:
			leapfile.write('solvateBox ligand TIP3PBOX {25 25 15}\n')		
			leapfile.write('addions ligand Na+ 2\n')
			leapfile.write('addions ligand Cl- 0\n')
			leapfile.write('saveamberparm ligand {}-{}.prmtop {}-{}.rst7\n'.format(ligandname,density,ligandname,density))
			#leapfile.write('clearVariables\n')
			return
		elif density > 0:
			translate = (100/density)**0.5/2
			copylist = []
			for i in range(0,5):
				for j in range(0,5):
					copyname = 'ligand{}'.format(i*5+j)
					copylist.append(copyname)
					leapfile.write('{} = copy ligand\n'.format(copyname))
					leapfile.write('translate {} {{ {} {} 0 }}\n'.format(copyname,i*translate,j*translate))
			leapfile.write('combined = combine {{ {} }}\n'.format(' '.join(copylist)))
			leapfile.write('solvateBox combined TIP3PBOX 30\n')
			leapfile.write('addions ligand Na+ 32\n')
			leapfile.write('addions ligand Cl- 0\n')
			leapfile.write('saveamberparm ligand {}-{}.prmtop {}-{}.rst7\n'.format(ligandname,density,ligandname,density))
			#leapfile.write('clearVariables\n')
			return
			
	gpupath='/home/mdmannin/mnt2/Ligands/'
	infilepath='/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/LeapScripts/'
	#gpupath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Leap/'
	#infilepath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Renamed/'
	ligandfilelist = glob.glob('{}*.mol2'.format(infilepath))
	leapfile = open('{}tleap.in'.format(leappath),'w')
	leapfile.write('source /home/mdmannin/amber16/dat/leap/cmd/oldff/leaprc.ff14SB \n')
	leapfile.write('loadamberparams /home/mdmannin/amber16/dat/leap/parm/frcmod.ionsjc_tip3p \nloadamberparams gaff.dat\n')
	print(ligandfilelist)
	for ligandfile, density in itertools.product(ligandfilelist,densitylist):
		ligandname = ligandfile.replace(infilepath,'').replace('_rn.mol2','')
        # leapfilename = '{}{}-{}.in'.format(leappath, ligand_name ,density)
		#print leapfilename
		for density in densitylist:
			builder(ligandfile, leapfile, ligandname, density)

if __name__ == "__main__":
    main()
