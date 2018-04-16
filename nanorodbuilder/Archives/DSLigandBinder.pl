use MdmDiscoveryScript;
use DSCommands;
use NucleicAcidDiscoveryScript;
use ProteinDiscoveryScript;
use Math::Trig;
#Extract parameters from batch file
my $inputfilename = '/home/mdmannin/git/nanorodbuilder/batch_input.txt';
open (my $fh, "<", $inputfilename)
		or die 'Input file not readable!\n';
		while (my $row = <$fh>) {
			my @newrow;
			chomp $row;
			if (index($row, 'LENGTH') != -1) {
				@newrow = split('=', $row);
				our @lengthlist = split(',', $newrow[1]);
			} elsif (index($row, 'RADIUS') != -1) {
				@newrow = split('=', $row);
				our @radiuslist = split(',', $newrow[1]);
			} elsif (index($row, 'RATIO') != -1) {
				@newrow = split('=', $row);
				our @ratiolist = split(',', $newrow[1]);
			} elsif (index($row, 'DENSITY') != -1) {
				@newrow = split('=', $row);
				our @densitylist = split(',', $newrow[1]);
			} elsif (index($row, 'LIGAND1') != -1) {
				@newrow = split('=', $row);
				our @ligand1list = split(',', $newrow[1]);
			} elsif (index($row, 'LIGAND2') != -1) {
				@newrow = split('=', $row);
				our @ligand2list = split(',', $newrow[1]);
			} elsif (index($row, 'SULFINDEX1') != -1) {
				@newrow = split('=', $row);
				our @sulf1indices = split(',', $newrow[1]);
			} elsif (index($row, 'SULFINDEX2') != -1) {
				@newrow = split('=', $row);
				our @sulf2indices = split(',', $newrow[1]);
			}
		}
#Declare variables for customization	
	my $ligandratio = $ratiolist[0];	#Sets ratio of Type 1: Type 2 ligands.
    my $ligandfolder = "/home/mdmannin/Desktop/Nanoparticles/Ligands/";
	my $XYZpath = "/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/NanorodXYZFiles/";
    my $cleanstep = 500;
#Define 45 degree rotation for diagonal faces
my $eighthrotate = Mdm::Matrix4x4::Create(
			cos(pi/4), 0, sin(pi/4), 0,
			0, 1, 0, 0,
			-sin(pi/4), 0, cos(pi/4), 0,
			0,	0,	0,	1
			);
#Placeholder for rotation matrix. 
my $rotationMatrix = Mdm::Matrix4x4::Create(
			1, 0, 0, 0,
			0, 0, 1, 0,
			0, 1, 0, 0,
			0,	0,	0,	1
			);
my $docrotate = Mdm::Matrix4x4::Create(
			cos(pi/2), 0, -sin(pi/2), 0,
			0, 1, 0, 0,
			sin(pi/2), 0, cos(pi/2), 0,
			0,	0,	0,	1
			);
my $AuSdist = 1.3;
for my $radius (@radiuslist) {
	for my $length (@lengthlist) {
		for my $ligandindex (0..$#ligand1list) {
		    for my $density (@densitylist) {
				my $ligand1name = $ligand1list[$ligandindex];
				my $ligand2name = $ligand2list[$ligandindex];
				my $sulfIndex1 = $sulf1indices[$ligandindex];
				my $sulfIndex2 = $sulf2indices[$ligandindex];
				my $ligandpath1 = $ligandfolder . $ligand1name . ".mol2";
				my $ligandpath2 = $ligandfolder . $ligand2name . ".mol2";
				my $rationame = int($ligandratio*100);
				#Set path names for gold core and shell atoms
				my $shellname = $XYZpath . $length . "x" . $radius . "nanorodshell.xyz";
				my $corename = $XYZpath . $length . "x" . $radius . "nanorodcore.xyz";
				my $savefile = "/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/dsmol2/" . $length . "x" . $radius . "-" . $rationame . 'ratio' . $density . 'dense' . $ligand1name . ".mol2";
				print $savefile . "/n";
				if (-f $savefile) {
					print("File of same name already exists. Moving to next build.\n");
					next;
				}
				print("Building nanorod. Length: " . $length . ". Radius: " . $radius . ". Charged/Uncharged Ratio: " . $ligandratio . ". Ligand Density: " . $density . ". With " . $ligand1name . " and " . $ligand2name . " ligands.\n"); 
				
				###Initialize arrays and create Document,Molecule,and Gold Chain/Residue.
				my @gold_surf;
				my @gold_core;
				my @goldatoms;
				my @ligands;
				my $document = Mdm::Document::Create();
				my $firstmolecule = $document->CreateMolecule();
				$firstmolecule->Name('nanorod' . $length . "x" . $radius . $ligand1name);
				my $chain = $firstmolecule->CreateChain();
				my $residue = $chain->CreateResidue(Mdm::unknownResidueType);
				$residue->Name('gold');			
				
				for (my $i=0; $i <= 7; $i++) {
					my $facetfile;
					if ($i < 4)
					{
						$facetfile = $XYZpath . $length . "x" . $radius . "diagonalboundgold" . $density . '-' . $i . ".xyz";
					}
					else
					{
						my $j = $i - 4;
						$facetfile = $XYZpath . $length . "x" . $radius . "normalboundgold" . $density . '-' . $j . ".xyz";
					}
					print $facetfile;
					FileInsert($document, $facetfile);    
					my $goldsurf = ();
					my $newfacet = $document->Molecules->Item(1);
					my $newatoms = $document->Molecules->Item(1)->Atoms;
					$document->DeselectAll();
					for $atom (@$newatoms) {
						my $copy = $firstmolecule->CreateAtom(Mdm::Au);
						my $XYZ = $atom->XYZ;
						$copy->XYZ($XYZ);
						$document->SetParent($residue,$copy);
						$copy->Select();
					}
					my $goldsurf = $document->SelectedObjects;
					$document->DeselectAll();
					$newfacet->Select();
					$document->DeleteObjects();
					#Facet rotation: Rotates entire document and inserts a new facet from file
					#$newfacet->NcsMatrix('eighthrotate') = $eighthrotate;
					#$document->SelectAll();
					#$document->TransformSelection($eighthrotate);
					$document->DeselectAll();
					if ($i<4) {
						#Facet rotation: Rotates entire document and inserts a new facet from file
						$firstmolecule->NcsMatrix("eighthrotate") = $eighthrotate;
						$firstmolecule->Select();
						$document->ApplyTransformationMatrix($firstmolecule, "eighthrotate", CreateCopy=>0);
					}
					#my $deletecount = 0;
					#for my $atom (@$goldsurf) {
						#$document->DeselectAll();
						#$atom->Select();
						#$document->SelectByRadius(0.2, "Atom");
						#$atom->Deselect();
						#if ($document->SelectedObjects->Count > 0)
						#{
							#$document->DeselectAll();
							#$atom->Select();
							#$document->DeleteObjects;
							#$deletecount = $deletecount + 1;
						#}	
					#}
					#print('Number of ligand bound atoms deleted: ' . $deletecount);
				for my $endAtom (@$goldsurf) {
						$document->DeselectAll();
						push(@goldatoms, $endatom);
						my $ligandtype = rand(1);
						my $rotatedligand;
						my $ligandname;
						my $sulfurA;
					#Insert ligand
						if ($ligandtype < $ligandratio) {
							FileInsert($document, $ligandpath1);
							$sulfurA=$document->Molecules->Item(1)->Atoms->Item($sulfIndex1);
							$rotatedligand = $document->Molecules->Item(1);
							$ligandname = $ligand1name;
							$rotatedligand->Name('phi');
						}
						else	{
							FileInsert($document, $ligandpath2);
							$sulfurA=$document->Molecules->Item(1)->Atoms->Item($sulfIndex2);
							$rotatedligand = $document->Molecules->Item(1);
							$ligandname = $ligand2name;
							$rotatedligand->Name('pho');
						}
						#Copy ligand
						$document->DeselectAll();
						$rotatedligand->Select();
						
						#####Fix messed up rotations.
						#if ($i>3) {
						#$rotatedligand->NcsMatrix("rotation") = $rotationMatrix;	
						#$document->ApplyTransformationMatrix($rotatedligand, 'rotation',CreateCopy=>False);
						#}
						$document->DeselectAll();
						#Store coordinates of sulfur and terminal atom in ligand molecule.
						my $xyz = $sulfurA->XYZ;
							my $SX=$xyz->X;
							my $SY=$xyz->Y;
							my $SZ=$xyz->Z;
						#Store position of gold atom.
						my $xyz=$endAtom->XYZ;
							my $endX=$xyz->X;
							my $endY=$xyz->Y;
							my $endZ=$xyz->Z;
						#Move ligand molecule so that sulfur is at gold atom
						$document->DeselectAll();
						$rotatedligand->Select();
						$document->TranslateSelection( { X => $AuSdist*$endX-$SX, Y => $AuSdist*$endY-$SY, Z => $AuSdist*$endZ-$SZ } );
						if (abs($endY) > $length*0.4) {
							$document->TranslateSelection( { X => 0, Y => $endY*0.2, Z => 0 } );
						}
						$document->DeselectAll();
						#Bind Au and S atoms and place ligand into parent molecule.
						$document->SetParent( $firstmolecule, $rotatedligand );
						$document->CreateBond($endAtom, $sulfurA);				
					}
					#my $firstmolecule = $document->Molecules->Item(0);
					if ($i>3) {
						#Facet rotation: Rotates entire document and inserts a new facet from file
						$firstmolecule->NcsMatrix('docrotate') = $docrotate;
						$document->SelectAll();
						$document->TransformSelection($docrotate);
						}
					if ($i<4) {
						#Facet rotation: Rotates entire document and inserts a new facet from file
						$firstmolecule->NcsMatrix("eighthrotate") = $eighthrotate;
						$firstmolecule->Select();
						$document->ApplyTransformationMatrix($firstmolecule, "eighthrotate", CreateCopy=>0);	
					}
				}
					#my $ligands = $document->Chains;
					#pop @ligands;
					##$document->CalculateCharges();
					#for my $ligand (@$ligands) {
						#$document->DeselectAll();
						#$ligand->Select();
						#$document->Clean($cleanstep);
					#}
					FileInsert($document, $shellname);				
					FileInsert($document, $corename);
					my $goldcore = $document->Molecules->Item(2);
					my $coreatoms = $goldcore->Atoms;
					my $goldsurface = $document->Molecules->Item(1);
					my $surfatoms = $goldsurface->Atoms;
					for $atom (@$coreatoms) {
						my $copy = $firstmolecule->CreateAtom(Mdm::Au);
						my $XYZ = $atom->XYZ;
						$copy->XYZ($XYZ);
						$document->SetParent($residue,$copy);
					}
					for $atom (@$surfatoms) {
						my $copy = $firstmolecule->CreateAtom(Mdm::Au);
						my $XYZ = $atom->XYZ;
						$copy->XYZ($XYZ);
						$document->SetParent($residue,$copy);
					}
					$document->DeselectAll();
					$goldcore->Select();
					$goldsurface->Select();
					$document->DeleteObjects();
					#for my $atom (@$originals) {
						#$document->DeselectAll();
						#$atom->Select();
						#$document->SelectByRadius(0.2, "Atom");
						#$atom->Deselect();
						#$deletecount = $deletecount + $document->SelectedObjects->Count;
						#$document->DeleteObjects;
					#}
					print('Total number of atoms deleted: ' . $deletecount);
					### Reclean ligands to prevent closecontact with surface gold.
					#$document->SelectLigands;
					#$document->Clean();
					#for my $ligand (@$ligands) {
						#$document->DeselectAll();
						#$ligand->Select();
						#$document->Clean();
					#}
				$document->Save( $savefile, "mol2" );
				$document->Close();
				print("Saving file:" . $savefile . ".\n");
				sleep 5;
			}
		}
	}
}
