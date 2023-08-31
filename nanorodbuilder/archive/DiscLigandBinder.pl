use MdmDiscoveryScript;
use DSCommands;
use NucleicAcidDiscoveryScript;
use ProteinDiscoveryScript;
use Math::Trig;
#Extract parameters from batch file
my $inputfilename = '/home/mdmannin/git/nanorodbuilder/batch_input_disc.txt';
open (my $fh, "<", $inputfilename)
		or die 'Input file not readable!\n';
		while (my $row = <$fh>) {
			my @newrow;
			chomp $row;
			if (index($row, 'THICK') != -1) {
				@newrow = split('=', $row);
				our @thicklist = split(',', $newrow[1]);
			} elsif (index($row, 'RADIUS') != -1) {
				@newrow = split('=', $row);
				our @radiuslist = split(',', $newrow[1]);
			} elsif (index($row, 'CURVE') != -1) {
				@newrow = split('=', $row);
				our @curvelist = split(',', $newrow[1]);
			} elsif (index($row, 'DENSITY') != -1) {
				@newrow = split('=', $row);
				our @densitylist = split(',', $newrow[1]);
			} elsif (index($row, 'LIGAND1') != -1) {
				@newrow = split('=', $row);
				our @ligand1list = split(',', $newrow[1]);
			} elsif (index($row, 'RATIO') != -1) {
				@newrow = split('=', $row);
				our @ratiolist = split(',', $newrow[1]);
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
		
sub PrintMatrix($)
{
    my ($matrix) = @_;
    my $size = 3;
    $size = 4 if ( $matrix->ClassName eq 'Mdm::Matrix4x4' );

    for ( my $row = 0 ; $row < $size ; ++$row )
    {
        for ( my $column = 0 ; $column < $size ; ++$column )
        {
            printf "    %8.3f", $matrix->GetEntry( $row, $column );
        }
        print "\n";
    }
}

#Declare variables for customization	
    my $ligandfolder = "/home/mdmannin/Desktop/Nanoparticles/Ligands/";
    my $XYZpath = "/home/mdmannin/Desktop/Nanoparticles/Nanodiscs/DiscXYZFiles/";
    my $cleanstep = 500;
#Define matrix for rotating ligand from [001] to [100] direction.
my $rimrotate = Mdm::Matrix4x4::Create(
			0, 0, 1, 0,
			0, 1, 0, 0,
			1, 0, 0, 0,
		    0, 0, 0, 1
			);
#Placeholder for rotation matrix. 
my $rotationMatrix = Mdm::Matrix4x4::Create(
			1, 0, 0, 0,
			0, 0, 1, 0,
			0, 1, 0, 0,
			0,	0,	0,	1
			);
my $invert = Mdm::Matrix4x4::Create(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, -1, 0,
			0,	0,	0,	1
			);
my $AuSdist = 1.3;

for my $radius (@radiuslist) {
	for my $thick (@thicklist) {
		for my $curve (@curvelist) {
			for my $ligandindex (0..$#ligand1list) {
				for my $ligandratio (@ratiolist) {
					for my $density (@densitylist) {
	#Declare variables from parameters.
		my $ligand1name = $ligand1list[$ligandindex];
		my $ligand2name = $ligand2list[$ligandindex];
		my $sulfIndex1 = $sulf1indices[$ligandindex];
		my $sulfIndex2 = $sulf2indices[$ligandindex];
		my $ligandpath1 = $ligandfolder . $ligand1name . ".mol2";
		my $ligandpath2 = $ligandfolder . $ligand2name . ".mol2";
		my $rationame = int($ligandratio*100);
	#Set path names for gold core and shell atoms
		my $shellname = $XYZpath . $thick . "x" . $radius . "nanodiscshell.xyz";
		my $corename = $XYZpath . $thick . "x" . $radius . "nanodisccore.xyz";
		my $savefile = "/home/mdmannin/Desktop/Nanoparticles/Nanodiscs/discmol2/" . $thick . "x" . $radius . "-" . $curve . "curve-" . $rationame . 'ratio' . $density . 'dense' . $ligand1name . ".mol2";
		if (-f $savefile) {
			print("File of same name already exists. Moving to next build.\n");
			next;
		}
		print("Building nanorod. Thickness: " . $thick . ". Radius: " . $radius . ". Radius of curvature: " . $curve . ". Charged/Uncharged Ratio: " . $ligandratio . ". Ligand Density: " . $density . ". With " . $ligand1name . " and " . $ligand2name . " ligands.\n"); 
		
	###Initialize arrays and create Document,Molecule,and Gold Chain/Residue.
		my @gold_surf;
		my @gold_core;
		my @goldatoms;
		my @ligands;
		my $document = Mdm::Document::Create();
		my $firstmolecule = $document->CreateMolecule();
		$firstmolecule->Name('nanorod' . $thick . "x" . $radius . $ligand1name);
		my $chain = $firstmolecule->CreateChain();
		my $residue = $chain->CreateResidue(Mdm::unknownResidueType);
		$residue->Name('gold');
	###Insert rim atoms for bonding.
		my $facetfile = $XYZpath . $thick . "x" . $radius . "-" . $curve . "curve-" . $density . "dense-rimboundgold.xyz";
		print $facetfile;
		FileInsert($document, $facetfile);    
		my $goldsurf = ();
		my $newfacet = $document->Molecules->Item(1);
		my $newatoms = $document->Molecules->Item(1)->Atoms;
		$document->DeselectAll();
	### Copies all Gold atoms into single residue.
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
			#Rotate ligand so that it is normal to gold surface. Ligands files should be oriented in the [001] direction.
				$document->DeselectAll();
				$rotatedligand->Select();
				$document->TransformSelection($rimrotate);
			### Create matrix to rotate around the z-axis.
				my $zaxis = Mdm::Point::Create( 0, 0, 1 );
				my $ligand3x3 = Mdm::Matrix3x3::Create();
				if ($endX == 0) { $endX = 0.001; }
				$ligand3x3->SetRotation( $zaxis, rad2deg(atan2($endY, $endX)) );
				my $ligandmatrix = Mdm::Matrix4x4::Create();
				$ligandmatrix->RotationMatrix($ligand3x3);
				$document->TransformSelection($ligandmatrix);
			my $xyz = $sulfurA->XYZ;
				my $SX=$xyz->X;
				my $SY=$xyz->Y;
				my $SZ=$xyz->Z;
			#Move ligand molecule so that sulfur is at gold atom
				$document->TranslateSelection( { X => $AuSdist*$endX-$SX, Y => $AuSdist*$endY-$SY, Z => $AuSdist*$endZ-$SZ } );
				$document->DeselectAll();
			#Bind Au and S atoms and place ligand into parent molecule.
				$document->SetParent( $firstmolecule, $rotatedligand );
				$document->CreateBond($endAtom, $sulfurA);				
		}						
	###Insert flat/lens atoms for bonding.
		my $facetfile = $XYZpath . $thick . "x" . $radius . "-" . $curve . "curve-" . $density . "dense-flatboundgold.xyz";
		FileInsert($document, $facetfile);    
		my $goldsurf = ();
		my $newfacet = $document->Molecules->Item(1);
		my $newatoms = $document->Molecules->Item(1)->Atoms;
		$document->DeselectAll();
	### Copies all Gold atoms into single residue.
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
			#Flip ligand if binding to gold atom with a negative Z coordinate.
				$document->DeselectAll();
				$rotatedligand->Select();
			### Create matrix to rotate around the z-axis.
				my $xaxis = Mdm::Point::Create( 1, 0, 0 );
				my $x3x3 = Mdm::Matrix3x3::Create();
				$x3x3->SetRotation( $xaxis, -rad2deg(atan2($endY, abs($endZ)-$thick/2+$curve)) );
				my $xaxismatrix = Mdm::Matrix4x4::Create();
				$xaxismatrix->RotationMatrix($x3x3);
				$document->TransformSelection($xaxismatrix);
			### Create matrix to rotate around the z-axis.
				my $yaxis = Mdm::Point::Create( 0, 1, 0 );
				my $y3x3 = Mdm::Matrix3x3::Create();
				$y3x3->SetRotation( $yaxis, rad2deg(atan2($endX, abs($endZ)-$thick/2+$curve)) );
				my $yaxismatrix = Mdm::Matrix4x4::Create();
				$yaxismatrix->RotationMatrix($y3x3);
				$document->TransformSelection($yaxismatrix);
				if ($endZ < 0) { $document->TransformSelection($invert); }				
			#Move ligand molecule so that sulfur is at gold atom
				my $xyz = $sulfurA->XYZ;
				my $SX=$xyz->X;
				my $SY=$xyz->Y;
				my $SZ=$xyz->Z;
				$document->TranslateSelection( { X => $AuSdist*$endX-$SX, Y => $AuSdist*$endY-$SY, Z => $AuSdist*$endZ-$SZ } );
				$document->DeselectAll();
			#Bind Au and S atoms and place ligand into parent molecule.
				$document->SetParent( $firstmolecule, $rotatedligand );
				$document->CreateBond($endAtom, $sulfurA);				
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
			#for my $ligand (@$ligands) {
				#$document->DeselectAll();
				#$ligand->Select();
				#$document->Clean();
			#}
			$document->Save( $savefile, "mol2" );

			print "Saving file:" . $savefile . ".\n";
			sleep 5;
					}
				}
			}
		}
	}
}
