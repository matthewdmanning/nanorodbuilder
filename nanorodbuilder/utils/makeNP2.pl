#!/usr/bin/perl -w
use MdmDiscoveryScript;
use DSCommands;
use NucleicAcidDiscoveryScript;
use Math::Trig;

#####--makeNP.pl--
#	Author: Jessica Nash (janash@ncsu.edu)
#           NC State University
#           Yingling Group
#           Dec. 2014
#			Modified for two ligand types - March 2015
#
#           This script uses Mdm Discovery Script (part of Discovery Studio Visualizer) to make
#           nanoparticles which can be input into leap
#
#
#
#####



#Declare variables
my @gs;            #initialization of empty array
my @thetaS;        #initialization of empty array    
my @rS;            #initialization of empty array
my $type1 = 36;        #Number of ligand_type 1 ligands to attach to surface.
my $last=60;       #total number of ligands to attach to surface
my @bond1;         #initialization of empty array for bonding of Au-S (first atom)    
my @bond2;         #initialization of empty array for bonding of Au-S (second atom)
my $plRad=3.0;         #distance between bonded gold atom and surface atom
my $npRad=8;       #Radius of nanoparticle in angstrom
my $sulfIndex=0; 	#Index of sulfur atom on ligand
my $termIndex=19;	#Index of terminal atom on ligand

#Set paths for gold cube and for ligand
our $filename =  '/home/mdmannin/Desktop/NanorodScripts/Nanoparticle_Scripts_from_Jessica/cube_bonded2.mol2';    #Path for gold cube pdb
our $fpath = '/home/jessica/Documents/NPs/PEG3/PEG3_c_red.mol2';              #Path for ligand ligand_type 1 mol2

our $fpath2 = '/home/jessica/Documents/NPs/PEG3/PEG3_un_red.mol2';              #path for ligand ligand_type 2 mol2



##------------Create spherical gold core from rectangular cube
    
#Open gold pdb file - create object gold_cube
my $gold_cube = DiscoveryScript::Open( $filename );
#Create object gold_atoms
my $gold_atoms = $gold_cube->Atoms;
#Create molecule object for gold cube
my $molecule = $gold_cube->Molecules->Item(0);
$gold_cube->DeselectAll();
#Create O at 0,0,0
$gold_cube->InsertFromSmiles('[O]');
my $my_oxygen=$gold_cube->Molecules->Item(1);

#Find middle atom
$my_oxygen->Select();
$gold_cube->SelectByRadius(1,'Atom');
my $middle_atom=$gold_cube->SelectedObjects->Atoms->Item(0);

#Make spherical nanoparticle
$gold_cube->DeselectAll();
$my_oxygen->Select();
my $new_selects=$gold_cube->SelectByRadius($npRad, 'Atom');
$gold_cube->InvertSelection();
$gold_cube->DeleteObjects();
#Delete oxygen atom used to find center
$my_oxygen->Select();
$gold_cube->DeleteObjects();


##Find and color surface atoms

#Set up color information for surface atoms
my $color_S = Mdm::Color::Create();
$color_S->Red   = 255;
$color_S->Green = 211;
$color_S->Blue  = 211;
my $count=0;
our @surfAtoms;

#Set up color information for surface atoms
my $color_B = Mdm::Color::Create();
$color_S->Red   = 0;
$color_S->Green = 211;
$color_S->Blue  = 255;

#Loop through each atom in gold core to find surface atoms
foreach my $atom (@$gold_atoms)
{
    my $xyz = $atom->XYZ;
    my $my_X=$xyz->X;
    my $my_Y=$xyz->Y;
    my $my_Z=$xyz->Z;
    $gold_cube->DeselectAll();
    $atom->Select();    
    my $new_selects=$gold_cube->SelectByRadius(4.0, 'Atom');
    my $neighbor_molecule=$gold_cube->SelectedObjects->Atoms;
    
    my $mol_search=@$neighbor_molecule;

    $gold_cube->DeselectAll(); 

    #If atom has less than 12 neighbors by distance, it is a surface atom
     if ($mol_search<12)
     {
     my $delBond=$atom->BondCount;
         #Delete bond if gold atom already has 8 bonds (max number in leap)
         if ($delBond>8)
         {
            my $bondes=$atom->Bond(0);
            $gold_cube->DeselectAll();
            $bondes->Select();
            $gold_cube->DeleteObjects();
         }
     
     $atom->Select();
     $gold_cube->SetColor($color_S);
     $gold_cube->SetAtomDisplayStyle( {DisplayStyle => Mdm::styleAtomBallAndStick} );
     $count=$count+1;
     $gold_cube->DeleteObjects();
     push(@surfAtoms,$atom);
     } 
}

for my $i ( 1 ..$last )
{
#initialize random index variable
my $random_index;  
if ($i<($type1+1)){
#insert ligand
FileInsert($gold_cube, $fpath);
} else {
    print "Inserting ligand_type 2";
	#insert ligand
	FileInsert($gold_cube, $fpath2);
	}


#Calculate translation needed to put sulfur at 0,0,0
my $ligandMolecule=$gold_cube->Molecules->Item($i);
$gold_cube->DeselectAll();
$ligandMolecule->Select();
my $ligAtoms=$ligandMolecule->Atoms;
my $COM_1 = $ligandMolecule->CenterOfMass;
  my $COM_X=$COM_1->X;
  my $COM_Y=$COM_1->Y;
  my $COM_Z=$COM_1->Z;

my $sulfurA=$gold_cube->Molecules->Item($i)->Atoms->Item($sulfIndex);
my $xyz = $sulfurA->XYZ;
 my $my_SX=$xyz->X;
 my $my_SY=$xyz->Y;
 my $my_SZ=$xyz->Z;
 push(@bond1,$sulfurA);
my $terminalA=$gold_cube->Molecules->Item($i)->Atoms->Item($termIndex);
my $xyz=$terminalA->XYZ;
 my $my_TX=$xyz->X;
 my $my_TY=$xyz->Y;
 my $my_TZ=$xyz->Z;
  
my $moveX=$COM_X-$my_SX;
my $moveY=$COM_Y-$my_SY;
my $moveZ=$COM_Z-$my_SZ;

#Move ligand molecule so that sulfur is at origin
$gold_cube->DeselectAll();
$ligandMolecule->Select();
$gold_cube->TranslateSelection( { X => $moveX-$COM_X, Y => $moveY-$COM_Y, Z => $moveZ-$COM_Z } );
my $endAtom;

$random_index = int(rand()*$count);  
$gold_cube->DeselectAll();
#Choose gold surface atom
$endAtom=$surfAtoms[$random_index];
push(@bond2,$endAtom);
#Remove chosen atom from surfAtoms so another ligand cannot be attached
splice(@surfAtoms,$random_index,1);
$count=$count-1;

#Store position of end atom
my $xyz= $endAtom->XYZ;
 my $my_EX=$xyz->X;
 my $my_EY=$xyz->Y;
 my $my_EZ=$xyz->Z;
 
#Create tethers - surface and terminal atom
$gold_cube->DeselectAll();
$endAtom->Select();
$terminalA->Select();
## created tethers from selection
my $firstTether=CreateTethers($gold_cube);

#Create tethers - sulfur and middle atom
$gold_cube->DeselectAll();
$middle_atom->Select();
$sulfurA->Select();
## created tethers from selection
my $secondTether=CreateTethers($gold_cube);

#Align
AlignMolecules($gold_cube);
push(@thetaS,0);
push(@gS,0);
push(@rS,0);

if ($my_EX!=0)
{
    $gS[$i-1]=atan($my_EY/$my_EX);
}

if (($my_EX<0) || ($my_EX<0)&&($my_EY<0))
{
    $gS[$i-1]=$gS[$i-1]+pi;
} 

    our $R=($my_EX**2+$my_EY**2+$my_EZ**2)**0.5;
    our $Th=acos($my_EZ/$R);
    $rS[$i-1]=$R;
    $thetaS[$i-1]=$Th;
    #print "The stored value is radius:$R\t theta:$Th\t phi:$gS[$i]\n";

    if ($i==$last)
    {
        $gold_cube->DeselectAll();
        # Remove the tethers
        my $monitors = $gold_cube->TupleMonitors;
            foreach my $monitor (@$monitors)
            {
                $monitor->Select;
            }
            $gold_cube->DeleteObjects;

        for my $i ( 1 ..$last )
        {
         
            my $sulfurA=$gold_cube->Molecules->Item(1)->Atoms->Item(0);
            my $xyz = $sulfurA->XYZ;
            my $my_SX=$xyz->X;
            my $my_SY=$xyz->Y;
            my $my_SZ=$xyz->Z;
            
            my $endAtom=$bond2[$i-1];
            my $xyz= $endAtom->XYZ;
            my $my_EX=$xyz->X;
            my $my_EY=$xyz->Y;
            my $my_EZ=$xyz->Z;
            
            my $distance=(($my_EX-$my_SX)**2+($my_EY-$my_SY)**2+($my_EZ-$my_SZ)**2)**0.5;
            my $pRad=$plRad+$distance;
            
           $gold_cube->DeselectAll();
           #Move ligand to surface   
           my $ligandMolecule=$gold_cube->Molecules->Item(1);
           $ligandMolecule->Select();
           my $Xmove=$pRad*sin($thetaS[$i-1])*cos($gS[$i-1]);
           my $Ymove=$pRad*sin($thetaS[$i-1])*sin($gS[$i-1]);
           my $Zmove=$pRad*cos($thetaS[$i-1]);
           $gold_cube->TranslateSelection( { X => $Xmove, Y => $Ymove, Z => $Zmove } );
            
            
            $gold_cube->ChangeObjectParent( $ligandMolecule, $gold_cube->Molecules->Item(0) );
            my $bond = $gold_cube->CreateBond($bond1[$i-1], $bond2[$i-1]);
            
            
            
          }
   
    }
}

