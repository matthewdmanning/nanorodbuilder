#!/usr/bin/perl -w
use MdmDiscoveryScript;
use DSCommands;
use NucleicAcidDiscoveryScript;
use ProteinDiscoveryScript;
use Math::Trig;
opendir(DIR, '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/BondedMol2/');
my @file_array = grep(/\.mol2$/, readdir(DIR));
closedir DIR;
#print(@file_array);
for my $file_name (@file_array) {
    print($file_name . "\n");
    $document = Mdm::Document::Create();
    FileInsert($document, $file_name);
    #Only if you want "ligands" relaxed; other you can use: SelectAll($mydocument)
    #my $ligands = GetLigands($mydocument);
    #Select($document, GetLigands($document));
    #$ligands->Select();
    my $core = $document->Molecules->Item(0);
    #$document->Select($core);
    $core->Select();
    $document->InvertSelection();
    $document->Clean(Iterations        => 500, IncludeIntermolecularForces => True, CleanOnlySelectedAtoms => True,
        IgnoreAutoCalculateChargesFlag => True);
    #$document->Clean();
    $file_name =~ s/".mol2"/"_rn.mol2"/;
    $document->Save($file_name, "mol2");
    $document->Close();
}