#!/bin/bash
# Script for moving files

echo "Building nanorods"
prinf nanorod_batch.txt
echo "Creating XYZ files"
python<nanorodbatch.py
echo "Attaching Ligands in Discovery Studio"
dsrun<DSLigandBinder.pl
echo "Converting atom name to be Amber-compatible."
python<atomtypereplacer.py
./renamingscript.sh
python<renamescript.py
echo "Preparing tLeap scripts for binding gold."
python<leapgoldbinder.py
echo "Executing tLeap scripts."
./~/Desktop/Nanoparticles/RoughNanorods/LeapScripts/masterscript.sh
