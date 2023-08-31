import subprocess

import nanorodbatch

nanorodbatch.nanorodbuilder()
# subprocess.call('~/.bashrc',shell=True, executable='/bin/bash')
subprocess.call(['~/BIOVIA/DiscoveryStudio2016/bin/./perl.sh', '~/git/nanorodbuilder/DSLigandBinder.pl'], shell=True)
