#; Script to Fix mol2 file. 
#; write by Ropón-Palacios G. 
#; date: April 16, 2022. 
#; use: pymol -c fixmol2.py 

import glob
import os
from pymol import cmd 

files = glob.glob('GOL[A-Z].pdb')
for f in files:
	name = f.split(".")[0]  
	cmd.load(f)
	cmd.remove("hydro")
	cmd.save(name+"notH.pdb")
 
	#; It funtcion add Hydrogens atoms at pH 7.4 with gastaiger charges. 
	os.system("obabel {mol} -O {name} -p 7.4".format(mol=name+"notH.pdb", name=name+".mol2"))
	os.system("sed \'s/{p}/GOL/g\' {mol} > tmp1.pdb".format(p=name+"notH.pdb", mol=name+".mol2"))
	os.system("sed \'s/GOL662/GOL/g\' tmp1.pdb > {ofile}".format(ofile=name+"_fixed.mol2"))
	#os.system("mv {} tmp/".format(f))
	os.system("mv {} tmp/".format(name+".mol2"))
	os.system("rm tmp1.pdb") 

print("[Info] Done!!")

