;# Script to move GOL to geometric center of protein.
;# write by Ropón-Palacios G. 
;# date: April 16, 2022. 
;# use:  - Open tkconsole into VMD and type : "source moveGOL.tcl". 
;#       - Into Terminal type: "vmd -dispdev text -e moveGOL.tcl"                                     

;# Fitting Template with model. 
set molid1 [mol new tetr_chain_temp.pdb type pdb waitfor all]
set molid2 [mol new aqp_notgol.pdb type pdb waitfor all]

set fitprot [atomselect $molid1 "protein and name CA"]
set fitgol  [atomselect $molid1 "resname GOL"]
set ref     [atomselect $molid2 "protein and name CA"]

set matrix1 [measure fit $fitprot $ref] 
$fitgol  move $matrix1   ; # Use matrix rot+trans of protein to move GOL molecules. 

set selgol [atomselect $molid1 "resname GOL"]
$selgol writepdb GOLs_from_Template.pdb 
mol delete all 

;# Merge mols. 
set id1 [mol new aqp_notgol.pdb  type pdb waitfor all]
set id2 [mol new GOLs_from_Template.pdb type pdb waitfor all]

package require topotools
set midlist [list $id1 $id2]

set mol [::TopoTools::mergemols $midlist]
animate write pdb aqp_withGOL.pdb  $mol

;# Place GOL to your center of mass at 20 Angstrom into Z+.
mol delete all 
mol new aqp_withGOL.pdb type pdb waitfor all 

set seltxt1 [atomselect top "resname GOL"]
set chainsHETA [lsort -unique [$seltxt1 get chain]]

foreach chain $chainsHETA { 
	set sel1 [atomselect top "chain $chain"]
	$sel1 set segid GOL$chain
	$sel1 set segname GOL$chain 

	set center1 [measure center $sel1]
	set movez [expr 20 - [lindex $center1 2]]  ; # Placed a ~13 Angstrom from Selective Filter. 

	set movelist [list 0 0 $movez]
	$sel1 moveby $movelist
         
	$sel1 writepdb GOL$chain.pdb 

}

set all [atomselect top all]
$all writepdb  aqp_withGOL_fixed.pdb
puts "fix file manually for center GOL into Pore!"
quit 
