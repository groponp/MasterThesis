;#  This Script create a Citosol Enviroment for membrane+protein.  
;# written by : Ropón-Palacios G.  
;# date : April 15, 2022.
;# use: vmd -dispdev text -e CitosolPacking.tcl 

;# Setting control parameters. 
;# ===========================
set solute 1                          ; # If it's 1 On and 0 Off. 
set ionconc 0.150                     ; # 0.150 is ionic concentration in mol/L. 
set initprot "aqp_withGOL_fixed.pdb"  ; # Initial PDB file containing prot+mem orient. 
set ofile  "aqp_pvivax"               ; # Prefix of ouputs files. 
set hmr 1                             ; # 1 => yes and 0 => no 
set oS  "MacOS"                       ; # MacOS or Linux
set soluteTopology "gol.rtf"          ; # Name of soloute topology. 
set soluteTopologyG "gol_g.rtf"  


;# Script for generate topology. 
;# ============================= 
package require psfgen 
topology ../../toppar/top_all36_prot.rtf           ; # Topology update for Proteins. 
topology ../../toppar/top_all36_na_nbfix.rtf       ; # Topology for Nucleic acids.
topology ../../toppar/toppar_water_ions_nbfix.str  ; # Topology for Waters (TIP3) and Ions.  
topology ../../toppar/top_all27_new_prot_na.inp    ; # Topology for Protein.  
topology ../../toppar/top_all36_cgenff.rtf         ; # General Topology for small molecules. 
topology ../../toppar/top_all36_lipid.rtf          ; # Lipids 

if {$solute == 1} {
topology $soluteTopology        ; # Glycerol Topology generated from CHARMM-GUI suite. 
topology $soluteTopologyG
	
} else {
puts "Not solute Topology given."

}

;# Protein selection.    
mol new $initprot type pdb waitfor all
set sel [atomselect top "protein"]
set segids [lsort -unique [$sel get segid]]
foreach segid $segids {
    puts "Adding protein seig $segid to psfgen"
    pdbalias residue HIS HSD 
    pdbalias atom ILE CD1 CD
    set seg ${segid} 
    set sel [atomselect top "protein and segid $segid"]
    $sel set segid $seg 
    $sel writepdb tmp.pdb 
    segment $seg {pdb tmp.pdb}
    if {$segid == "PROA"} {patch GLUP PROA:28}    ; # The resids 28 are protonated to pH 7.0. 
    if {$segid == "PROB"} {patch GLUP PROB:28}  
    if {$segid == "PROC"} {patch GLUP PROC:28}  
    if {$segid == "PROD"} {patch GLUP PROD:28}    
    coordpdb tmp.pdb 
    guesscoord                                    ; # Add H and missing atoms. 
}  

;# Membrane selection.  
set sel [atomselect top "resname POPE"]
set segids [lsort -unique [$sel get segid]]
foreach segid $segids {
    puts "Adding protein segid $segid to psfgen"
    set seg ${segid}
    set sel [atomselect top "resname POPE and segid $segid"]
    $sel set segid $seg
    $sel writepdb tmp.pdb
    segment $seg {pdb tmp.pdb}
    coordpdb tmp.pdb
}


if {$solute == 1} {

mol delete all
mol delete $sel 
set files [glob  "GOL*_f*.mol2"]
foreach mol $files {
	set molid [mol new $mol type mol2 waitfor all]
	puts "Loading $mol"
	set seltxt [atomselect top "resname GOL"]
	set seg [lindex [split $mol "_"] 0]
	$seltxt set segname $seg
	$seltxt writepdb ${seg}_vmd.pdb 
	
	segment $seg { pdb $seg.pdb }
	coordpdb $seg.pdb
	guesscoord   
	
	mol delete all
	mol delete $seltxt 
	mol delete $seg
	guesscoord   
}
} 

regenerate angles dihedrals 

;# Write files.
writepsf ${ofile}_top.psf ; # Output PSF.
writepdb ${ofile}_top.pdb ; # Output PDB.

;# This secction create box, solvate, remove bad waters, ionize and calculate box vectors and origin.
;# =================================================================================================

proc make_citosol { ofile ionconc } {

        ;# Calculate Size of membrane. 

	set id [mol new ${ofile}_top.psf type psf waitfor all]
	mol addfile ${ofile}_top.pdb type pdb waitfor all
	set seltxt [atomselect $id "protein or (resname POPE and type PL)"]
	set minmax [measure minmax $seltxt]

	set zlen [vecsub [lindex $minmax 1] [lindex $minmax 0]]
	set zboxpad [expr ([lindex $zlen 2] + 40)/2] ; # add 20 Angstrom more to Z(+/-) axis. 
	
	set xmin [lindex [lindex $minmax 0] 0] 
	set ymin [lindex [lindex $minmax 0] 1]
	set zmin [expr -1.0*$zboxpad] 
        set boxmin [list $xmin $ymin $zmin]
	
	set xmax [lindex [lindex $minmax 1] 0]
	set ymax [lindex [lindex $minmax 1] 1]
	set zmax [expr 1.0*$zboxpad] 
	set boxmax [list $xmax $ymax $zmax]

	;# Solvate system using size calculate above. 	

	package require solvate 
	solvate ${ofile}_top.psf ${ofile}_top.pdb -o ${ofile}_solvated  -b 1.5 -s WT -minmax [list $boxmin $boxmax] 

	;# Remove bad waters.

	mol delete all
	mol delete $id 

	set id [mol new ${ofile}_solvated.psf type psf waitfor all]
	mol addfile ${ofile}_solvated.pdb type pdb waitfor all 

	set all [atomselect $id all]
	set badw1 [atomselect top "water and same residue as within 3 of protein"]
	set badw2 [atomselect top "segid WT1 to WT99 and same residue as abs(z) < 25"]
	$badw1 num
	$badw2 num

	$badw1 set beta 1
	$badw2 set beta 1

	set allbadwater [atomselect top "name OH2 and beta > 0"]
	set seglistwater [$allbadwater get segid]
	set reslistwater [$allbadwater get resid]

	mol delete all
	package require psfgen
	resetpsf
	readpsf ${ofile}_solvated.psf
	coordpdb ${ofile}_solvated.pdb

	foreach segid $seglistwater resid $reslistwater {
   		delatom $segid $resid
	}

	writepsf ${ofile}_removeW.psf
	writepdb ${ofile}_removeW.pdb

	#; Calculate Cell box vectors. 

	mol delete all 
	set topid [mol new ${ofile}_removeW.psf type psf waitfor all] 
	mol addfile ${ofile}_removeW.pdb type pdb waitfor all
	set wat [atomselect $topid "waters"] ;## only membrane proteins 
	set all  [atomselect $topid "all"] 
	set minmax [measure minmax $wat]
	set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
	set center [measure center $all]

	set file [open "cellbox.inp" w+]
	puts $file "cellBasisVector1 [lindex $vec 0] 0 0"
	puts $file "cellBasisVector2 0 [lindex $vec 1] 0"
	puts $file "cellBasisVector3 0 0 [lindex $vec 2]"
	puts $file "CellOrigen $center"
	close $file
	
	puts "Done!!!!"

	;# Ionize system to an concentration of 0.150 mol/L. 
	
	mol delete all 
	package require autoionize 
	autoionize -psf ${ofile}_removeW.psf -pdb ${ofile}_removeW.pdb -cation SOD -anion CLA -sc $ionconc -o ${ofile}_150mM -seg ION
}

;# Call procedures. 
;#=====================================================
make_citosol  $ofile $ionconc

;# call Hyrogen Mass Repartition (HMR) if you need! 
;# ====================================================

if {$hmr == 1} {

puts "Start to Reparting Hydrogen!!!"

if { $oS == "MacOS"} {
	set vmd_exec "/Applications/VMD\ 1.9.4a48-Catalina-Rev7.app/Contents/MacOS/startup.command" 
	exec $vmd_exec -dispdev text -e do_hmr.tcl -args ${ofile}_150mM.psf ${ofile}_150mM.pdb 
	
} elseif {$oS == "Linux"} {
	exec vmd -dispdev text -e do_hmr.tcl -args ${ofile}_150mM.psf ${ofile}_150mM.pdb 
} else {
	puts "Error not operative system given"
}

} 

quit 

