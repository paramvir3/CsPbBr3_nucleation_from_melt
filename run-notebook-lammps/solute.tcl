#!/usr/bin/tclsh

if {[catch {package require topotools 1.1} ver]} {
	vmdcon -error "$ver. This script requires at least TopoTools v1.1. Exiting..."
		quit
}

if {[catch {package require pbctools 2.3} ver]} {
	vmdcon -error "$ver. This script requires at least pbctools v2.3. Exiting..."
		quit
}

set fname slab.pdb
# check for presence of coordinate file
if {! [file exists $fname]} {
	vmdcon -error "Required file '$fname' not available. Exiting..."
		quit
}

proc distance {new_pos x y z} {
    set bc [pbc get]
    set lx [lindex $bc 0 0]
    set ly [lindex $bc 0 1]
    set lz [lindex $bc 0 2]
    set x1 [lindex $new_pos 0]
    set y1 [lindex $new_pos 1]
    set z1 [lindex $new_pos 2]
    set dx [expr $x - $x1]
    set dy [expr $y - $y1]
    set dz [expr $z - $z1]
    if {$dx >   $lx * 0.5} { set dx [expr $dx - $lx] }
    if {$dx <= -$lx * 0.5} { set dx [expr $dx + $lx] }
    if {$dy >   $ly * 0.5} { set dy [expr $dy - $ly] }
    if {$dy <= -$ly * 0.5} { set dy [expr $dy + $ly] }
    if {$dz >   $lz * 0.5} { set dz [expr $dz - $lz] }
    if {$dz <= -$lz * 0.5} { set dz [expr $dz + $lz] }
    set dr [expr ($dx)**2+($dy)**2+($dz)**2]
    return $dr
}

mol new $fname autobonds no waitfor all
############### SOLUTE ##################
set selCs [atomselect top {name Cs}]
$selCs set type Cs
$selCs set mass 132.905
$selCs set charge 1.0
set numCs [$selCs num]
set posCs [$selCs get {x y z}]
set indexCs [$selCs get index]
set selPb [atomselect top {name Pb}]
$selPb set type Pb
$selPb set mass 204.20
$selPb set charge 2.0
set numPb [$selPb num]
set posPb [$selPb get {x y z}]
set indexPb [$selPb get index]
set selBr [atomselect top {name Br}]
$selBr set type Br
$selBr set mass 79.904 
$selBr set charge -1.0
set numBr [$selBr num]
set posBr [$selBr get {x y z}]
set indexBr [$selBr get index]

#######################################
set sel [atomselect top all]
set natoms [molinfo top get numatoms]
mol reanalyze top
pbc get
# we use a high-level tool from to multiply the system.
#TopoTools::replicatemol top 4 4 4

# and write out the result as a lammps data file.
topo writelammpsdata data.CPBr full
#topo writegmxtop topol.top

# for easier testing and visualization, we
# also write out copies in .pdb and .psf format.
animate write pdb CPBr.pdb
animate write psf CPBr.psf
# final sanity check the whole system has to be neutral.
$sel delete
set sel [atomselect top all]
set totq [vecsum [join [$sel get charge] { }]]
if {[expr {abs($totq)}] > 0.0005} {
	vmdcon -warning "Total system not neutral: $totq"
}

# done. now exit vmd
quit
