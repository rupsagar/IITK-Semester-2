# UNITS: N, m, s
wipe all;

# MEMBER PARAMETERS
set DCol 0.5;
set BCol 0.3;
set HCol 3;
set NumFibY 20;
set NumFibZ 1;

# MATERIAL PARAMETERS
set Fy 415.0e6;
set Es 210.0e9;
set Bs 0.003;
set R0 20;
set cR1 0.925;
set cR2 0.15;

# SOME CALCULATED PARAMETERS
set Y [expr $DCol/2.0];
set Z [expr $BCol/2.0];

# DEFINE TAGS
set ColSecTag 1;
set ReinfMatTag 1;
set ColTransfTag 1;

# DEFINE MODEL
model BasicBuilder -ndm 2 -ndf 3;
set dir OSResults;
file mkdir $dir;

node 1 0 0;
node 2 0 $HCol;
set TopNode 2;

fix 1 1 1 1;
fix 2 0 0 0;

# GEOMETRIC TRANSFORMATION
geomTransf Linear $ColTransfTag;

# DEFINE STEEL MATERIAL
#set ReinfMat Steel01
#uniaxialMaterial $ReinfMat $ReinfMatTag $Fy $Es $Bs
set ReinfMat Steel02
uniaxialMaterial $ReinfMat $ReinfMatTag $Fy $Es $Bs $R0 $cR1 $cR2 0.05 1 0.07 1

# DEFINE FIBERS
section fiberSec $ColSecTag {patch rect $ReinfMatTag $NumFibY $NumFibZ [expr -$Y] [expr -$Z] $Y $Z};

# DEFINE ELEMENTS
set ElementType forceBeamColumn;
set numIntgrPts 5;
set IntType Lobatto;
element $ElementType 1 1 2 $ColTransfTag $IntType $ColSecTag $numIntgrPts;
set RecorderElement 1;

# DEFINE RECORDERS
recorder Node -file $dir/DispTop.txt -node $TopNode -dof 1 2 3 disp;
recorder Element -file $dir/ForceCol.txt -ele $RecorderElement force;
recorder Element -file $dir/ForceColSec1.txt -ele $RecorderElement section 1 force;
recorder Element -file $dir/DefoColSec1.txt -ele $RecorderElement section 1 deformation;

# LOAD CONTROLLED ANALYSIS - GRAVITY LOAD
pattern Plain 1 Linear {load $TopNode 0 -1 0};
constraints Plain;
numberer Plain;
system BandGeneral;
test EnergyIncr 1.0e-16 10;
algorithm Newton;
set PColMax 10000.0e3;
set PIncr 50.0e3;
set NumLoadStep [expr int($PColMax/$PIncr)];
integrator LoadControl $PIncr;
analysis Static;
analyze $NumLoadStep;
loadConst -time 0.0;

# DISPLACEMENT CONTROLLED ANALYSIS - TOP DISPLACEMENT
pattern Plain 2 Linear {load $TopNode 1 0 0};
constraints Plain;
numberer Plain;
system BandGeneral;
test EnergyIncr 1.0e-16 10;
algorithm Newton;
set DispCtrlDOF 1;
set Disp "0 0.4 -0.2 0.6 -0.4 -0.2 0"
set DispIncr 0.005;
set NumDisp [llength $Disp]
for {set Index 1;} {$Index<=[expr $NumDisp-1]} {incr Index 1} {
    set DispChng [expr [lindex $Disp $Index]-[lindex $Disp [expr $Index-1]]]
    if {$DispChng<0} {
    set IncrDecr -$DispIncr
        } else {
            set IncrDecr $DispIncr
    }
    integrator DisplacementControl $TopNode $DispCtrlDOF $IncrDecr
    analysis Static
    analyze [expr int($DispChng/$IncrDecr)]        
}