# strongback.tcl
# Test script for the strongback model
# Units: kip, ft, sec

#################################### Setup #####################################
set lib "/home/petertalley/Github/OpenSees-ShearBuildings/lib"

source [file join $lib {updateRayleighDamping.tcl}]

set storyMass [list 20.1242 20.1242 7.54658]
set g 32.2

#################################### Model #####################################
model BasicBuilder -ndm 2 -ndf 3

#----------------------------------- Nodes ------------------------------------#
# Ground
node 10  0  0
node 20  0  0
# Springs
node  1  0 20
node  2  0 35
node  3  0 50
# Trusses
node 11  0 20 -mass [lindex $storyMass 0] [lindex $storyMass 0] [lindex $storyMass 0]
node 12  0 35 -mass [lindex $storyMass 1] [lindex $storyMass 1] [lindex $storyMass 1]
node 13  0 50 -mass [lindex $storyMass 2] [lindex $storyMass 2] [lindex $storyMass 2]
# Strongback
node 21  0 20
node 22  0 35
node 23  0 50

#-------------------------------- Constraints ---------------------------------#
fix 10 1 1 1  ;# Trusses are fixed at the base
fix 20 1 1 0  ;# Strongback is pinned at the base

# Truss nodes need to be fixed in rotation
fix 11 0 0 1
fix 12 0 0 1
fix 13 0 0 1

# Attach springs to trusses
equalDOF 10 1 1 2 3
equalDOF 11 2 1 2 3
equalDOF 12 3 1 2 3

# Attach trusses to strongback -- not sure if this is the best way
equalDOF 11 21 1
equalDOF 12 22 1
equalDOF 13 23 1

#--------------------------------- Materials ----------------------------------#

# uniaxialMaterial Bilin $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S
#   $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg
#   $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg
#   $D_Plus $D_Neg <$nFactor>
uniaxialMaterial Bilin 1 3847.50 0.03 0.03 461.700 -461.700 10 10 10 10 1 1 1 1 1.00 1.00 1.500 1.500 0.3 0.3 52.4 52.4 1 1 0
uniaxialMaterial Bilin 2 3738.81 0.03 0.03 336.493 -336.493 10 10 10 10 1 1 1 1 0.75 0.75 1.125 1.125 0.3 0.3 39.3 39.3 1 1 0
uniaxialMaterial Bilin 3 1304.24 0.03 0.03 117.381 -117.381 10 10 10 10 1 1 1 1 0.75 0.75 1.125 1.125 0.3 0.3 39.3 39.3 1 1 0

# uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
uniaxialMaterial Elastic 11 769500
uniaxialMaterial Elastic 12 445500
uniaxialMaterial Elastic 13 121500

#---------------------------------- Elements ----------------------------------#

# Springs:
# element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2 ...
#   <-doRayleigh $rFlag> <-orient $x1 $x2 $x3 $yp1 $yp2 $yp3>
element zeroLength 1 1 11 -mat 1 -dir 1 -doRayleigh 1
element zeroLength 2 2 12 -mat 2 -dir 1 -doRayleigh 1
element zeroLength 3 3 13 -mat 3 -dir 1 -doRayleigh 1

# Trusses:
# element corotTruss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
element corotTruss 11 10 11 1 11
element corotTruss 12 11 12 1 12
element corotTruss 13 12 13 1 13

# Strongback itself:
# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag <-mass $massDens> <-cMass>
geomTransf Corotational 1
set A 1
set E 1e6
set Iz 1e6
element elasticBeamColumn 21 20 21 $A $E $Iz 1
element elasticBeamColumn 22 21 22 $A $E $Iz 1
element elasticBeamColumn 23 22 23 $A $E $Iz 1

################################# Gravity Loads ################################
pattern Plain 0 Linear {
  load 11 0 -[expr [lindex $storyMass 0]*$g] 0
  load 12 0 -[expr [lindex $storyMass 1]*$g] 0
  load 13 0 -[expr [lindex $storyMass 2]*$g] 0
}

system UmfPack
constraints Transformation
numberer RCM
test NormDispIncr 1e-6 10 1 2
algorithm KrylovNewton
integrator LoadControl 0.1
analysis Static

analyze 10
loadConst -time 0.0
wipeAnalysis

################################### Pushover ###################################
