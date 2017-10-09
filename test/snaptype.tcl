# Units: kip, ft, sec

source [file join {../lib/} {updateRayleighDamping.tcl}]

#################################### Model #####################################
model BasicBuilder -ndm 1 -ndf 1

#----------------------------------- Nodes ------------------------------------#
node 0 0
node 1 0 -mass 1
node 2 0 -mass 2
node 3 0 -mass 3

#-------------------------------- Constraints ---------------------------------#
fix 0 1

#--------------------------------- Materials ----------------------------------#
uniaxialMaterial Elastic 1 1 1 1
uniaxialMaterial Elastic 2 1 1 1
uniaxialMaterial Elastic 3 1 1 1

#---------------------------------- Elements ----------------------------------#
element zeroLength 1 0 1 -mat 1 -dir 1 -doRayleigh 1
element zeroLength 2 1 2 -mat 2 -dir 1 -doRayleigh 1
element zeroLength 3 2 3 -mat 3 -dir 1 -doRayleigh 1

