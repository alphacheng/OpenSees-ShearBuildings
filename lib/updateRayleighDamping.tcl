proc updateRayleighDamping { modeA ratioA modeB ratioB } {
  # ###################################################################
  # updateRayleighDamping $modeA $ratioA $modeB $ratioB
  # ###################################################################
  # Runs an eigenvalue analysis and set proportional damping based on
  # the current state of the structure
  #
  # Input Parameters:
  # modeA, modeB - modes that will have prescribed damping ratios
  # ratioA, ratioB - damping ratios prescribed at the specified modes

  # Get natural frequencies at the desired modes
  if { $modeA > $modeB } {
    set maxMode $modeA
  } else {
    set maxMode $modeB
  }

  set eigs    [eigen -fullGenLapack $maxMode]
  set freqA   [expr sqrt([lindex $eigs [expr $modeA-1]])]
  set freqB   [expr sqrt([lindex $eigs [expr $modeB-1]])]

  # Compute the damping factors
  set tempVal [expr 2.0/($freqA*$freqA-$freqB*$freqB)]
  set aM      [expr $tempVal*$freqA*$freqB*($ratioB*$freqA-$ratioA*$freqB)]
  set aK      [expr $tempVal*($ratioA*$freqA-$ratioB*$freqB)]

  # Set the damping
  rayleigh $aM 0.0 0.0 $aK
}
