
spinsys {
  channels 13C 15N
  nuclei   13C 15N
  dipole   1 2 895 0 0 0
}

par {
  variable index   1
  np               64
  spin_rate        15000
  proton_frequency 400e6
  start_operator   I1x
  detect_operator  I1x
  method           direct
  crystal_file     rep144
  gamma_angles     100
  sw               spin_rate/2
  variable tsw     1e6/sw
  variable tr      1e6/spin_rate
  variable rf      100000
  variable pi     0.5e6/rf
  variable dly     0.5*tr-pi
  variable dly_ref 0.5*(tr-pi)
}

proc pulseq {} {
  global par
  reset
  
  
  store 0
  store 1
  
  delay $par(dly)
  pulse $par(pi) 0 x $par(rf) x
  delay $par(dly)
  pulse $par(pi) 0 x $par(rf) y
  store 2
  
  reset
  pulse $par(pi) 0 x $par(rf) y
  delay $par(dly)
  pulse $par(pi) 0 x $par(rf) x
  delay $par(dly)
  store 3
  
  reset
  delay $par(dly_ref)
  pulse $par(pi) $par(rf) x 0 x
  delay $par(dly_ref)
  store 4
  

  reset
  acq
  for {set i 1} {$i < $par(np)} {incr i} {
    prop 0
    prop 2
    store 0
    reset
    
    prop 1
    prop 3
    store 1
    reset
    
    prop 0
    prop 4
    prop 1
    acq
  }
}
proc main {} {
  global par

  set f [fsimpson]
  fsave $f $par(name).fid
}
