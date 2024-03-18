
spinsys {
  channels 13C 15N
  nuclei   13C 15N
  dipole   1 2 895 10 20 30
  shift    1 10p 1p 0.5 50 20 10
}

par {
  variable index   1

  np               128
  spin_rate        10000
  proton_frequency 400e6
  start_operator   I1x
  detect_operator  I1p
  method           direct
  crystal_file     rep144
  gamma_angles     100
  sw               spin_rate/2
  variable tsw     1e6/sw
  verbose          1101
  variable rfF1    75000
  variable rfF2    50000
  variable t180F1  0.5e6/rfF1
  variable t180F2  0.5e6/rfF2
  variable Del_t180 0.5*(t180F2-t180F1)     
  variable tr2     0.5e6/spin_rate-t180F2
  variable tr1     0.5e6/spin_rate-0.5*t180F2
}

proc pulseq {} {
  global par
  store 4
  store 5

  reset
  delay $par(tr2)
  pulse $par(t180F2) 0 x $par(rfF2) x
  delay $par(tr2)
  pulse $par(t180F2) 0 x $par(rfF2) y
  store 1

  reset
  pulse $par(t180F2) 0 x $par(rfF2) y
  delay $par(tr2)
  pulse $par(t180F2) 0 x $par(rfF2) x
  delay $par(tr2)
  store 2

  reset
  delay $par(tr1)
  pulse $par(Del_t180) 0 x $par(rfF2) x
  pulse $par(t180F1) $par(rfF1) x $par(rfF2) x
  pulse $par(Del_t180) 0 x $par(rfF2) x
  delay $par(tr1)
  store 3

  acq
  for {set i 1} {$i < $par(np)} {incr i} {
    reset
    prop 5
    prop 2
    store 5

    reset
    prop 4
    prop 1
    store 4

    prop 3
    prop 5
    acq
  }
}
proc main {} {
  global par

  set f [fsimpson]
  fsave $f $par(name).fid
}
