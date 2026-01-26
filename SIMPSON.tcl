
spinsys {
   channels 13C 1H
   nuclei 13C 1H
   dipole 1 2 -1500 0 0 0
}

par {
   crystal_file     alpha0beta0
   variable rfC     62500
   variable rfH     57500
   variable tsw     5
   variable index   2


   np               512
   proton_frequency 400e6
   start_operator   I2x
   detect_operator  I1x
   method           direct
   gamma_angles     1
   spin_rate        5000
   sw               spin_rate
   variable tsf     200e-6
   verbose          1
}

proc pulseq {} {
   global par

   acq_block {
      pulse 200e-6 $par(rfC) x $par(rfH) x
   }
}

proc main {} {
   global par

   set f [fsimpson]
   fsave $f $par(name).csdf
}
