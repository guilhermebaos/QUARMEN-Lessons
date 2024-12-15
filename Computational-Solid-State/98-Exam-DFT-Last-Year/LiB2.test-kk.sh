#!/bin/bash

for kk in 2 4 6 8 10
do

cat > LiB2.scf.kk.in <<EOF

&CONTROL
  calculation = 'scf'
  restart_mode = 'from_scratch'
  outdir = './tmp'
  pseudo_dir = './PSEUDO'
  prefix = 'LiB2'
  tstress = .true.
  tprnfor = .true.
/
&SYSTEM
  ibrav = 4
  celldm(1) = 5.826
  celldm(3) = 1.142
  nat = 3
  ntyp = 2
  ecutwfc = 60
  occupations = 'tetrahedra'
/
&ELECTRONS
  conv_thr = 1.0e-6
/
ATOMIC_SPECIES
  Li  6.94  Li_ONCV_PBE_sr.upf
  B  10.81  B_ONCV_PBE_sr.upf
ATOMIC_POSITIONS crystal
  Li   0.00   0.00   0.00
  B   1/3   2/3   1/2
  B   2/3   1/3   1/2
K_POINTS automatic
  $kk $kk $kk 0 0 0
EOF

# Run the SCF calculation
$QESP/pw.x < LiB2.scf.kk.in > output/LiB2.scf.kk-$kk.out

# Put the energies in a file
grep -H '^!' output/LiB2.scf.kk-*.out | awk -F ':' '{print $1, $2}' | awk '{print $1, $6}' > LiB2.scf.kk.out

done