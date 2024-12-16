#!/bin/bash

for kk in 2 4 6 8 10 12
do

cat > AlO.scf.kk.in <<EOF

&CONTROL
  calculation = 'scf'
  restart_mode = 'from_scratch'
  outdir = './tmp'
  pseudo_dir = './PSEUDO'
  prefix = 'AlO'
  tstress = .true.
  tprnfor = .true.
/
&SYSTEM
  ibrav = 2
  celldm(1) = 6.576
  nat = 2
  ntyp = 2
  ecutwfc = 80
  occupations = 'tetrahedra'
/
&ELECTRONS
  conv_thr = 1.0e-6
/
ATOMIC_SPECIES
  Al  26.982  Al_ONCV_PBE_sr.upf
  O  15.999  O_ONCV_PBE_sr.upf
ATOMIC_POSITIONS crystal
  Al   0.00   0.00   0.00
  O   0.50   0.50   0.50
K_POINTS automatic
  $kk $kk $kk 0 0 0
EOF

# Run the SCF calculation
$QESP/pw.x < AlO.scf.kk.in > output/AlO.scf.kk-$kk.out

# Put the energies in a file
grep -H '^!' output/AlO.scf.kk-*.out | awk -F ':' '{print $1, $2}' | awk '{print $1, $6}' > AlO.scf.kk.out

done