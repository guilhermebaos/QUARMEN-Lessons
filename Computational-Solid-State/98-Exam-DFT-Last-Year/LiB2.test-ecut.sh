#!/bin/bash

for ECUT in 30 40 50 60 80 100
do

cat > LiB2.scf.ecut.in <<EOF

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
  ecutwfc = $ECUT
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
  6 6 6 0 0 0
EOF

# Run the SCF calculation
$QESP/pw.x < LiB2.scf.ecut.in > output/LiB2.scf.ecut-$ECUT.out

# Put the energies in a file
grep -H '^!' output/LiB2.scf.ecut-*.out | awk -F ':' '{print $1, $2}' | awk '{print $1, $6}' > LiB2.scf.ecut.out

done