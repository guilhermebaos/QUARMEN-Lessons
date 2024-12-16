#!/bin/bash

for ECUT in 20 30 40 50 60 80 100 120
do

cat > AlO.scf.ecut.in <<EOF

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
  ecutwfc = $ECUT
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
  8 8 8 0 0 0
EOF

# Run the SCF calculation
$QESP/pw.x < AlO.scf.ecut.in > output/AlO.scf.ecut-$ECUT.out

# Put the energies in a file
grep -H '^!' output/AlO.scf.ecut-*.out | awk -F ':' '{print $1, $2}' | awk '{print $1, $6}' > AlO.scf.ecut.out

done