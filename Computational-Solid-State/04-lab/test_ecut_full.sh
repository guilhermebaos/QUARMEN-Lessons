#!/bin/bash

# Check if energy_summary.dat exists and delete it if so
if [ -f energy_summary.dat ]; then
  rm energy_ecut.dat
fi


for ECUT in 20 30 40 50 60 80 100 120 160
do
cat > neon.scf.in <<EOF
&CONTROL
  calculation = 'scf'
  restart_mode = 'from_scratch'
  outdir = './tmp'
  pseudo_dir = './PSEUDO'
  prefix = 'neon'
  tstress = .true.
  tprnfor = .true.
/
&SYSTEM
  ibrav = 2
  celldm(1) = 8.43
  nat = 1
  ntyp = 1
  ecutwfc = $ECUT
  occupations = 'fixed'
/
&ELECTRONS
  conv_thr = 1.0e-6
/
ATOMIC_SPECIES
  Ne  20.18  Ne.LDA.upf
ATOMIC_POSITIONS crystal
  Ne   0.00   0.00   0.00
K_POINTS automatic
  4 4 4 0 0 0
EOF
# Run the SCF calculation
pw.x < neon.scf.in > neon.$ECUT.scf.out
# Extract the total energy and append to the summary file
  ETOT=$(grep '!' neon.$ECUT.scf.out | cut -d'=' -f2 | cut -d'R' -f1)
  echo "$ECUT $ETOT" >> energy_ecut.dat

done
