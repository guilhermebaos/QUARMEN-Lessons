#!/bin/bash
for dx in -0.10 -0.08 -0.06 -0.04 -0.02 0.00 0.02 0.04 0.06 0.08 0.10
do
cat > si.movex.in <<EOF
&CONTROL
  calculation = 'scf'
  restart_mode = 'from_scratch'
  outdir = './tmp'
  pseudo_dir = './PSEUDO'
  prefix = 'cubicm'
/
&SYSTEM
  ibrav = 2
  celldm(1) = 10.26
  nat = 2
  ntyp = 1
  ecutwfc = 40 
  occupations = 'tetrahedra_opt'
/
&ELECTRONS
  conv_thr = 1.0e-6
/
ATOMIC_SPECIES
  Si  28.76  Si.LDA.upf
ATOMIC_POSITIONS 
Si   $dx   0.00   0.00
Si   0.25   0.25   0.25
K_POINTS automatic
4 4 4  0 0 0
EOF

# Run the SCF calculation
$QESP/pw.x < si.movex.in > si.movex.$dx.out

done
