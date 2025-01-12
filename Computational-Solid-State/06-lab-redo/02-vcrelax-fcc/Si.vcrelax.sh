#!/bin/bash

for PT in 0 100 200 300 400 500 600 700 800
do

cat > Si.vcrelax.$PT.in <<EOF

&CONTROL
  calculation = 'vc-relax'
  restart_mode = 'from_scratch'
  outdir = './tmp'
  pseudo_dir = './PSEUDO'
  prefix = 'vcrelax'
  tstress = .true.
  tprnfor = .true.
/
&SYSTEM
  ibrav = 2
  celldm(1) = 10.21
  nat = 1
  ntyp = 1
  ecutwfc = 60 
  occupations = 'tetrahedra_opt'
/
&ELECTRONS
  conv_thr = 1.0e-6
/
&IONS
  ion_dynamics = 'bfgs', trust_radius_ini = 0.1 
/
&CELL
  cell_dynamics='bfgs'
  press=$PT
/
ATOMIC_SPECIES
  Si  28.76  Si.LDA.upf
ATOMIC_POSITIONS 
  Si   0.00   0.00   0.00
K_POINTS automatic
  6 6 6  0 0 0
EOF

done