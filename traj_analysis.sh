# ANALYSIS
gmx trjconv -f md.xtc -s NC_FME_*.pdb -pbc nojump -o trajout.xtc << EOF
0
EOF

gmx make_ndx -f NC_FME_*.pdb -o index_r.ndx << EOF
r 166
q
EOF

gmx trjconv -s NC_FME_*.pdb -f trajout.xtc -n index_r.ndx -o fitted.xtc -fit rot+trans << EOF
14
0
q
EOF

/scratch/project_465000291/tomek/gromaps/bin/gmx_mpi  maptide -f fitted.xtc -s NC_FME_*.pdb -spacing 0.1 -margin 4.0 -mo average_final << EOF
0
EOF
