gmx trjconv -f extended.part0002.xtc -e 200000 -o second.xtc
# ANALYSIS
gmx trjconv -f second.xtc -s NC_FME_*.pdb -pbc nojump -o trajout.xtc << EOF
0
EOF

item='abau'
for i in 10 20 30 40 60 ; do
scp wlodarsk@lumi.csc.fi:/scratch/project_465000291/tomek/ribosome_tunnels/methionine/${item}/${i}/fitted_all.xtc ${item}/${i}/
scp wlodarsk@lumi.csc.fi:/scratch/project_465000291/tomek/ribosome_tunnels/methionine/${item}/${i}/ref_1/fitted_all.xtc ${item}/${i}/ref_1/
scp wlodarsk@lumi.csc.fi:/scratch/project_465000291/tomek/ribosome_tunnels/methionine/${item}/${i}/ref_2/fitted_all.xtc ${item}/${i}/ref_2/
scp wlodarsk@lumi.csc.fi:/scratch/project_465000291/tomek/ribosome_tunnels/methionine/${item}/${i}/average_final_all.ccp4 ${item}/${i}/
scp wlodarsk@lumi.csc.fi:/scratch/project_465000291/tomek/ribosome_tunnels/methionine/${item}/${i}/ref_1/average_final_all.ccp4 ${item}/${i}/ref_1/
scp wlodarsk@lumi.csc.fi:/scratch/project_465000291/tomek/ribosome_tunnels/methionine/${item}/${i}/ref_2/average_final_all.ccp4 ${item}/${i}/ref_2/
done


name='4ybb'
for item in 10 20 30 40 60 ; do
cd ${item}
chimerax --nogui --cmd "open average_final_all.ccp4; save average_final_all.mrc format mrc; exit"
chimerax --nogui --cmd "open ref_1/average_final_all.ccp4; save ref_1/average_final_all.mrc format mrc; exit"
chimerax --nogui --cmd "open ref_2/average_final_all.ccp4; save ref_2/average_final_all.mrc format mrc; exit"
cd ..
done
