
##Generate the topology
#gmx pdb2gmx -f ZIKA.pdb -o ZIKA_processed.gro -water spce -ignh
## 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)

#gmx editconf -f ZIKA_processed.gro -o ZIKA_newbox.gro -c -d 2.5 -bt triclinic
#gmx solvate -cp ZIKA_newbox.gro -cs spc216.gro -o ZIKA_solv.gro -p topol.top
#gmx grompp -f ions.mdp -c ZIKA_solv.gro -p topol.top -o ions.tpr
#gmx genion -s ions.tpr -o ZIKA_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
##Chose Grup  3:SOL
##Group     3 (            SOL) has 383157 elements

#gmx grompp -f minim.mdp -c ZIKA_solv_ions.gro -p topol.top -o em.tpr

#gmx mdrun -v -deffnm em

#gmx energy -f em.edr -o potential.xvg

##change in the mdp file 
##RNA Water_and_ions
#gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

#gmx mdrun -deffnm nvt

#gmx energy -f nvt.edr -o temperature.xvg
## 18 0

#gmx energy -f npt.edr -o density.xvg
## 24 0
