# GROMACS tutorial for RNA molecular dynamics simulation

All the files you need to start MD on 3'end RNA of Zika Virus

This example will guide newbies through the process of setting up a simulation system containing an RNA molecule (3'end Zika virus genome
) in a box of water, with ions. Each step will contain an explanation of input and output, using typical settings for general use.

This tutorial assumes you are using a GROMACS version in the 2018.1

## Create topology

Verified that all the necessary atoms are present. Always check your pdb file for MISSING entries, as these entries indicate either atoms or whole residues that are not present in the crystal structure. Terminal regions may be absent, and may not present a problem for dynamics. Incomplete internal sequences or nucleic residues that have missing atoms will cause pdb2gmx to fail. These missing atoms/residues must be modeled in using other software packages: rosetta, robetta or molprobity. 

Please note that pdb2gmx is not magic. It cannot generate topologies for arbitrary molecules, just the residues defined by the force field (in the *.rtp files - generally proteins, nucleic acids, and a very finite amount of cofactors, like NAD(H) and ATP).
 
The purpose of pdb2gmx is to generate three files:

- The topology for the molecule.
- A position restraint file.
- A post-processed structure file.

The topology (topol.top by default) contains all the information necessary to define the molecule within a simulation. This information includes nonbonded parameters (atom types and charges) as well as bonded parameters (bonds, angles, and dihedrals). We will take a more detailed look at the topology once it has been generated.

Execute pdb2gmx by issuing the following command

`
gmx pdb2gmx -f ZIKA.pdb -o ZIKA_processed.gro -water spce -ignh
`
