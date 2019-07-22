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
The structure will be processed by pdb2gmx, and you will be prompted to choose a force field:

```
GROMACS:      gmx pdb2gmx, version 2018.1
Executable:   /usr/bin/gmx
Data prefix:  /usr
Working dir:  /home/chris/Zika/Zika_MD
Command line:
  gmx pdb2gmx -f ZIKA.pdb -o ZIKA_processed.gro -water spce -ignh
Select the Force Field:
From '/usr/share/gromacs/top':
 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
 2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
 3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
 4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
 5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
 6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
 7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
 8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
 9: GROMOS96 43a1 force field
10: GROMOS96 43a2 force field (improved alkane dihedrals)
11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)
```

