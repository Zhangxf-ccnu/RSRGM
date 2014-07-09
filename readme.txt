README file for Matlab source code supporting the paper "Exploring overlapping functional units with various
structure in protein interaction networks".

Contents of this archive
------------------------
This archive contains several Matlab scripts used for revealing fucntional units in PPI networks using the 
algorithm RSRGM described in the above paper. 

1) RSRGM.m: Matlab script for the main algorithm of RSRGM described in the Fig. 2 of the above paper.

2) multi_RSRGM.m: Matlab script that repeats the entire calculation of RSRGM multiple times and choose the result
that gives the lowest value of the obejective function of (7). We also write the cohesive protein complexes revealed
by RSRGM to file  'cohesive_protein_complexes.txt'  and  the non-cohesive functional units revealed by RSRGM to file
'non_cohesive_functional_units.txt'. For both two files, each row corresponds an identified functional unit.

3) demo_RSRGM.m: A simple Matlab script to test RSRGM. When a PPI network is choosed, it can be run in a straightforward 
manner within a Matlab window.


This archive also contains a folder named as "data" which includes the four PPI networks used in the study. The four PPI 
networks are saved with Matlab .mat format. For each network, it is save as a structure that contains information of PPI 
network. Filed of "adjacent matrix" is the adjacent matrix PPI network and filed of "protein_list" is the name list of 
correspding proteins.

Please do not hesitate to contact Dao-Qing Dai at stsddq@mail.sysu.edu.cn (or Xiao-Fei Zhang at zhangfx9@mail2.sysu.edu.cn)
to seek any clarifications regarding any contents or operation of the archive.
