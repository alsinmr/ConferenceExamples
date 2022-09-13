#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 13:19:30 2022

@author: albertsmith
"""

import pyDR
from pyDR.PCA import PCA


topo='/Volumes/My Book/HETs/HETs_3chain.pdb'
traj='/Volumes/My Book/HETs/MDSimulation/HETs_MET_4pw.xtc'
select=pyDR.MolSelect(topo=topo,traj_files=traj)

select.select_bond(Nuc='15N',segids='B')

pca=PCA(select)
pca.select_atoms('name N C CA')

select.tf=10000 #Truncate the trajectory for now


