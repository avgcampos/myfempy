# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 13:07:28 2022

@author: viniv
"""
import numpy as np
from myfempy.felib.struct.plane32 import plane_t3_stif, plane_t3_stress
from myfempy.solver.setup import load_apply, bond_condi
from myfempy.solver.static import solve
from myfempy.io.ioctrl import export_vtkmesh

# coord = np.array([[1, 0, 0, 0],\
#                   [2, 1, 0, 0],\
#                   [3, 1, 1, 0],\
#                   [4, 0, 1, 0]])
    
# inci = np.array([[1,2,1,1,1,2,4],\
#                  [2,2,1,1,2,3,4]])

inci = np.array([[1,2,1,1,1,6,5],
                 [2,2,1,1,1,2,6],
                 [3,2,1,1,2,7,6],
                 [4,2,1,1,2,3,7],
                 [5,2,1,1,3,8,7],
                 [6,2,1,1,3,4,8],
                 [7,2,1,1,5,10,9],
                 [8,2,1,1,5,6,10],
                 [9,2,1,1,6,11,10],
                 [10,2,1,1,6,7,11],
                 [11,2,1,1,7,12,11],
                 [12,2,1,1,7,8,12]])
    
    
coord = np.array([[1,0,0,0],
                 [2,2,0,0],
                 [3,4,0,0],
                 [4,6,0,0],
                 [5,0,1.5,0],
                 [6,2,1.5,0],
                 [7,4,1.5,0],
                 [8,6,1.5,0],
                 [9,0,3,0],
                 [10,2,3,0],
                 [11,4,3,0],
                 [12,6,3,0]])

tabmat = np.array([[200E9, 0.33, 100E9, 0, 0, 0, 7800]])

tabgeo = np.array([[0, 0, 0, 0, 0.05]])

# forcenodeaply = np.array([[2, 1, 1000, 1],\
#                           [3, 1, 1000, 1]])

# boncdnodeaply = np.array([[0, 1],\
#                           [0, 4]])

forcenodeaply = np.array([[12, 2, -1500, 1],\
                          [8, 2, -1500, 1],\
                          [4, 2, -1500, 1]])

boncdnodeaply = np.array([[0, 1],\
                          [0, 5],\
                          [0, 9]])

KEYELEM = "plane32"
NDOF = 2
datamesh = {'nodedof' : NDOF,
            'lencoord': len(coord),
            'leninci': len(inci),
            'fulldof': NDOF*len(coord),
            'keyelem': KEYELEM,
            'matdefi':'isotropic',
            'matbeha':'planestress',
            'stresscal': 'avr'}

KG = plane_t3_stif(datamesh, inci, coord, tabmat, tabgeo)

freedof, fixedof = bond_condi(datamesh, boncdnodeaply)

forcelist = load_apply(datamesh, forcenodeaply)

solverset = {'step' : 0,
             'type' : 'iterative',
             'tol': 2E-4}

U = solve(datamesh, forcelist, freedof, KG, solverset)

stress_list, strain_list = plane_t3_stress(U, inci, coord, tabmat, datamesh)

data_result = np.concatenate((stress_list,strain_list),axis=1)
postdata = {'filename' : 'postprocess',
            'vtkcelltype' : 5,              # vtk cell type --> triangule
            'datalist' : data_result,
            'datatitle' :  ['STRESS_VM','STRESS_XX','STRESS_YY','STRESS_XY','STRAIN_VM','STRAIN_XX','STRAIN_YY','STRAIN_XY'],
            'datatype' : datamesh['stresscal']}

export_vtkmesh(postdata, datamesh, coord, inci)