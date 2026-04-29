import sys
# setting path
sys.path.append('../myfempy')

from myfempy import newAnalysis
from myfempy import HomogenizationPlane
from myfempy import PlaneStress

import numpy as np
import matplotlib.pyplot as plt
from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(HomogenizationPlane)

E_solid = 0.91
v = 0.3
r_solid = 1

mat1 = {
    "NAME": "solid",
    # "VXY": v,
    # "EXX": E_solid,       # N/mm^2 --> MPa
    "RHO": r_solid,       # kg/mm^2
    }


mat2 = {
    "NAME": "void",
    # "VXY": v,
    # "EXX": E_solid*1E-4,
    "RHO": r_solid*1E-4,
    }


geo = {"NAME": "geo",
       "THICKN": 1}

# ===============================================================================
#                    SET NEW MATERIAL <ORTHOTROPIC ELASTIC>
# ===============================================================================
class UserNewMaterial(PlaneStress):
    def __init__(self):
        super().__init__()
        
    def getMaterialSet():
        matset = {
            "mat": "planestress",
            "type": "homogenization",
        }
        return matset
    
    def getElasticTensor(tabmat, inci, element_number, a=None, b=None):
        # GUEDES (1989)
        if inci[element_number, 2] == 1:
            E1111 = 30
            E2222 = 30
            E1122 = 10
            E1212 = 10
        elif inci[element_number, 2] == 2:
            E1111 = 3
            E2222 = 3
            E1122 = 1
            E1212 = 1

        D = np.array([[E1111, E1122, 0],
                      [E1122, E2222, 0],
                      [0, 0, E1212]])
        return D

# # MODEL SET
La = 1   # mm
points = [
    [0, 0, 0],
    [La/2, 0, 0],
    [La, 0, 0],
    [0, La/2, 0],
    [La/2, La/2, 0],
    [La, La/2, 0],
    [0, La, 0],
    [La/2, La, 0],
    [La, La, 0]
]
         
lines = [
    [1, 2], # linha 1
    [2, 3], # linha 2
    [1, 4], # linha 3
    [2, 5], # linha 4
    [3, 6],
    [4, 7],
    [4, 5],
    [5, 8],
    [5, 6],
    [6, 9],
    [7, 8],
    [8, 9],
         ]

plane = [
    [1, 4, 7, 3],
    [2, 5, 9, 4],
    [9, 10, 12, 8],
    [7, 8, 11, 6],
         ]

esize = 0.01

modeldata = {

    "MESH": {
        'TYPE': 'gmsh',
        'filename': 'beso_myfempy_gmsh',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        'meshconfig': {
            'mesh': 'quad4',   #quad4
            'sizeelement': 2 * esize,
            'meshmap': {'on': True,
                        'edge': 'all',
            }
        }
    },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        # 'INTGAUSS': 8,
    },

    "MATERIAL": {
        "MAT": 'planestress',
        "TYPE": 'usermaterial',
        "CLASS": UserNewMaterial,
        # 'TYPE': 'isotropic',
        "PROPMAT": [mat1, mat2, mat1, mat2],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo, geo, geo, geo],
    },
}
fea.Model(modeldata)

# inci = fea.getInci()
# # coord = fea.getCoord()
# tabmat = fea.getTabmat()
# # tabgeo = fea.getTabgeo()
# # intgauss = fea.getIntGauss()

# print(tabmat)

# print(inci[:,2])

# sys.exit()

bc_X0_XX = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    'STEP': 1
    }

bc_X1_XX = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgex',
    'LOC': {'x': La, 'y': 999, 'z': 0},
    'STEP': 1
    }

bc_Y0_XX = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'STEP': 1
    }

bc_Y1_XX = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': La, 'z': 0},
    'STEP': 1
    }

bc_X0_YY = bc_X0_XX.copy()
bc_X0_YY['STEP'] = 2

bc_X1_YY = bc_X1_XX.copy()
bc_X1_YY['STEP'] = 2

bc_Y0_YY = bc_Y0_XX.copy()
bc_Y0_YY['STEP'] = 2

bc_Y1_YY = bc_Y1_XX.copy()
bc_Y1_YY['STEP'] = 2

bc_X0_XY = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    'STEP': 3
    }

bc_X1_XY = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgex',
    'LOC': {'x': La, 'y': 999, 'z': 0},
    'STEP': 3
    }

bc_Y0_XY = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'STEP': 3
    }

bc_Y1_XY = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': La, 'z': 0},
    'STEP': 3
    }

strainzero = {
    'TYPE': 'strainzero',
    'VAL': [[1,0,0], [0,1,0], [0,0,1]],    # mm/s^2 
    'DOF': 'none',
    'DIR':'none',
    'MESHNODE': 0,
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [strainzero],
               "BOUNDCOND": [
                    bc_X0_XX, bc_X1_XX, bc_Y0_XX, bc_Y1_XX,
                    bc_X0_YY, bc_X1_YY, bc_Y0_YY, bc_Y1_YY,
                    bc_X0_XY, bc_X1_XY, bc_Y0_XY, bc_Y1_XY
                             ],
    },
}
fea.Physic(physicdata)

loadaply = fea.getLoadApply()
print('forca total [N]', np.sum(loadaply[:,2]))

# sys.exit()

previewset = {'RENDER': {'filename': 'cell_preview', 'show': True, 'scale': 4, 'savepng': True, 'lines': False,
                        #  'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
# fea.PreviewAnalysis(previewset)

# sys.exit()
# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':True,
             'RHOH': True,
             }
solverdata = fea.Solve(solverset)

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'test_shakedown', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True, 'inci': True, 'coord': True}},
            }
postprocdata = fea.PostProcess(postprocset)

# print(solverdata["solution"]['U'])

CH = solverdata["solution"]['CH']
print('Homoge. Elastic Tensor\n', CH)


RH = solverdata["solution"]['RHOH']
print('Homoge. Density\n', RH)


# 2. Calcular a matriz de flexibilidade S (S = inv(C))
S = np.linalg.inv(CH)

# 3. Extrair as propriedades de engenharia
# No estado plano de tensões:
# S11 = 1/E1, S22 = 1/E2, S33 = 1/G12
# S12 = -nu12/E1 => nu12 = -S12 * E1
# S21 = -nu21/E2 => nu21 = -S21 * E2

E1 = 1 / S[0, 0]
E2 = 1 / S[1, 1]
G12 = 1 / S[2, 2]

nu12 = -S[0, 1] * E1
nu21 = -S[1, 0] * E2

# 4. Exibir resultados
print(f"--- Propriedades Calculadas ---")
print(f"Módulo de Young E1:  {E1:10.2f}")
print(f"Módulo de Young E2:  {E2:10.2f}")
print(f"Módulo de Cisalhamento G12: {G12:10.2f}")
print(f"Poisson nu12:        {nu12:10.4f}")
print(f"Poisson nu21:        {nu21:10.4f}")

# Verificação de simetria (Consistência termodinâmica)
print(f"\nVerificação (nu12/E1 == nu21/E2):")
print(f"{nu12/E1:.2e} == {nu21/E2:.2e}")