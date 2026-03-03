import sys
# setting path
sys.path.append('../myfempy')

from myfempy import newAnalysis
from myfempy import HomogenPlane
from myfempy import PlaneStress

import numpy as np
import matplotlib.pyplot as plt
from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(HomogenPlane)

mat = {
    "NAME": "material_dummy",
    "VXY": 0.33,
    "EXX": 70E3,    
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
        
        # BENDSOE (1988) and HASSANI (1998)
        D11 = 30
        D22 = 30
        D12 = 10
        D66 = 10

        D = np.array([[D11, D12, 0],
                      [D12, D22, 0],
                      [0, 0, D66]])
        return D

# # MODEL SET
La = 1   # mm
Lbx = 0.4   # mm 
Lby = 0.6   # mm 
points = [
    [0, 0, 0],
    [La, 0, 0],
    [La, La, 0],
    [0, La, 0],
    [(La-Lbx)/2, (La-Lby)/2, 0],
    [(La-Lbx)/2 + Lbx, (La-Lby)/2, 0],
    [(La-Lbx)/2 + Lbx, (La-Lby)/2 + Lby, 0],
    [(La-Lbx)/2, (La-Lby)/2 + Lby, 0],
]
         
lines = [
    [1, 2], # linha 1
    [2, 3], # linha 2
    [3, 4], # linha 3
    [4, 1], # linha 4
    [5, 6],
    [6, 7],
    [7, 8],
    [8, 5],
         ]

plane = [[1, 2, 3, 4],
         [-5, -6, -7, -8],
         ]

esize = 0.05

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
        "PROPMAT": [mat],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

# inci = fea.getInci()
# coord = fea.getCoord()
# tabmat = fea.getTabmat()
# tabgeo = fea.getTabgeo()
# intgauss = fea.getIntGauss()

# print(inci)

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

# loadaply = fea.getLoadApply()
# bcaply = fea.getBCApply()
# print(loadaply[:,2])

# freedof, fixedof, constdof = fea.getConstrains(bcaply)
# print('free',freedof)
# print('fixed',fixedof)
# print('const',constdof)
# print('forca total [N]', np.sum(loadaply[:,2]))

# sys.exit()

previewset = {'RENDER': {'filename': 'cell_preview', 'show': True, 'scale': 4, 'savepng': True, 'lines': False,
                        #  'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# sys.exit()

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':True,
             }
solverdata = fea.Solve(solverset)

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'test_shakedown', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True, 'inci': True, 'coord': True}},
            }
postprocdata = fea.PostProcess(postprocset)

# print(solverdata["solution"]['U'])

CH = solverdata["solution"]['CH']
print('Homoge. Elastic Tensor\n', CH)
