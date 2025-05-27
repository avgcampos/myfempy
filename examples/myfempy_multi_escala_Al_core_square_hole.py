import sys
# setting path
sys.path.append('../myfempy')

import numpy as np

from myfempy import newAnalysis
from myfempy import HomogenizationPlane


# ===============================================================================
#                    HOMOGENIZAÇAO MICRO ESCALA <BASE CELL>
# ===============================================================================
fea_micro_CH = newAnalysis(HomogenizationPlane)

# # Matrix Al
Ea = 71E3 # MPa @ https://en.wikipedia.org/wiki/Young%27s_modulus
va = 0.33 # @ https://en.wikipedia.org/wiki/Poisson%27s_ratio 

mat_micro = {
    "NAME": "Aluminium Alloy",
    "EXX": Ea,    
    "VXY": va,
    }

geo = {"NAME": "geo",
       "THICKN": 1}

# # MODEL SET
La = 10   # mm
Lbx = 4   # mm 
Lby = 6   # mm 
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
        #  [5, 6, 7, 8],
         ]

esize = 0.5

# modeldata = {

#     "MESH": {
#         'TYPE': 'gmsh',
#         'filename': 'beso_myfempy_gmsh',
#         'pointlist': points,
#         'linelist': lines,
#         'planelist': plane,
#         'meshconfig': {
#             'mesh': 'quad4',   #quad4
#             'sizeelement': 2 * esize,
#             'meshmap': {'on': True,
#                         'edge': 'all',
#             }
#         }
#     },

#     "ELEMENT": {
#         'TYPE': 'structplane',
#         'SHAPE': 'quad4',
#     },

#     "MATERIAL": {
#         "MAT": 'planestress',
#         "TYPE": 'isotropic',
#         "PROPMAT": [mat_micro],
#     },
    
#     "GEOMETRY": {
#         "GEO": 'thickness',
#         "PROPGEO": [geo],
#     },
# }
# fea_micro_CH.Model(modeldata)

# bc_X0_XX = {
#     'TYPE': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     'STEP': 1
#     }

# bc_X1_XX = {
#     'TYPE': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'edgex',
#     'LOC': {'x': La, 'y': 999, 'z': 0},
#     'STEP': 1
#     }

# bc_Y0_XX = {
#     'TYPE': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 0, 'z': 0},
#     'STEP': 1
#     }

# bc_Y1_XX = {
#     'TYPE': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': La, 'z': 0},
#     'STEP': 1
#     }

# bc_X0_YY = bc_X0_XX.copy()
# bc_X0_YY['STEP'] = 2

# bc_X1_YY = bc_X1_XX.copy()
# bc_X1_YY['STEP'] = 2

# bc_Y0_YY = bc_Y0_XX.copy()
# bc_Y0_YY['STEP'] = 2

# bc_Y1_YY = bc_Y1_XX.copy()
# bc_Y1_YY['STEP'] = 2

# bc_X0_XY = {
#     'TYPE': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     'STEP': 3
#     }

# bc_X1_XY = {
#     'TYPE': 'fixed',
#     'DOF': 'uy',
#     'DIR': 'edgex',
#     'LOC': {'x': La, 'y': 999, 'z': 0},
#     'STEP': 3
#     }

# bc_Y0_XY = {
#     'TYPE': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 0, 'z': 0},
#     'STEP': 3
#     }

# bc_Y1_XY = {
#     'TYPE': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': La, 'z': 0},
#     'STEP': 3
#     }

# strainzero = {
#     'TYPE': 'strainzero',
#     'VAL': [[1,0,0], [0,1,0], [0,0,1]],    # mm/s^2 
#     }

# physicdata = {
#     "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
#                "LOAD": [strainzero],
#                "BOUNDCOND": [
#                     bc_X0_XX, bc_X1_XX, bc_Y0_XX, bc_Y1_XX,
#                     bc_X0_YY, bc_X1_YY, bc_Y0_YY, bc_Y1_YY,
#                     bc_X0_XY, bc_X1_XY, bc_Y0_XY, bc_Y1_XY
#                              ],
#     },
# }
# fea_micro_CH.Physic(physicdata)

# previewset = {'RENDER': {'filename': 'cell_preview', 'show': True, 'scale': 4, 'savepng': True, 'lines': False,
#                          },
#               }
# fea_micro_CH.PreviewAnalysis(previewset)

# # #-------------------------------- SOLVER -------------------------------------#
# solverset = {"STEPSET": {'type': 'table',
#                         'start': 0,
#                         'end': 1,
#                         'step': 1},
#              'SYMM':True,
#              }
# solverdata = fea_micro_CH.Solve(solverset)

# postprocset = {"SOLVERDATA": solverdata,
#                 "COMPUTER": {'structural': {'displ': True, 'stress': True}},
#                 "PLOTSET": {'show': True, 'filename': 'micro_solution', 'savepng': True},
#                 # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
#                 "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
#             }
# postprocdata = fea_micro_CH.PostProcess(postprocset)
# CH = solverdata["solution"]['CH']
# print('Homoge. Elastic Tensor\n', CH)

# sys.exit()

# ===============================================================================
#                   FEA MACRO ESCALA
# ===============================================================================
from myfempy import PlaneStress
from myfempy import SteadyStateLinear

fea = newAnalysis(SteadyStateLinear)

mat_MACRO = {
    "NAME": "material_1",
    "VXY": 1,
    "EXX": 1,    
    }

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
        CH = np.array([[3.30799824e+04, 7.95908595e+03, 0.000000000000],
                       [7.95908595e+03, 4.46143001e+04, 0.000000000000],
                       [0.000000000000, 0.000000000000, 7.11995190e+03]])
        return CH
    

# MODEL SET
LX = 1000
LY = 100

esize = 10

points = [
    [0, 0, 0],
    [LX, 0, 0],
    [LX, LY, 0],
    [0, LY, 0]
]
         
lines = [[1, 2],
         [2, 3],
         [3, 4],
         [4, 1],
         ]

plane = [[1, 2, 3, 4]]


modeldata = {
    "MESH": {
        'TYPE': 'gmsh',
        'filename': 'myfempy_gmsh',
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
    },

    "MATERIAL": {
        "MAT": 'planestress',

        "TYPE": 'usermaterial',
        "CLASS": UserNewMaterial,
        "PROPMAT": [mat_MACRO],

        # 'TYPE': 'isotropic',
        # 'PROPMAT': [mat_micro]  # Aluminio puro!
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

f1 = {
    'TYPE': 'forceedge',
    'DOF': 'fy',
    'DIR': 'edgex',
    'LOC': {'x': LX, 'y': 999, 'z': 0},
    'VAL': [10],
    }


bc1 = {
    'TYPE': 'fixed',  
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }


physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [f1],
               "BOUNDCOND": [bc1],
    },
}
fea.Physic(physicdata)

loadaply = fea.getLoadApply()
print(np.sum(loadaply[:,2]))

previewset = {'RENDER': {'filename': 'plane_stress', 'show': True, 'scale': 4, 'savepng': True, 'lines': True,
                         },
              }
fea.PreviewAnalysis(previewset)


# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)


postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'macro_solution_homog', 'savepng': True},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)

print(solverdata['solution']['U'])