
from myfempy import newAnalysis
from myfempy import StaticLinearCyclicSymmPlane

import numpy as np
from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================

fea = newAnalysis(StaticLinearCyclicSymmPlane)

mat1 = {
    "NAME": "material_1",
    "VXX": 0.3,
    "EXX": 2E5,    # MPa
    "MAT": 'isotropic',
    "DEF": 'planestress'
    }

mat2 = {
    "NAME": "material_2",
    "VXX": 0.3,
    "EXX": 2E5,    # MPa
    "MAT": 'isotropic',
    "DEF": 'planestress'
    }

geo = {"NAME": "geo1",
       "THICKN": 10.0}


# MODEL SET 1/4 disc
RI = 50   # mm 
RE = 200

points = [
    [RI, 0, 0],
    [RE, 0, 0],
    [0, RI, 0],
    [0, RE, 0],
    [RI*np.cos(np.deg2rad(45)), RI*np.sin(np.deg2rad(45)), 0],
    [RE*np.cos(np.deg2rad(45)), RE*np.sin(np.deg2rad(45)), 0],
    [0, 0, 0]
]
         
lines = [
    [1, 2], # linha 1
    [5, 6], # linha 2
    [3, 4], # linha 3
         ]

# circle = [[RE, [0, 0, 0], ['0', 'Pi/2']], # linha 3
#         [RI, [0, 0, 0], ['0', 'Pi/2']]] # linha 4

arcs = [[1, 7, 5],  # arco 4
        [2, 7, 6],  # arco 5
        [5, 7, 3],  # arco 6
        [6, 7, 4],  # arco 7
         ]

plane = [[1, 5, 2, 4],
         [2, 7, 3, 6]]





# # MODEL SET 1/8 disc
# RI = 50   # mm 
# RE = 200

# points = [
#     [RI*np.cos(np.deg2rad(45)), RI*np.sin(np.deg2rad(45)), 0],
#     [RE*np.cos(np.deg2rad(45)), RE*np.sin(np.deg2rad(45)), 0],
#     [0, RI, 0],
#     [0, RE, 0],
#     [0, 0, 0]
# ]
         
# lines = [
#     [1, 2], # linha 1
#     [3, 4], # linha 2
#          ]

# arcs = [[3, 5, 1],  # arco 3
#         [4, 5, 2],  # arco 4
#         ]

# plane = [[1, 4, 2, 3],
#         ]

modeldata = {
   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'cs_disc',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        # 'circle': circle,  # 2 points degree begin --> end ()
        'arc': arcs,             # 3 points 
        'meshconfig': {
            'mesh': 'quad4',
            'sizeelement': 20,
            'meshmap': {'on': True,
                        'edge': [[1, 2, 3], [4, 5, 6, 7]], #'all',
                         "numbernodes": [1*12, 1*8],
            },
            }
    },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 8,
    },

    "MATERIAL": {
        "MAT": 'planestress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat1, mat2],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo, geo],
    },
}
fea.Model(modeldata)

# inci = fea.getInci()

# print(inci[0,:])
# print(inci[-1,:])

# regions = fea.getRegions()
# print(regions)

f1 = {
    'TYPE': 'forcenode',
    'DOF': 'fx',
    'DIR': 'point',
    'TAG': 3,
    'VAL': [100.0],
    }

f1bx = {
    'TYPE': 'forcenode',
    'DOF': 'fx',
    'DIR': 'point',
    'TAG': 6,
    'VAL': [100.0*np.cos(np.deg2rad(45))],
    }

f1by = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'point',
    'TAG': 6,
    'VAL': [-100.0*np.sin(np.deg2rad(45))],
    }

f2 = {
    'TYPE': 'forceedge',
    'DOF': 'fx',
    'DIR': 'line',
    'TAG': 4,
    'VAL': [10],
    }

f2b = {
    'TYPE': 'forceedge',
    'DOF': 'pressure',
    'DIR': 'line',
    'TAG': 5,
    'VAL': [10],
    }

f3b = {
    'TYPE': 'forceedge',
    'DOF': 'pressure',
    'DIR': 'line',
    'TAG': 7,
    'VAL': [10],
    }

bc_fixo1 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'line',
    'TAG': 4,
    }

bc_fixo2 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'line',
    'TAG': 6,
    }

bc1 = {
    'TYPE': 'fixed', 
    'DOF': 'ux',
    'DIR': 'line',
    'TAG': 2,
    }

bc2 = {
    'TYPE': 'fixed',  
    'DOF': 'uy',
    'DIR': 'line',
    'TAG': 1,
    }

bc_cs_l = {
    'TYPE': 'csymm',  # csleft
    'DOF': 'left',
    'DIR': 'line',
    'TAG': 3,
    }


bc_cs_r = {
    'TYPE': 'csymm',  # csleft
    'DOF': 'right',
    'DIR': 'line',
    'TAG': 1,
    }


physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [f2b, f3b],
               "BOUNDCOND": [bc_cs_l, bc_cs_r, bc_fixo1, bc_fixo2],
    },
}
fea.Physic(physicdata)

# loadaply = fea.getLoadApply()
# bcaply = fea.getBCApply()
# print(loadaply)
# print(bcaply)

# freedof, fixedof, constdof = fea.getConstrains(bcaply)
# print('free',freedof)
# print('fixed',fixedof)
# print('const',constdof)
# print('forca total [N]', np.sum(loadaply[:,2]))


previewset = {'RENDER': {'filename': 'cs_disc_preview', 'show': True, 'scale': 4, 'savepng': True, 'lines': False,
                         'plottags': {'line': True}
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

# print(solverdata['solution']['U'])

postprocset = {"SOLVERDATA": solverdata,
                "PLOTSET": {'show': True, 'filename': 'solution_quarter_cs', 'savepng': True},
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                # "TRACKER": {'point': {'x': 0.5*LX, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True, 'inci': True, 'coord': True}},
            }
postprocdata = fea.PostProcess(postprocset)
