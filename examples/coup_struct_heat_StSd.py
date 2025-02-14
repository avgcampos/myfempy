
import numpy as np
import matplotlib.pyplot as plt

from myfempy import newAnalysis
from myfempy import SteadyStateLinear, SteadyStateLinearIterative

from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================

fea = newAnalysis(SteadyStateLinear) # StaticLinearIterative FASTER THEN StaticLinear

mat1 = {
    "NAME": "aco",
    "KXX": 0.0605,	 # W/mm·°C
    "KYY": 0.0605,
    "VXX": 0.3,
    "EXX": 200E3,       # N/mm^2 --> MPa
    "RHO": 7.85E-06,    # kg/mm^3
    "CTE": 12E-6    # (mm/mm) /°C
    }

geo = {
    "NAME": "espessura",
    "THICKN": 0.5, # mm
    # "DIM": [b, h, t, d],
    }

# MODEL SET
LX = 20 # 80   # mm
LY = 20 # 60
nelx = 20 # 40 # 80 # 160
nely = 20 # 30 # 60 # 120

# esize = 1 #LX/nelx

# MODEL SET
# example 10.4 Logan
# nodes = [[1, 0, 0, 0],
#          [2, LX, 0, 0],
#          [3, LX, LY, 0],
#          [4, 0, LY, 0]]

# conec = [[1, 1, 1, 1, 2, 3, 4]]

# gmsh config
# example 10.1 Logan
# points = [
#     [0, 0, 0],
#     [8, 0, 0],
#     [5, 4, 0],
#     [0, 4, 0]
# ]
         
# lines = [[1, 2],
#          [2, 3],
#          [3, 4],
#          [4, 1],
#          ]

# plane = [[1, 2, 3, 4]]

modeldata = {
    # "MESH": {
    #     'TYPE': 'add',
    #     'COORD': nodes,
    #     'INCI': conec,
    #     },

    "MESH": {
        'TYPE': 'legacy',
        'LX': LX,
        'LY': LY,
        'NX': nelx,
        'NY': nely,
        },

#    "MESH": {
#         'TYPE': 'gmsh',
#         'filename': 'test_line_force_x',
#         # "meshimport": 'object_dir',
#         'pointlist': points,
#         'linelist': lines,
#         'planelist': plane,
#         # 'arc': arcs,
#         'meshconfig': {
#             'mesh': 'quad4',   #quad4
#             'sizeelement': 2*esize,
#             # 'extrude': 20,
#             'meshmap': {'on': True,
#                         'edge': 'all', #'all'
#                         "numbernodes": 4,
#             }
#             }
#     },

    "ELEMENT": {
        'TYPE': 'heatplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'heatplane',
        "TYPE": 'isotropic',
        "PROPMAT": [mat1],
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

# element_number = 0
# st = time()
# ke = fea.getElemStifLinearMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# print('time KE ', time() - st)
# print(np.array2string(ke, separator=', '))
# # sys.exit()

# kg = fea.getGlobalMatrix(inci, coord, tabmat, tabgeo, intgauss)
# kg = kg['stiffness']
# print(kg.todense())

# import matplotlib.pylab as plt
# import scipy

# print('kg is sym ', scipy.linalg.issymmetric(kg.todense(), atol = 1e-09))
# # print('kg_cy is close to kg ',np.allclose(kg.todense(), kg_cy.todense()))

# plt.figure(2)
# plt.spy(kg.todense(), markersize=4)
# plt.show()

# sys.exit()
q = 1.0 #1
Q = 0.2         # W/mm3
h = 0.000005    # W/mm2.ºC 0.000005 (air)
Tinf = 20

heatflux = {
    'TYPE': 'heatfluxedge',
    'DOF': 'heatflux',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'VAL': [0], #q
    }

heatgen = {
    'TYPE': 'heatgeneration',
    'DOF': 'heatflux',
    'DIR': 'node',
    # 'LOC': {'x': 0, 'y': 999, 'z': 0},
    'VAL': [Q],
    }

convc_s1 = {
    'TYPE': 'convectionedge',
    'DOF': 'convection',
    'DIR': 'edgex',
    'LOC': {'x': LX, 'y': 999, 'z': 0},
    'VAL': [h*Tinf],
    }

convc_s2 = {
    'TYPE': 'convectionedge',
    'DOF': 'convection',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    'VAL': [h*Tinf],
    }

convc_s3 = {
    'TYPE': 'convectionedge',
    'DOF': 'convection',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': LY, 'z': 0},
    'VAL': [h*Tinf],
    }

bc1 = {
    'TYPE': 'insulated',  
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }

bc2 = {
    'TYPE': 'insulated',
    'DOF': 'full',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    }

bc3 = {
    'TYPE': 'insulated',
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': LX, 'y': 999, 'z': 0},
    }

bcNH1 = {
    'TYPE': 'temperature',
    'DOF': 't',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    'VAL': [100.0],
    }

bcNH2 = {
    'TYPE': 'temperature',
    'DOF': 't',
    'DIR': 'edgex',
    'LOC': {'x': LX, 'y': 999, 'z': 0},
    'VAL': [100.0, 0.0, 0.0],
    }

bcNH3 = {
    'TYPE': 'temperature',
    'DOF': 't',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'VAL': [100.0],
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "thermal",           # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [heatflux],
               "BOUNDCOND": [bcNH1],
    }
}
fea.Physic(physicdata)

# loadaply = fea.getLoadApply()
# # bcaply = fea.getBCApply()
# print(loadaply)
# # print(bcaply)
# print(fea.getLoadArray(loadaply))
# print(fea.getConstrains(bcaply))
# print(fea.getDirichletNH(bcaply))
# print('forca total [N]', np.sum(loadaply[:,2]))
# print('\n peso real', 9806.6*7.85E-06*LX*LY)

previewset = {'RENDER': {'filename': 'heat_transfer', 'show': True, 'scale': 10, 'savepng': True, 'lines': False,
                        #  'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# nodes = fea.getNodesFromRegions(1, 'plane')
# print(nodes)
# elem = fea.getElementFromNodesList(nodes)
# print(elem)

# sys.exit()

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

# print(solverdata['solution']['U'])

# sys.exit()

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'thermal': {'temp': True, 'heatflux': True}},
                "PLOTSET": {'show': True, 'filename': 'test_TSI_heat', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)


### new sim
modeldata = {
    "MESH": {
        'TYPE': 'legacy',
        'LX': LX,
        'LY': LY,
        'NX': nelx,
        'NY': nely,
        },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'planestress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat1],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

fc = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': LX/2, 'y': LY, 'z': 0},
    'VAL': [0.0],
    }

bcb = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    }

bcl = {
    'TYPE': 'fixed',
    'DOF': 'ux',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }

bcr = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': LX, 'y': 999, 'z': 0},
    }


physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [fc],
               "BOUNDCOND": [bcl, bcb],
    },
    "COUPLING": {"TYPE": 'thermalstress', # fsi asi
                 "POST": [postprocdata]
    },
    }
fea.Physic(physicdata)

# loadaply = fea.getCouplingInterface()
# print(loadaply)

previewset = {'RENDER': {'filename': 'structural', 'show': True, 'scale': 10, 'savepng': True, 'lines': False,
                        #  'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
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
                "PLOTSET": {'show': True, 'filename': 'test_TSI_struct', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)