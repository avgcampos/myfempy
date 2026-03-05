'''

Simulação estática em sólido com geração de malha com gmsh

Necessario instalação prévia do gmsh (não nativo do myfempy)

'''
# from os import environ
# environ['OMP_NUM_THREADS'] = '1'        #win
import numpy as np
import matplotlib.pyplot as plt

from myfempy import newAnalysis
from myfempy import SteadyStateLinear, SteadyStateLinearIterative

from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(SteadyStateLinear) # StaticLinearIterative FASTER THEN StaticLinear
# MODEL SET
mat1 = {
    "NAME": "aco",
    "VXY": 0.33,
    "EXX": 71E3,
    # "RHO": 1,    # kg/mm^2
    }

geo2 = {"NAME": "secao1",}

# MODEL SET
# mesh tetra 4
nodes = [[1, 0, 0, 0],
         [2, 1, 0, 0],
         [3, 0, 1, 0],
         [4, 0, 0, 1]]

# nodes = [[1, 500.,  500.,    0.],
#         [2,   0.,  500.,  250.],
#         [3,   0., 1000.,    0.],
#         [4,   0.,    0.,    0.]]
conec = [[1, 1, 1, 1, 2, 3, 4]]


# mesh hexa 8
nodes = [
    [1, 0, 0, 0],  # Nó 1
    [2, 1, 0, 0],  # Nó 2
    [3, 1, 1, 0],  # Nó 3
    [4, 0, 1, 0],  # Nó 4
    [5, 0, 0, 1],  # Nó 5
    [6, 1, 0, 1],  # Nó 6
    [7, 1, 1, 1],  # Nó 7
    [8, 0, 1, 1]   # Nó 8
]   
conec = [[1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8]]

# # gmsh config
LX = 200
LY = 20
LZ = 20

points = [
    [0, 0, 0],
    [LX, 0, 0],
    [LX, LY, 0],
    [0, LY, 0],
]
         
lines = [[1, 2],
         [2, 3],
         [3, 4],
         [4, 1],
         ]

plane = [[1, 2, 3, 4]]

modeldata = {
    # "MESH": {
    #     'TYPE': 'manual',
    #     'COORD': nodes,
    #     'INCI': conec,
    #     },

   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'data_mesh',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        'meshconfig': {
            'mesh': 'hexa8',
            # "numbernodes": 21,
            'sizeelement': 10,
            'extrude': LZ,
            'meshmap': {'on': True,
                        'edge': 'all',
                        "numbernodes": [20],
            }
            }
    },

    "ELEMENT": {
        'TYPE': 'structsolid',
        'SHAPE': 'hexa8',
        # 'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'solidelastic',
        "TYPE": 'isotropic',
        "PROPMAT": [mat1],
    },

    "GEOMETRY": {
        "GEO": 'solid',
        "PROPGEO": [geo2],
    },
}
fea.Model(modeldata)

inci = fea.getInci()
coord = fea.getCoord()
# tabmat = fea.getTabmat()
# tabgeo = fea.getTabgeo()
# intgauss = fea.getIntGauss()
# print(intgauss)

# vol = fea.getElementVolume(inci, coord, tabgeo)
# print('Volume ',np.sum(vol))

# reg = fea.getRegions()
# print(reg)


# sys.exit()

# element_number = 0
# st = time()
# ke = fea.getElemStifLinearMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# print('time KE ', time() - st)
# print(np.array2string(ke, separator=', '))
# print(np.allclose(ke, ke.T, rtol=1e-05, atol=1e-08))
# print(np.all(np.linalg.eigvals(ke)))


# st = time()
# kg = fea.getGlobalMatrix(inci, coord, tabmat, tabgeo, intgauss, SYMM=False)
# t_kg = time()-st
# print('time numpy sym ', t_kg)
# kg = kg['stiffness'].todense()
# print(kg.todense())
# print('kg is sym ', scipy.linalg.issymmetric(kg, atol = 1e-09))

# plt.figure(2)
# plt.spy(kg, markersize=4)
# plt.show()

# sys.exit()

f1 = {
    'TYPE': 'forcenode',
    'DOF': 'fx',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': 0, 'z': 0},
    'VAL': [1],
    }

f1 = {
    'TYPE': 'forcenode',
    'DOF': 'fx',
    'DIR': 'surfyz',
    'LOC': {'x': LX, 'y': 999, 'z': 999},
    'VAL': [1],
    }

f2 = {
    'TYPE': 'forcesurf',
    'DOF': 'fx',
    'DIR': 'surfyz',
    'LOC': {'x': LX, 'y': 999, 'z': 999},
    'VAL': [1],
    }

f2b = {
    'TYPE': 'forcesurf',
    'DOF': 'fy',
    'DIR': 'plane',
    'TAG': 4,
    'VAL': [-1.0],
    }

f3 = {
    'TYPE': 'forcesurf',
    'DOF': 'pressure',
    'DIR': 'surfyz',
    'LOC': {'x': LX, 'y': 999, 'z': 999},
    'VAL': [1],
    }

f3pt = {
    'TYPE': 'forcesurf',
    'DOF': 'pressure',
    'DIR': 'plane',
    'TAG': 4,
    'VAL': [-1.0],
    }

bc_node = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'node',
    'LOC': {'x': 1, 'y': 1, 'z': 1},
    }

bc = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'surfyz',
    'LOC': {'x': 0, 'y': 999, 'z': 999},
    }


physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [f2b],
               "BOUNDCOND": [bc],
    },
}
fea.Physic(physicdata)

# regions = fea.getRegions()
# print(regions)

# loadaply = fea.getLoadApply()
# # bcaply = fea.getBCApply()
# print(loadaply) 
# # print(bcaply)
# # print(fea.getLoadArray(loadaply))
# # print(fea.getDirichletNH(bcaply))
# print('forca total [N]', np.sum(loadaply[:,2]))

# sys.exit()

previewset = {'RENDER': {'filename': 'beam', 'show': True, 'scale': 2, 'savepng': True, 'lines': True,
                         'plottags': {'point': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# sys.exit()

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 1,
                        'step': 1},
             'SYMM':False,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

print(max(abs(solverdata['solution']['U'])))

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'test_shakedown', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)

# print(postprocdata["solution"])