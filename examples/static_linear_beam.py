
import sys
# setting path
sys.path.append('../myfempy')

from myfempy import newAnalysis
from myfempy import SteadyStateLinear, SteadyStateLinearIterative

from time import time

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(SteadyStateLinear) # StaticLinearIterative FASTER THEN StaticLinear
LX = 10 
nelx = 32
# MODEL SET
mat1 = {
    "NAME": "aco",
    "EXX": 1000000,       # N/mm^2 --> MPa
    "GXY": 0.5*1000000,
    # "RHO": 1,    # kg/mm^2
    }

geo1 = {
    "NAME": "secao1",
    "AREACS": 1,
    "INERYY": 0.01,
    "INERZZ": 0.01,
    "INERXX": 0.02,
    }

# gmsh config
points = [
    [0, 0, 0],
    [1000, 0, 0],
    # [1000, 1000, 0],
    # [0, 1000, 0]
]
         
lines = [[1, 2],
        #  [3, 4],
        #  [4, 1],
         ]

modeldata = {
    # "MESH": {
    #     'TYPE': 'add',
    #     'COORD': nodes,
    #     'INCI': conec,
    #     },

    # "MESH": {
    #     'TYPE': 'legacy',
    #     'LX': LX,
    #     'NX': nelx,
    #     },

   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'data_mesh',
        'pointlist': points,
        'linelist': lines,
        'meshconfig': {
            'mesh': 'line2',
            "numbernodes": 101,
            }
    },

    "ELEMENT": {
        'TYPE': 'structbeam',
        'SHAPE': 'line2',
        # 'INTGAUSS': 8,
    },

    "MATERIAL": {
        "MAT": 'uniaxialstress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat1, mat1, mat1],
    },
    
    "GEOMETRY": {
        "GEO": 'userdefined',
        "PROPGEO": [geo1, geo1, geo1],
    },
}
fea.Model(modeldata)

inci = fea.getInci()
coord = fea.getCoord()
tabmat = fea.getTabmat()
tabgeo = fea.getTabgeo()
intgauss = fea.getIntGauss()

# print(tabgeo)

# sys.exit()

# vol = fea.getElementVolume(inci, coord, tabgeo)
# print(inci)

element_number = 0
# # st = time()
ke = fea.getElemStifLinearMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# me = fea.getElemMassConsistentMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# print('time KE ', time() - st)
# print(ke)
# print(me)
# print(np.array2string(me, separator=', '))
# print(np.allclose(ke, ke.T, rtol=1e-05, atol=1e-08))
# print(np.all(np.linalg.eigvals(ke) > 0))

# sys.exit()


# # element_number = 0
# # st = time()
# # ke = fea.getElemStifLinearMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# # print('time KE ', time() - st)
# print(np.array2string(ke, separator=', '))
# sys.exit()

# st = time()
# kg = fea.getGlobalMatrix(inci, coord, tabmat, tabgeo, intgauss)
# t_kg = time()-st
# print('time numpy sym ', t_kg)
# kg = kg['stiffness']
# print(kg.todense())

# # st = time()
# # # kg_cy = fea.experimental_asmbCython(inci, coord, tabmat, tabgeo, intgauss)
# # t_kg_cy = time()-st
# # print('time cython sym ', t_kg_cy)
# # # print(kg_cy.todense())

# # print(f'{t_kg/t_kg_cy=}')

# import matplotlib.pylab as plt
# import scipy

# print('kg is sym ', scipy.linalg.issymmetric(kg.todense(), atol = 1e-09))
# print('kg_cy is close to kg ',np.allclose(kg.todense(), kg_cy.todense()))

# plt.figure(2)
# plt.spy(kg.todense(), markersize=4)
# plt.show()

# sys.exit()

# peso = {
#     'TYPE': 'forcebody',
#     'DOF': 'fy',
#     'DIR': 'node',
#     # 'LOC': {'x': LX, 'y': 999, 'z': 0},
#     'VAL': [-9806.6],    #mm/s^2 
#     }


# f1 = {
#     'TYPE': 'forceedge',
#     'DOF': 'fx',
#     'DIR': 'edgex',
#     'LOC': {'x': LX, 'y': 999, 'z': 0},
#     'VAL': [100],
#     }

# f1c = {
#     'TYPE': 'forceedge',
#     'DOF': 'pressure',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': LY, 'z': 0},
#     'VAL': [100],
#     }

# f1d1 = {
#     'TYPE': 'forceedge',
#     'DOF': 'pressure',
#     'DIR': 'line',
#     'TAG': 2,
#     'VAL': [2000],
#     }

f1 = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': 1000, 'y': 0, 'z': 0},
    'VAL': [-10],
    }

f1b = {
    'TYPE': 'forcebeam',
    'DOF': 'fy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    'VAL': [-100],
    }

f2 = {
    'TYPE': 'forcenode',
    'DOF': 'tz',
    'DIR': 'node',
    'LOC': {'x': 10, 'y': 10, 'z': 0},
    'VAL': [5000],
    }

bc1 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': 0, 'z': 0},
    }

bc2 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'node',
    'LOC': {'x': LX, 'y': 0, 'z': 0},
    }


physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [f1],
               "BOUNDCOND": [bc1],
    },
}
fea.Physic(physicdata)

# loadaply = fea.getLoadApply()
# bcaply = fea.getBCApply()
# print(loadaply)
# # print(bcaply)
# print(fea.getLoadArray(loadaply))
# # print(fea.getConstrains(bcaply))
# # print(fea.getDirichletNH(bcaply))
# print('forca total [N]', np.sum(loadaply[:,2]))
# print('\n peso real', 9806.6*7.85E-06*LX*LY*0.1)

# sys.exit()

previewset = {'RENDER': {'filename': 'beam', 'show': True, 'scale': 20, 'savepng': True, 'lines': False,
                        #  'plottags': {'line': True}
                        # 'cs': True,
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

print(solverdata['solution']['U'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'test_shakedown', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)

# print(postprocdata["solution"])