'''

Simulação estática

'''

import numpy as np

from myfempy import newAnalysis
from myfempy import SteadyStateLinear, SteadyStateLinearIterative
# ===============================================================================
#                                   FEA
# ===============================================================================

fea = newAnalysis(SteadyStateLinear)

mat1 = {
    "NAME": "aco",
    "VXY": 0.3,
    "EXX": 200E6,       # N/mm^2 --> MPa
    "RHO": 7.85E-06,    # kg/mm^2
    }


mat2 = {
    "NAME": "alu",
    "VXY": 0.33,
    "EXX": 70E6,
    }


geo = {
    "NAME": "espessura",
    "THICKN": 1,
    # "DIM": [b, h, t, d],
    }

# MODEL SET
LX = 1.0 # 80
LY = 1.0 # 60
nelx = 2 # 40 # 80 # 160
nely = 2 # 30 # 60 # 120

esize = 0.004 #LX/nelx

# MODEL SET
# example 10.4 Logan
nodes = [[1, 0, 0, 0],
         [2, LX, 0, 0],
         [3, LX, LY, 0],
         [4, 0, LY, 0]]

conec = [[1, 1, 1, 1, 2, 3, 4]]

# gmsh config
# example 10.1 Logan
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
    # "MESH": {
    #     'TYPE': 'manual',
    #     'COORD': nodes,
    #     'INCI': conec,
    #     },

    # "MESH": {
    #     'TYPE': 'legacy',
    #     'LX': LX,
    #     'LY': LY,
    #     'NX': nelx,
    #     'NY': nely,
    #     },

   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'test_line_force_x',
        # "meshimport": 'object_dir',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        # 'arc': arcs,
        'meshconfig': {
            'mesh': 'tria3',   #quad4 tria3
            'sizeelement': 2*esize,
            # 'extrude': 20,
            'meshmap': {'on': True,
                        'edge': 'all', #[[1,3], [2,4]], #'all'
                        # "numbernodes": [20, 10],
            }
            }
    },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'tria3',
        # 'INTGAUSS': 1,
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
# coord = fea.getCoord()
# tabmat = fea.getTabmat()
# tabgeo = fea.getTabgeo()
# intgauss = fea.getIntGauss()
# print(type(intgauss))

# # vol = fea.getElementVolume(inci, coord, tabgeo)
# # print(np.sum(vol))

# np.set_printoptions(threshold=sys.maxsize)
# string_representation = np.array2string(coord, separator=', ')
# print(string_representation)
# print(coord)


# modelfem = fea.getModel()

# inci[0, 2] = 2
# inci[1, 2] = 2
# inci[2, 2] = 2

# modelfem.inci = inci

# print(modelfem.inci)

# element_number = 0
# st = time()
# ke = fea.getElemStifLinearMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# print('time KE ', time() - st)
# print(np.array2string(ke, separator=', '))
# print(np.allclose(ke, ke.T, rtol=1e-05, atol=1e-08))
# print(np.all(np.linalg.eigvals(ke) > 0))

# sys.exit()

# # element_number = 0
# # st = time()
# # ke = fea.getElemStifLinearMat(inci, coord, tabmat, tabgeo, intgauss, element_number)
# # print('time KE ', time() - st)
# print(np.array2string(ke, separator=', '))
# sys.exit()

# # st = time()
# kg = fea.getGlobalMatrix(modelfem)
# # t_kg = time()-st
# # print('time numpy sym ', t_kg)
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
# # print('kg_cy is close to kg ',np.allclose(kg.todense(), kg_cy.todense()))

# plt.figure(2)
# plt.spy(kg.todense(), markersize=4)
# plt.show()

# sys.exit()

peso = {
    'TYPE': 'forcebody',
    'DOF': 'fy',
    'DIR': 'node',
    # 'LOC': {'x': LX, 'y': 999, 'z': 0},
    'VAL': [-9806.6],    #mm/s^2 
    }


f1 = {
    'TYPE': 'forceedge',
    'DOF': 'fy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': LY, 'z': 0},
    'VAL': [-100],
    }

f1c = {
    'TYPE': 'forceedge',
    'DOF': 'pressure',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': LY, 'z': 0},
    'VAL': [-100],
    }

f1d1 = {
    'TYPE': 'forceedge',
    'DOF': 'pressure',
    'DIR': 'line',
    'TAG': 2,
    'VAL': [2000],
    }

f1d2 = {
    'TYPE': 'forceedge',
    'DOF': 'pressure',
    'DIR': 'line',
    'TAG': 3,
    'VAL': [2000],
    }

# f1b = {
#     'TYPE': 'forceedge',
#     'DOF': 'fy',
#     'DIR': 'line',
#     'TAG': 3,
#     'VAL': [-10],
#     }

f2 = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': LX, 'y': LY/2, 'z': 0},
    'VAL': [-100],
    }

f2b = {
    'TYPE': 'forcenode',
    'DOF': 'fy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': LY, 'z': 0},
    'VAL': [-100],
    }

f3 = {
    'TYPE': 'forcenode',
    'DOF': 'fx',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': 0, 'z': 0},
    'VAL': [1000],
    }


fs = {
    'TYPE': 'forcenode',
    'DOF': 'spring2ground',
    'DIR': 'node',
    'LOC': {'x': LX, 'y': 0, 'z': 0},
    'VAL': [100000.0],
    }


# bc2 = {
#     'TYPE': 'fixed',
#     'DOF': 'full',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     }


bc1 = {
    'TYPE': 'fixed',  
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }


bc2 = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 0, 'z': 0},
    }


bc3 = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'edgey',
    'LOC': {'x': 999, 'y': 10, 'z': 0},
    }


bcdnh = {
    'TYPE': 'displ',
    'DOF': 'ux',
    'DIR': 'edgex',
    'LOC': {'x': LX, 'y': 999, 'z': 0},
    'VAL': [0.025],
    }


# f_mbb = {
#     'TYPE': 'forcenode',
#     'DOF': 'fy',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 0, 'z': 0},
#     'VAL': [-100.0],
#     }

# bc_mbb1 = {
#     'TYPE': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     }

bc_mbb2 = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'node',
    'LOC': {'x': LX, 'y': LY, 'z': 0},
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural", # 'fluid' 'thermal'; "COUPLING": 'fsi'
               "LOAD": [f2b, fs],
               "BOUNDCOND": [bc1, bcdnh],
    },
}
fea.Physic(physicdata)

loadaply = fea.getLoadApply()
bcaply = fea.getBCApply()
print(loadaply)
# # print(bcaply)
# # print(fea.getLoadArray(loadaply))
# # print(fea.getDirichletNH(bcaply))
print('forca total [N]', np.sum(loadaply[:,2]))
# print('\n peso real', 9806.6*7.85E-06*LX*LY*0.1)]
# print(type(loadaply))

# sys.exit()

previewset = {'RENDER': {'filename': 'plane_stress', 'show': True, 'scale': 4, 'savepng': True, 'lines': True,
                        #  'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea.PreviewAnalysis(previewset)

# nodes = fea.getNodesFromRegions(1, 'plane')
# print(nodes)
# elem = fea.getElementFromNodesList(nodes)
# print(elem)



# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 1,
                        'step': 1},
            #  'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'test_shakedown', 'savepng': True},
                # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "REPORT": {'log': True,
                        'get':{
                            'nelem': True,
                            # 'nnode': True,
                            # 'inci': True,
                            # 'coord':True,
                            # 'tabmat':True,
                            # 'tabgeo':True,
                            # 'bc_list':True,
                            # 'lo_list':True,
                            # # 'numpy_decimals': 12,   # int()
                            # 'u_list': True,
                    }
            }}
postprocdata = fea.PostProcess(postprocset)

# print(solverdata['solution']['U'])
print(postprocdata["SOLUTION"])
