from myfempy import newAnalysis
from myfempy import SteadyStateLinearIterative, SteadyStateLinear
from myfempy.io.iocsv import writer2csv

import numpy as np
# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(SteadyStateLinearIterative)
mat = {
    "NAME": "mat",
    "VXY": 0.3,
    "EXX": 1000,       # N/mm^2 --> MPa
    }

geo = {
    "NAME": "geo",
    "THICKN": 10.0,
    # "DIM": [b, h, t, d],
    }

# gmsh config
points = [
    [0, 0, 0],    # ponto A
    [64, 0, 0],   # ponto B
    [64, 8, 0],   # ponto C   
    [0, 24, 0],   # ponto D
]


lines = [[1, 2],  # linha 1
         [2, 3],  # linha 2
         [3, 4], # linha 3
         [4, 1], # linha 4
         ]

plane = [[1, 2, 3, 4],
         ]

NFINE = 2

modeldata = {

   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'fem_tbar',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        'meshconfig': {
            'mesh': 'quad4',   #quad4
            'sizeelement': 1,
            'meshmap': {'on': True,
                        'edge': [[2, 4], [1, 3]], #'all'
                        "numbernodes": [NFINE*2, NFINE*4],
            }
            }
    },

    "ELEMENT": {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 1,
    },

    "MATERIAL": {
        "MAT": 'planestress',
        "TYPE": 'isotropic',
        "PROPMAT": [mat],
    },
    
    "GEOMETRY": {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

force = {
    'TYPE': 'forceedge',
    'DOF': 'fx',
    'DIR': 'line',
    'TAG': 2,
    'VAL': [0.625],
    }

bc_symm_x = {
    'TYPE': 'fixed',
    'DOF': 'uy',
    'DIR': 'line',
    'TAG': 1,
    }

bc_fixed_wall = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'line',
    'TAG': 4,
    }

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [force],
               "BOUNDCOND": [bc_symm_x, bc_fixed_wall],
    },
    }
fea.Physic(physicdata)

loadaply = fea.getLoadApply()
# print(loadaply)
# print(fea.getLoadArray(loadaply))
print('forca total [N]', np.sum(loadaply[:,2]))

previewset = {'RENDER': {'filename': 'fem_tbar', 'show': True, 'scale': 3, 'savepng': True, 'lines': False,
                         'plottags': {
                             'line': True
                            # 'point': True
                             }
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

print(np.max(np.abs(solverdata['solution']['U'])))

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'fem_tbar', 'savepng': True},
                # "PLOT": {'data': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)

# ===============================================================================
#                              RESULTADO ANALITICO
# ===============================================================================
# Parâmetros (substitua pelos valores reais)
F = 100.0    # força
l = 64.0     # comprimento
A0 = 16*10   # área inicial
E = 1000     # módulo de elasticidade

# Definição da função ux(x)
def ux(xn):
    return -((F * l) / (2 * A0 * E)) * np.log(1 - ((2 * xn) / (3 * l)))

# Intervalo de valores de x
x_vals = np.linspace(0, l, 200)
ux_vals = ux(x_vals)
print(np.max(np.abs(ux_vals)))

data = [x_vals, ux_vals]
label = ['x', 'u_x(x)']

# comando interno de io do myfempy para salvar arquivos CSV .txt para visualizar no ParaView
writer2csv('out/resultado_analitico.txt', data, label)  

