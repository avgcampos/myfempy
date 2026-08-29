from myfempy import newAnalysis
from myfempy import SteadyStateLinearIterative
from myfempy.io.iocsv import writer2csv

import numpy as np
# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(SteadyStateLinearIterative) 

mat = {
    "NAME": "mat",
    "VXY": 0.3,
    "EXX": 1000,
    }

geo = {
    "NAME": "geo",
    "THICKN": 10.0,
    }


# gmsh config
points = [
    [0, 0, 0],      # ponto 0
    [8, 0, 0],      # ponto A
    [128, 0, 0],    # ponto B
    [128, 32, 0],   # ponto C
    [0, 32, 0],     # ponto D
    [0, 8, 0],      # ponto E
]
         

lines = [[2, 3],        # line 1
         [3, 4],        # line 2
         [4, 5],        # line 3
         [5, 6],        # line 4
         ]

arcs = [
    [2, 1, 6],          # arco line 5
    ]  

plane = [
    [1, 2, 3, 4, 5],    # plano 1
        ]

NFINE = np.arange(20, 0, -1)
STR_XX_E = np.ones_like(NFINE)
NNODES = np.ones_like(NFINE)

for delem in range(len(NFINE)):

    modeldata = {

    "MESH": {
            'TYPE': 'gmsh',
            'filename': 'stress_conc',
            'pointlist': points,
            'linelist': lines,
            'planelist': plane,
            'arc': arcs,  
            'meshconfig': {
                'mesh': 'tria3', 
                'sizeelement': NFINE[delem],
                'meshmap': {'on': True,
                            'edge': 'all', #[[5]], #'all'
                            # "numbernodes": [50],
                }
                }
        },

        "ELEMENT": {
            'TYPE': 'structplane',
            'SHAPE': 'tria3',
            # 'INTGAUSS': 4,
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

    coord = fea.getCoord()

    force = {
        'TYPE': 'forceedge',
        'DOF': 'fx',
        'DIR': 'line',
        'TAG': 2,
        'VAL': [75]
    }

    bc_left = {
        'TYPE': 'fixed',
        'DOF': 'ux',
        'DIR': 'line',
        'TAG': 4,
    }

    bc_bottom = {
        'TYPE': 'fixed',
        'DOF': 'uy',
        'DIR': 'line',
        'TAG': 1
    }

    physicdata = {
        "PHYSIC": {"DOMAIN": "structural",
                "LOAD": [force],
                "BOUNDCOND": [bc_left, bc_bottom],
        },
        }
    fea.Physic(physicdata)

    # loadaply = fea.getLoadApply()
    # # print(loadaply)
    # # # print(fea.getLoadArray(loadaply))
    # print('forca total [N]', np.sum(loadaply[:,2]))

    # previewset = {'RENDER': {'filename': 'stress_conc', 'show': True, 'scale': 4, 'savepng': True, 'lines': False,
    #                          'plottags': {'line': True}
    #                          },
    #             #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
    #               }
    # fea.PreviewAnalysis(previewset)

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
                    "PLOTSET": {'show': True, 'filename': 'stress_conc', 'savepng': True},
                    # "TRACKER": {'point': {'x': 0, 'y': 0, 'z': 0, 'dof':1}},
                    "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
                }
    postprocdata = fea.PostProcess(postprocset)

    STR_XX_E[delem] = max(postprocdata['STRESS_XX'])
    NNODES[delem] = len(coord)

STR_XX_E_TEORICA = 242*np.ones_like(STR_XX_E)

data = [NNODES, STR_XX_E_TEORICA]
label = ['MESH', 'STR_XX_E_TEORICA']

# comando interno de io do myfempy para salvar arquivos CSV .txt para visualizar no ParaView
writer2csv('out/resultado_analitico.txt', data, label)  

data = [NNODES, STR_XX_E]
label = ['MESH', 'STR_XX_E']

# comando interno de io do myfempy para salvar arquivos CSV .txt para visualizar no ParaView
writer2csv('out/resultado_fem.txt', data, label)  