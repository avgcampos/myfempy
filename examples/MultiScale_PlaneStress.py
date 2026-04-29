import sys
# setting path
sys.path.append('../myfempy')

import numpy as np

from myfempy import newAnalysis
from myfempy import SteadyStateLinear, DynamicEigenLinear
from myfempy import PlaneStress

# ===============================================================================
#                                   FEA STATIC
# ===============================================================================
fea = newAnalysis(SteadyStateLinear)

# #-------------------------------- PRE -------------------------------------#
mat = {
    "NAME": "dummy",
    }

geo = {
    "NAME": "espessura",
    "THICKN": 10, # mm
    }

# ===============================================================================
#      SET NEW MATERIAL <ELASTIC TENSOR FROM FILE: honeycomb.propmat.py>
# ===============================================================================

import honeycomb_propmat
D = honeycomb_propmat.solverdata["solution"]['CH']
rho_micro = honeycomb_propmat.solverdata["solution"]['RHOH']

def limpar_ruidos_e_simetrizar(C, tol=1e-2):
    """
    Limpa ruídos numéricos e força simetria em um tensor elástico ortotrópico em estado plano.
    
    Parâmetros:
    - C: matriz numpy 3x3
    - tol: tolerância para considerar valores como zero
    
    Retorna:
    - C_limpo: matriz numpy 3x3 ortotrópica
    """
    C_limpo = C.copy()
    
    # Zera valores muito pequenos
    C_limpo[np.abs(C_limpo) < tol] = 0.0
    
    # Força simetria: média dos termos C_ij e C_ji
    for i in range(C.shape[0]):
        for j in range(i+1, C.shape[1]):
            valor_medio = 0.5 * (C_limpo[i, j] + C_limpo[j, i])
            C_limpo[i, j] = valor_medio
            C_limpo[j, i] = valor_medio
    
    # Para ortotropia em estado plano, eliminamos acoplamentos com cisalhamento
    C_limpo[0, 2] = 0.0
    C_limpo[1, 2] = 0.0
    C_limpo[2, 0] = 0.0
    C_limpo[2, 1] = 0.0
    
    return C_limpo

CH = limpar_ruidos_e_simetrizar(D, tol=1e-2)

class HoneyCombPropMat(PlaneStress):
    def __init__(self):
        super().__init__()
        
    def getMaterialSet():
        matset = {
            "mat": "planestress",
            "type": "orthotropic",
        }
        return matset
    
    def getElasticTensor(tabmat, inci, element_number, a=None, b=None):
        return CH
    
    # def getFailureCriteria(stress):
    #     stress
    #     return super().getFailureCriteria()
# ===============================================================================

# MODEL SET
LX = 200
LY = 40

# gmsh config
points = [
        [0, 0, 0],
        [LX, 0, 0],
        [LX, LY, 0],
        [0, LY, 0]
        ]
         
lines = [
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 1],
        ]

plane = [
        [1, 2, 3, 4]
        ]

modeldata = dict()
modeldata["MESH"] = {
        'TYPE': 'gmsh',
        'filename': 'plane',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        'meshconfig': {
            'mesh': 'quad4',   #quad4
            'sizeelement': 2,
            'meshmap': {'on': True,
                        'edge': 'all', # [[1,3], [2,4]],
                        # "numbernodes": [100, 10],
            }
        }
    }

modeldata["ELEMENT"] = {
        'TYPE': 'structplane',
        'SHAPE': 'quad4',
        'INTGAUSS': 4,
    }

modeldata["MATERIAL"] = {
        "MAT": 'planestress',
        "TYPE": 'usermaterial',
        "PROPMAT": [mat],
        "CLASS": HoneyCombPropMat,

    }

modeldata["GEOMETRY"] =  {
        "GEO": 'thickness',
        "PROPGEO": [geo],
    }

fea.Model(modeldata)

f1 = {
    'TYPE': 'forceedge',
    'DOF': 'fy',
    'DIR': 'edgex',
    'LOC': {'x': 200, 'y': 999, 'z': 0},
    'VAL': [-0.25],
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
print('forca total [N]', np.sum(loadaply[:,2]))

previewset = {'RENDER': {'filename': 'model_preview', 'show': True, 'scale': 2, 'savepng': True, 'lines': False,
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
             'SYMM':False,
             }
solverdata = fea.Solve(solverset)

# print(solverdata['solution']['U'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'displ': True, 'stress': True}},
                "PLOTSET": {'show': True, 'filename': 'macro_static', 'savepng': True},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)


# ===============================================================================
#                                   FEA MODAL
# ===============================================================================
fea_modal = newAnalysis(DynamicEigenLinear)

mat = {
    "NAME": "dummy",
    'RHO': rho_micro,
}

# #-------------------------------- PRE -------------------------------------#

modeldata["MATERIAL"] = {
        "MAT": 'planestress',
        "TYPE": 'usermaterial',
        "PROPMAT": [mat],
        "CLASS": HoneyCombPropMat,
    }

fea_modal.Model(modeldata)

bc1 = {
    'TYPE': 'fixed',
    'DOF': 'full',
    'DIR': 'edgex',
    'LOC': {'x': 0, 'y': 999, 'z': 0},
    }  

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [],
               "BOUNDCOND": [bc1],
    },
    }
fea_modal.Physic(physicdata)

previewset = {'RENDER': {'filename': 'model_preview', 'show': True, 'scale': 2, 'savepng': True, 'lines': False,
                         'plottags': {'line': True}
                         },
            #   'LABELS': {'show': True, 'lines': True, 'scale': 1},
              }
fea_modal.PreviewAnalysis(previewset)

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'mode',  # mode, freq, time ...
                        'start': 0,
                        'end': 6,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea_modal.Solve(solverset)

print(solverdata['solution']['FREQ'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {'structural': {'modes': True}},
                "PLOTSET": {'show': True, 'filename': 'macro_modal', 'savepng': True},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)