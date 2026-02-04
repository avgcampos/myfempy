import sys
# setting path
sys.path.append('../myfempy')

import numpy as np

from myfempy import newAnalysis
from myfempy import SteadyStateLinear, SteadyStateLinearIterative
from myfempy import Material, PlaneStress

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(SteadyStateLinear)

# #-------------------------------- PRE -------------------------------------#
mat = {
    "NAME": "fibra_vidro", # Fonte: CALLISTER
    "EXX": 45E3,       # N/mm^2 --> MPa
    "EYY": 10E3,       # N/mm^2 --> MPa
    "VXY": 0.3,
    "VYZ": 0.4,
    }

geo = {
    "NAME": "espessura",
    "THICKN": 5, # mm
    }

# ===============================================================================
#                    SET NEW MATERIAL <ORTHOTROPIC ELASTIC>
# ===============================================================================
class UserNewMaterial(PlaneStress):
    def __init__(self):
        super().__init__()
        
    def getMaterialSet():
        matset = {
            "mat": "planestress",
            "type": "orthotropic",
        }
        return matset
    
    def getElasticTensor(Model, element_number):
        
        # # material elasticity
        EXX = Model.tabmat[int(Model.inci[element_number, 2]) - 1]["EXX"]
        EYY = Model.tabmat[int(Model.inci[element_number, 2]) - 1]["EYY"]
        # material poisson ratio
        VXY = Model.tabmat[int(Model.inci[element_number, 2]) - 1][ "VXY"]
        VYZ = Model.tabmat[int(Model.inci[element_number, 2]) - 1][ "VYZ"]  
        
        # EXX = 45E3       # N/mm^2 --> MPa
        # EYY = 12E3       # N/mm^2 --> MPa
        # VXY = 0.23
        # VYX = 0.66
                
        S00 = EXX/(1-VXY*VYZ)   
        S01 = (VYZ*EXX)/(1-VXY*VYZ)
        S10 = (VXY*EYY)/(1-VXY*VYZ)
        S11 = EYY/(1-VXY*VYZ)
        S22 = (EXX*EYY)/(EXX+EYY+2.0*EYY*VXY)     
                
        D = np.zeros((3, 3), dtype=np.float64)
        D[0, 0] = S00
        D[0, 1] = S01
        D[1, 0] = S10
        D[1, 1] = S11
        D[2, 2] = S22
        return D
    
    # def getFailureCriteria(stress):
    #     stress
    #     return super().getFailureCriteria()
# ===============================================================================

# config model geometry --> mesh(gmsh)
points = [
    [0, 0, 0],
    [200, 0, 0],
    [0, 10, 0],
    [200, 10, 0],
    [0, 20, 0],
    [200, 20, 0]
]

lines = [[1, 2],
         [2, 4],
         [3, 4],
         [3, 1],
         [6, 4],
         [5, 6],
         [5, 3],
         ]

plane = [[1, 2, 3, 4], [3, 5, 6, 7]]

# #-------------------------------- THERMAL SIM -------------------------------------#

modeldata = dict()
modeldata["MESH"] = {
        'TYPE': 'gmsh',
        'filename': 'test_line_force_x',
        'pointlist': points,
        'linelist': lines,
        'planelist': plane,
        # 'arc': arcs,
        'meshconfig': {
            'mesh': 'quad4',   #quad4
            'sizeelement': 1,
            'meshmap': {'on': True,
                        'edge': [[1,3,6], [2,4,5,7]], #'all'
                        "numbernodes": [20, 10],
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
        "PROPMAT": [mat, mat],
        "CLASS": UserNewMaterial,

    }

modeldata["GEOMETRY"] =  {
        "GEO": 'thickness',
        "PROPGEO": [geo, geo],
    }

fea.Model(modeldata)


f1 = {
    'TYPE': 'forcenode',
    'DOF': 'fx',
    'DIR': 'edgex',
    'LOC': {'x': 200, 'y': 999, 'z': 0},
    'VAL': [10.0],
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
                "PLOTSET": {'show': True, 'filename': 'solution', 'savepng': True},
                "OUTPUT": {'log': True, 'get': {'nelem': True, 'nnode': True}},
            }
postprocdata = fea.PostProcess(postprocset)
