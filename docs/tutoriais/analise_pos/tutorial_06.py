'''
myfempy Tutorial 06

Geração da malha com gmsh/ importação de modelos cad (.stp/ .step) e malha (.vtk) externa

Necessario instalação prévia do gmsh (não nativo do myfempy)

'''

# Imports

from myfempy import ModelGen
from myfempy import Solver
from myfempy import PostProcess
from myfempy import postproc_plot
from myfempy import preview_plot

import sys

# Definição do material, geometria 

mat = {
    "NAME": "material",
    "VXX": 0.25,
    "EXX": 200E6,
    "MAT": 'isotropic',
    "DEF": 'planestress'
    }

geo = {
    "NAME": "geometria",
    "THICKN": 1.0
    }


#---------------------------------------------


# # import cad .stp

# meshdata = {"GMSH": {'filename': 'tutorial_06a',
#                      'cadimport': {'object': 'tutorial_06_cad.stp'},
#                      'meshconfig': {'mesh': 'tetr4', 'elem': 'solid41', 'sizeelement': 10}},
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_06a',
#                          'show': True, 'scale': 10,
#                          'savepng': True,
#                          'lines': True,
#                          'plottags': {'point': False,
#                                       'edge': False}},
#               }

# preview_plot(previewset, modelinfo)


# ----------------------------------------------------

# import mesh .vtk

meshdata = {"GMSH": {'filename': 'tutorial_06b',
                     'meshimport':{'object':'tutorial_06_mesh'},
                     'meshconfig': {'mesh': 'tetr4', 'elem': 'solid41', 'sizeelement': 10}},
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            }

modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_06b',
                         'show': True, 'scale': 10,
                         'savepng': True,
                         'lines': True,
                         'plottags': {'point': False,
                                      'surf': True}},
              }

preview_plot(previewset, modelinfo)