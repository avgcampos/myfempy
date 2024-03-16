'''
myfempy Tutorial 01

Geração da malha manual

'''

# Imports

from myfempy import ModelGen
from myfempy import Solver
from myfempy import PostProcess
from myfempy import postproc_plot
from myfempy import preview_plot

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

# force = {
#     'DEF': 'forceedge',
#     'DOF': 'fx',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     'VAL': [1.0],
#     }

# bondcond = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 0, 'z': 0},
#     }

# Retangulo estado plano (malha quad 4) 100 x 50 mm

elementos = [[1,"plane41","material","geometria",[1, 2, 3, 4]]]


coordenadas = [[1, 0, 0, 0],
               [2, 100, 0, 0],
               [3, 100, 50, 0],
               [4, 0, 50, 0]]


meshdata = {"ADD121": elementos,
            "NODELIST": coordenadas,
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            # "FORCES": [force],
            # "BOUNDCOND": [bondcond],
            # "QUADRATURE":{'meth':'gaussian','npp':8},
            }

modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_01a', 'show': True, 'scale': 1, 'savepng': True, 'lines': True},
              }

preview_plot(previewset, modelinfo)

# Retangulo estado plano (malha tria 4) 50 x 50 mm/ teste de qualidade 'Aspect Ratio'

elementos = [[1,"plane31","material","geometria",[1, 2, 4]],
             [2,"plane31","material","geometria",[2, 3, 4]]]


coordenadas = [[1, 0, 0, 0],
               [2, 50, 0, 0],
               [3, 50, 40, 0],
               [4, 0, 50, 0]]


meshdata = {"ADD121": elementos,
            "NODELIST": coordenadas,
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            # "FORCES": [force],
            # "BOUNDCOND": [bondcond],
            # "QUADRATURE":{'meth':'gaussian','npp':8},
            }

modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_01b', 'show': True, 'scale': 1, 'savepng': True, 'lines': True},
              'QUALITY': {'show': True, 'method': 1, 'scale': 0.1},
              }

preview_plot(previewset, modelinfo)

# Solido 100 x 100 x 100 mm

mat2 = {
    "NAME": "material2",
    "VXX": 0.25,
    "EXX": 200E6,
    "MAT": 'isotropic',
    "DEF": 'solid'
    }

geo2 = {
    "NAME": "geometria2",
    "THICKN": 0.0
    }

elementos = [[1,"solid81","material2","geometria2",[1, 2, 3, 4, 5, 6, 7, 8]]]

coordenadas = [[1, 0, 0, 0],
               [2, 100, 0, 0],
               [3, 100, 100, 0],
               [4, 0, 100, 0],
               [5, 0, 0, 100],
               [6, 100, 0, 100],
               [7, 100, 100, 100],
               [8, 0, 100, 100]]

meshdata = {"ADD121": elementos,
            "NODELIST": coordenadas,
            "PROPMAT": [mat2],
            "PROPGEO": [geo2],
            # "FORCES": [force],
            # "BOUNDCOND": [bondcond],
            # "QUADRATURE":{'meth':'gaussian','npp':8},
            }

modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_01c', 'show': True, 'scale': 1, 'savepng': True, 'lines': True},
              'LABELS': {'show': True, 'lines': True, 'scale': 10},
              }

preview_plot(previewset, modelinfo)