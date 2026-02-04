'''
myfempy Tutorial 02

Geração da malha por meio da opção legacy

'''

# Imports

from myfempy import ModelGen
from myfempy import Solver
from myfempy import PostProcess
from myfempy import postproc_plot
from myfempy import preview_plot
from math import pi

import sys

# Definição do material, geometria 

mat = {
    "NAME": "material",
    "VXX": 0.25,
    "EXX": 200E6,
    "MAT": 'isotropic',
    "DEF": 'axial'
    }


geo = {
    "NAME": "geo",
    "SEC":"I",
    'DIM':{'b':100,'h':150,'t':5,'d':5}}

# Linha (malha beam 2) 500 mm

meshdata = {"LEGACY": {'lx': 400, 'mesh': 'line2', 'elem': 'beam21', 'nx': 2},
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            # "FORCES": [force],
            # "BOUNDCOND": [bondcond],
            # "QUADRATURE": {'meth': 'gaussian', 'npp': 4},
            # "DOMAIN":'structural'
            }


modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_02a', 'show': True, 'scale': 1, 'savepng': True, 'lines': True, 'cs': True},
              }

preview_plot(previewset, modelinfo)


sys.exit()
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


# Retangulo estado plano (malha quad 4) 100 x 50 mm

meshdata = {"LEGACY": {'lx': 100, 'ly': 50, 'mesh': 'quad4', 'elem': 'plane41', 'nx': 10, 'ny': 5},
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            # "FORCES": [force],
            # "BOUNDCOND": [bondcond],
            # "QUADRATURE": {'meth': 'gaussian', 'npp': 4},
            # "DOMAIN":'structural'
            }


modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_02b', 'show': True, 'scale': 1, 'savepng': True, 'lines': True},
              }

preview_plot(previewset, modelinfo)


# Retangulo estado plano (malha tria 3) 100 x 50 mm

meshdata = {"LEGACY": {'lx': 100, 'ly': 50, 'mesh': 'tria3', 'elem': 'plane31', 'nx': 10, 'ny': 5},
            "PROPMAT": [mat],
            "PROPGEO": [geo],
            # "FORCES": [force],
            # "BOUNDCOND": [bondcond],
            # "QUADRATURE": {'meth': 'gaussian', 'npp': 4},
            # "DOMAIN":'structural'
            }


modelinfo = ModelGen.get_model(meshdata)

previewset = {'RENDER': {'filename': 'tutorial_02c', 'show': True, 'scale': 1, 'savepng': True, 'lines': True},
              }

preview_plot(previewset, modelinfo)