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

# # Retangulo estado plano (malha manual quad 4) 100 x 50 mm com forcas e bc. nodais

# f1 = {
#     'DEF': 'forcenode',
#     'DOF': 'fy',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 0, 'z': 0},
#     'VAL': [-1000.0],
#     }

# f2 = {
#     'DEF': 'forcenode',
#     'DOF': 'fy',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 50, 'z': 0},
#     'VAL': [-1000.0],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 0, 'z': 0},
#     }

# bc2 = {
#     'DEF': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 50, 'z': 0},
#     }


# elementos = [[1,"plane41","material","geometria",[1, 2, 3, 4]]]


# coordenadas = [[1, 0, 0, 0],
#                [2, 100, 0, 0],
#                [3, 100, 50, 0],
#                [4, 0, 50, 0]]


# meshdata = {"ADD121": elementos,
#             "NODELIST": coordenadas,
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [f1,f2],
#             "BOUNDCOND": [bc1,bc2],
#             "QUADRATURE":{'meth':'gaussian','npp':4},
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_01a', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
#               'LABELS': {'show': True, 'lines': True, 'scale': 10},
#               }

# preview_plot(previewset, modelinfo)

#----------------------------------------------

# # Retangulo estado plano (malha manual quad 4) 100 x 50 mm com forcas distribuida e bc. nodais

# f1 = {
#     'DEF': 'forceedge',
#     'DOF': 'fy',
#     'DIR': 'edgex',
#     'LOC': {'x': 100, 'y': 999, 'z': 0},
#     'VAL': [-1000.0],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 0, 'z': 0},
#     }

# bc2 = {
#     'DEF': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 50, 'z': 0},
#     }


# elementos = [[1,"plane41","material","geometria",[1, 2, 3, 4]]]


# coordenadas = [[1, 0, 0, 0],
#                [2, 100, 0, 0],
#                [3, 100, 50, 0],
#                [4, 0, 50, 0]]


# meshdata = {"ADD121": elementos,
#             "NODELIST": coordenadas,
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [f1],
#             "BOUNDCOND": [bc1,bc2],
#             "QUADRATURE":{'meth':'gaussian','npp':4},
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_01b', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
#               'LABELS': {'show': True, 'lines': True, 'scale': 10},
#               }

# preview_plot(previewset, modelinfo)


#----------------------------------------------

# # Retangulo estado plano (malha legacy quad 4 10x5) 100 x 50 mm com forcas e bc distribuida

# f1 = {
#     'DEF': 'forceedge',
#     'DOF': 'fy',
#     'DIR': 'edgey',
#     'LOC': {'x': 999, 'y': 50, 'z': 0},
#     'VAL': [-1000.0],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'edgex',
#     'LOC': {'x': 0, 'y': 999, 'z': 0},
#     }

# bc2 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 0, 'z': 0},
#     }

# meshdata = {"LEGACY": {'lx': 100, 'ly': 50, 'mesh': 'quad4', 'elem': 'plane41', 'nx': 10, 'ny': 5},
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [f1],
#             "BOUNDCOND": [bc1,bc2],
#             "QUADRATURE": {'meth': 'gaussian', 'npp': 4},
#             # "DOMAIN":'structural'
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_01c', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
#               'LABELS': {'show': True, 'lines': True, 'scale': 10},
#               }

# preview_plot(previewset, modelinfo)


# --------------------------------

# # Retangulo estado plano (malha manual quad 4) 100 x 50 mm com massa e mola concentrada

# f1 = {
#     'DEF': 'forcenode',
#     'DOF': 'spring2ground',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 0, 'z': 0},
#     'VAL': [100.0],
#     }

# f2 = {
#     'DEF': 'forcenode',
#     'DOF': 'damper2ground',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 50, 'z': 0},
#     'VAL': [100.0],
#     }

# f3 = {
#     'DEF': 'forcenode',
#     'DOF': 'masspoint',
#     'DIR': 'node',
#     'LOC': {'x': 100, 'y': 50, 'z': 0},
#     'VAL': [1.0],
#     }

# bc1 = {
#     'DEF': 'fixed',
#     'DOF': 'all',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 0, 'z': 0},
#     }

# bc2 = {
#     'DEF': 'fixed',
#     'DOF': 'ux',
#     'DIR': 'node',
#     'LOC': {'x': 0, 'y': 50, 'z': 0},
#     }


# elementos = [[1,"plane41","material","geometria",[1, 2, 3, 4]]]


# coordenadas = [[1, 0, 0, 0],
#                [2, 100, 0, 0],
#                [3, 100, 50, 0],
#                [4, 0, 50, 0]]


# meshdata = {"ADD121": elementos,
#             "NODELIST": coordenadas,
#             "PROPMAT": [mat],
#             "PROPGEO": [geo],
#             "FORCES": [f1,f2,f3],
#             "BOUNDCOND": [bc1,bc2],
#             "QUADRATURE":{'meth':'gaussian','npp':4},
#             }

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'tutorial_01a', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
#             #   'LABELS': {'show': True, 'lines': True, 'scale': 10},
#               }

# preview_plot(previewset, modelinfo)


# ----------------------------

# Retangulo estado plano dois materiais (malha manual quad 4) 200 x 100 mm com forcas e bc. nodais

mat2 = {
    "NAME": "material2",
    "VXX": 0.25,
    "EXX": 200E6,
    "MAT": 'isotropic',
    "DEF": 'planestress'
    }

geo2 = {
    "NAME": "geometria2",
    "THICKN": 1.0
    }

f1 = {
    'DEF': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': 100, 'z': 0},
    'VAL': [-100.0],
    }

f2 = {
    'DEF': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': 100, 'y': 100, 'z': 0},
    'VAL': [-100.0],
    }

f3 = {
    'DEF': 'forcenode',
    'DOF': 'fy',
    'DIR': 'node',
    'LOC': {'x': 200, 'y': 100, 'z': 0},
    'VAL': [-100.0],
    }

bc1 = {
    'DEF': 'fixed',
    'DOF': 'all',
    'DIR': 'node',
    'LOC': {'x': 200, 'y': 0, 'z': 0},
    }

bc2 = {
    'DEF': 'fixed',
    'DOF': 'ux',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': 0, 'z': 0},
    }

bc3 = {
    'DEF': 'fixed',
    'DOF': 'ux',
    'DIR': 'node',
    'LOC': {'x': 0, 'y': 100, 'z': 0},
    }


elementos = [[1,"plane41","material","geometria",[1, 2, 3, 4]],
             [2,"plane41","material2","geometria2",[2, 5, 6, 3]]]


coordenadas = [[1, 0, 0, 0],
               [2, 100, 0, 0],
               [3, 100, 100, 0],
               [4, 0, 100, 0],
               [5, 200, 0, 0],
               [6, 200, 100, 0]]


meshdata = {"ADD121": elementos,
            "NODELIST": coordenadas,
            "PROPMAT": [mat,mat2],
            "PROPGEO": [geo,geo2],
            "FORCES": [f1,f2,f3],
            "BOUNDCOND": [bc1,bc2,bc3],
            "QUADRATURE":{'meth':'gaussian','npp':4},
            }

modelinfo = ModelGen.get_model(meshdata)

print(modelinfo['inci'])

previewset = {'RENDER': {'filename': 'tutorial_01a', 'show': True, 'scale': 10, 'savepng': True, 'lines': True},
              'LABELS': {'show': True, 'lines': True, 'scale': 10},
              }

preview_plot(previewset, modelinfo)

