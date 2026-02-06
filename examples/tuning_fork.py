from myfempy import newAnalysis
from myfempy import DynamicEigenLinear

# ===============================================================================
#                                   FEA
# ===============================================================================
fea = newAnalysis(DynamicEigenLinear, 'sim')
# MODEL SET
mat = {
    "NAME": "Aluminum_Alloy",
    "VXY": 0.33,
    "EXX": 71E6,        # MPa (N/mm2)
    "RHO": 2.77E-6,    # kg/mm3
    }

geo = {"NAME": "Solid"}

modeldata = {
   "MESH": {
        'TYPE': 'gmsh',
        'filename': 'mesh_gmsh',
        # !!! COLOCAR O CAD NO DIRETORIO DA SIMULACAO !!!
        'cadimport': {'object': 'tuning_fork.STEP'}, 
        'meshconfig': {
            'mesh': 'tetr4',
            'sizeelement': 2,
            },
    },

    "ELEMENT": {
        'TYPE': 'structsolid',
        'SHAPE': 'tetr4',
        # 'INTGAUSS': 4,
    },

    "MATERIAL": {
        "MAT": 'solidelastic',
        "TYPE": 'isotropic',
        "PROPMAT": [mat],
    },

    "GEOMETRY": {
        "GEO": 'solid',
        "PROPGEO": [geo],
    },
}
fea.Model(modeldata)

physicdata = {
    "PHYSIC": {"DOMAIN": "structural",
               "LOAD": [],
               "BOUNDCOND": [],
    },
}
fea.Physic(physicdata)

previewset = {'RENDER': {'filename': 'preview', 'show': True, 'scale': 4, 'savepng': True, 'lines': True,
                        #  'plottags': {'plane': True}
                         },
              }
fea.PreviewAnalysis(previewset)

# #-------------------------------- SOLVER -------------------------------------#
solverset = {"STEPSET": {'type': 'table',  # mode, freq, time ...
                        'start': 0,
                        'end': 12,
                        'step': 1},
             'SYMM':True,
            #  'MP':True,
            }
solverdata = fea.Solve(solverset)

print(solverdata['solution']['FREQ'])

postprocset = {"SOLVERDATA": solverdata,
                "COMPUTER": {
                    'structural': {'modes': True},
                    },
                "PLOTSET": {'filename': 'output', 'savepng': True},
                "REPORT": {'log': True, 'get': {'nelem': True, 'nnode': True}},           
                }
postprocdata = fea.PostProcess(postprocset)