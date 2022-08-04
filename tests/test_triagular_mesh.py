'''
TEST WITH TRIANGULAR MESH 
'''
from myfempy.mesh.genmesh import ModelGen
from myfempy.core.solver import gen_static_solution
from myfempy.postprc.postcomp import PostComputer
from myfempy.plots.postplot import postproc_plot
from myfempy.plots.prevplot import preview_plot

# #----------------------------- PRE-PROCESS -----------------------------------#
# mat = {
    # "NAME": "steel",
    # "VXX": 0.33,
    # "EXX": 200E3,
    # "MAT": 'isotropic',
    # "DEF": 'planestress'
# }

# geo = {"NAME": "geo1", "THICKN": 100}

# force = {'DEF': 'forceedge',
         # 'DOF': 'fx',
         # 'DIR': 'edgex',
         # 'LOC': {'x': 4000, 'y': 999, 'z': 0},
         # 'VAL': [100],
         # }

# bondcond = {'DEF': 'fixed',
            # 'DOF': 'all',
            # 'DIR': 'edgex',
            # 'LOC': {'x': 0, 'y': 999, 'z': 0},
            # }

# #-------------------------- GEN WITH LEGACY MESH ----------------------------#

# meshdata = {"LEGACY": {'lx': 4000, 'ly': 1000, 'mesh': 'tria3', 'elem': 'plane31', 'nx': 8, 'ny': 2},
            # "PROPMAT": [mat],
            # "PROPGEO": [geo],
            # "FORCES": [force],
            # "BOUNDCOND": [bondcond],
            # # "QUADRATURE": {'meth': 'gaussian', 'npp': 3}
            # }

# #-------------------------------- BY GMSH -----------------------------------#

# # points = [[0, 0, 0],
# #           [4000, 0, 0],
# #           [4000, 1000, 0],
# #           [0, 1000, 0]]

# # lines = [[1, 2],
# #          [2, 3],
# #          [3, 4],
# #          [4, 1]]

# # plane = [[1, 2, 3, 4]]

# # meshdata = {"GMSH": {'filename': 'zz_malha_triagular',
# #                      'pointlist': points,
# #                      'linelist': lines,
# #                      'planelist': plane,
# #                      'meshconfig': {
# #                          'mesh': 'quad4',
# #                          'elem': 'plane41',
# #                          'sizeelement': 100,
# #                          'meshmap': {'on': True,
# #                                      'edge': 'all'}}},
# #             "PROPMAT": [mat],
# #             "PROPGEO": [geo],
# #             "FORCES": [force],
# #             "BOUNDCOND": [bondcond],
# #             "QUADRATURE": {'meth': 'gaussian', 'npp': 4}
# #             }

# #-------------------------------- GEN MESH -----------------------------------#

# modelinfo = ModelGen.get_model(meshdata)

# previewset = {'RENDER': {'filename': 'malha_triagular_pre', 'show': False, 'scale': 4, 'savepng': True, 'lines': True, 'plottags': {'point': False}},
              # 'QUALITY': {'show': False, 'method': 1, 'scale': 0.1},
              # 'LABELS': {'show': False, 'lines': True, 'scale': 1},
              # }

# preview_plot(previewset, modelinfo)

# #-------------------------------- SOLVER -------------------------------------#
# solverset = {"SOLVER": 'SLD',
             # 'TOL': 2E-3,
             # "STEPSET": {'type': 'table',  # mode, freq, time ...
                         # 'start': 0,
                         # 'end': 1,
                         # 'step': 1},
             # "TRACKER": {'show': False, 'result2plot': 'displ', 'max': []}
             # }


# solution = gen_static_solution(solverset, modelinfo)

# #----------------------------- POST-PROCESS ----------------------------------#

# postprocset = {
    # "SOLUTION": solution['U'],
    # "COMPUTER": {'displ': True, 'stress': True, 'average': True},
    # # 'step':2
    # "PLOTSET": {'show': False, 'result2plot': {'displ': []}, 'filename': 'malha_triagular_post', 'savepng': True},
    # "TRACKER": {'show': False, 'result2plot': 'stress', 'point': {'x': 6, 'y': 1.5, 'z': 0}}
# }

# postporc_result = PostComputer(modelinfo).main(postprocset)

#----------------------------- VIEW SOLUTION ---------------------------------#
# postproc_plot(postprocset, postporc_result, modelinfo)
