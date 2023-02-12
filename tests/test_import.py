'''
TEST TO IMPORT MODULES FROM MYFEMPY
'''

import myfempy
from myfempy.mesh.genmesh import ModelGen
from myfempy.core.solver import Solver
from myfempy.postprc.postcomp import PostProcess
from myfempy.plots.postplot import postproc_plot
from myfempy.plots.prevplot import preview_plot


print("myfempy version: ",myfempy.__version__)


print("test successful")