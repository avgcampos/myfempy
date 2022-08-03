#!/usr/bin/env python
"""
myfempy -- MultiphYsics Finite Element Method with PYthon
Copyright (C) 2022 Antonio Vinicius Garcia Campos
"""
from myfempy.tools.tools import *
from myfempy.tools.path import *
from myfempy.postprc.postset import *
from myfempy.postprc.postcomp import *
from myfempy.postprc.displcalc import *
from myfempy.plots.prevplot import *
from myfempy.plots.postplot import *
from myfempy.plots.plotxy import *
from myfempy.plots.plotmesh import *
from myfempy.plots.physics import *
from myfempy.plots.meshquality import *
from myfempy.mesh.legacy import *
from myfempy.mesh.gmsh import *
from myfempy.mesh.genmesh import *
from myfempy.io.iovtk import *
from myfempy.io.iomsh import *
from myfempy.io.filters import *
from myfempy.felib.struct.truss21 import *
from myfempy.felib.struct.spring21 import *
from myfempy.felib.struct.solid81 import *
from myfempy.felib.struct.solid41 import *
from myfempy.felib.struct.plane41 import *
from myfempy.felib.struct.plane31 import *
from myfempy.felib.struct.frame22 import *
from myfempy.felib.struct.frame21 import *
from myfempy.felib.struct.beam21 import *
from myfempy.felib.physics.loadsconstr import *
from myfempy.felib.physics.getnode import *
from myfempy.felib.physics.force2node import *
from myfempy.felib.materials.solid import *
from myfempy.felib.materials.planestress import *
from myfempy.felib.materials.planestrain import *
from myfempy.felib.materials.lumped import *
from myfempy.felib.materials.axial import *
from myfempy.felib.fsi import *
from myfempy.felib.fluid import *
from myfempy.felib.quadrature import *
from myfempy.felib.physicset import *
from myfempy.felib.materset import *
from myfempy.felib.felemset import *
from myfempy.felib.crossec import *
from myfempy.core.vibralinear import *
from myfempy.core.staticlinear import *
from myfempy.core.solverset import *
from myfempy.core.solver import *
from myfempy.core.assembler import *
#-------------- VERSION
from myfempy import version
# --------------
__author__ = "Antonio Vinicius Garcia Campos"
__copyright__ = "Copyright @ 2022, Antonio Vinicius Garcia Campos"
__credits__ = ["Antonio Vinicius Garcia Campos", "3D EasyCAE"]
__license__ = "GPL"
__status__ = "Development"
__website__ = "https://github.com/easycae-3d/myfempy"
__version__ = version.__version__
