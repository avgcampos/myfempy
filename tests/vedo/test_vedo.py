# -*- coding: utf-8 -*-
"""MYFEMPY          
SHOW
"""

import numpy as np
import vedo as vd

# "depth peeling" may improve the rendering of transparent objects
# vd.settings.useDepthPeeling = True
# vd.settings.multiSamples = 0  # needed on OSX vtk9

mesh = vd.Mesh('solid_stress_vm.vtk').lineWidth(0.1).flat()

# generate a numpy array for mesh quality
# mesh.addQuality(measure=17)

# hist = vd.pyplot.histogram(mesh.celldata["Quality"], xtitle='mesh quality', bc='w')
# # make it smaller and position it, useBounds makes the cam
# # ignore the object when resetting the 3d qscene
# hist.scale(1.6).pos(0,0,0).useBounds(False)

# labs = mesh.labels('id',
#                    cells=True,
#                    scale=1.0,
#                    c='black')

# create numeric labels of active scalar on top of cells
labs = mesh.labels(cells=True,
                    precision=3,
                    scale=0.010,
                    c='black')


# pids = mesh.boundaries(returnPointIds=True)
# bpts = mesh.points()[pids]

# pts = vd.Points(bpts, r=20, c='red')

pts1 = [200, 0, 100]
pts2 = [200, 20, 100]

arrows = vd.shapes.Arrow(startPoint=pts2, endPoint=pts1, s=None, c='black', alpha=1, res=8)


mycmap="jet"
mesh.cmap(mycmap, on='cells', n=16).addScalarBar()

vd.show(mesh, labs, arrows, __doc__, viewup='y', bg='white', axes=3, screensize=(1280, 720)).close()
