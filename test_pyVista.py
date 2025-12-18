# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 23:26:45 2025

@author: antvi
"""

import pyvista as pv
# from pyvirtualdisplay import Display

# display = Display(visible=0, size=(1280, 1024))
# display.start()

p = pv.Plotter()
# m = pv.read("out/plane_stress.vtk")
m = pv.read("out/beso_myfempy_topevo_results_step-53.vtk")
p.add_mesh(m, show_edges=False)
pts = m.points
p.show(window_size=[512, 384], cpos="xy")

# display.stop()