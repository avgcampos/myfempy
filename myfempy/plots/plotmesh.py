#!/usr/bin/env python
__doc__ = """
Physics Vtk Plot.
"""
import vedo as vd
from myfempy.tools.tools import get_version, get_logo


def post_show_mesh(file2plot, plotset):
    win = vd.Plotter(title="POST-PROCESS", sharecam=False, screensize=(1280, 720))
    mesh = vd.Mesh(file2plot + ".vtk").lineWidth(1).flat()
    if plotset["edge"] == False:
        mesh = vd.Mesh(file2plot + ".vtk")
    else:
        pass
    cname = vd.colorMap(range(21), "jet")
    mesh.cmap(cname, on=plotset["apply"]).addScalarBar(
        title=plotset["text_plot"], c="w"
    )
    text = vd.Text2D(
        "MYFEMPY "
        + get_version()
        + " < "
        + plotset["text_plot"]
        + ' >\nPress "q" to continue...',
        s=1,
        font="Arial",
        c="white",
    )
    win.show(text, mesh, viewup="y", bg="black", axes=4)
    # win.close()
