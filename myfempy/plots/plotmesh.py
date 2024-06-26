#!/usr/bin/env python
__doc__ = """
Physics Vtk Plot.
"""
import vedo as vd

# from myfempy.utils.utils import get_logo, get_version


def post_show_mesh(file2plot: str, plotset: dict):
    """vedo code"""

    win = vd.Plotter(
        title="POST-PROCESS", sharecam=False, screensize=(1280, 720), interactive=True
    )
    if plotset["edge"] == True:
        mesh = vd.UnstructuredGrid(file2plot + ".vtk").lineWidth(1).flat().c("blue7")
    else:
        # mesh = vd.Mesh(file2plot + ".vtk").c("blue7")
        mesh = vd.UnstructuredGrid(file2plot + ".vtk")
    cname = vd.colorMap(range(21), "jet")
    mesh.cmap(cname, on=plotset["apply"]).addScalarBar(
        title=plotset["text_plot"], c="w"
    )
    text = vd.Text2D(
        "MYFEMPY "
        # + get_version()
        + " < " + plotset["text_plot"] + ' >\nPress "q" to continue...',
        s=1,
        font="Arial",
        c="white",
    )
    win.show(text, mesh, viewup="y", bg="black", axes=4)
    # win.close()
