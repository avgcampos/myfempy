#!/usr/bin/env python
from myfempy.utils.utils import get_version
import numpy as np
import vedo as vd

__doc__ = """
Mesh Quality Calc.
"""


class MeshProp:
    """_summary_"""

    def __init__(self, plotset):
        self.plotset = plotset

    def mesh_numbering(self):
        """_summary_"""
        win = vd.Plotter(title="PRE-PROCESS", sharecam=False, screensize=(1280, 720))
        mesh = (
            vd.Mesh(self.plotset["RENDER"]["filename"] + ".vtk").lineWidth(0.1).flat()
        )
        if self.plotset["LABELS"]["lines"] == False:
            mesh = vd.Mesh(self.plotset["RENDER"]["filename"] + ".vtk")
        else:
            pass
        mesh.cmap("RdYlBu", on="cells", n=4)
        text = vd.Text2D(
            "MYFEMPY " + get_version() + " < mesh numb. > ",
            s=1,
            font="Arial",
            c="white",
        )
        nodes = self.plotset["inci"][:, 4 : 4 + self.plotset["nodecon"]].reshape(
            (self.plotset["nnode"] * self.plotset["nodecon"],)
        )
        noduni, idx = np.unique(nodes, return_index=True)
        labs0 = mesh.labels(
            content=nodes[np.sort(idx)].astype(int),  # 'id'
            cells=False,
            scale=self.plotset["LABELS"]["scale"],
            font="Arial",
            c="white",
        )
        labs1 = mesh.labels(
            "cellid",  # content=self.plotset['inci'][:,0].astype(int), # 'cellid'
            cells=True,
            scale=self.plotset["LABELS"]["scale"],
            font="Arial",
            c="green",
        )
        win.show(text, mesh, labs0, viewup="y", bg="black", axes=4)

    def mesh_quality(self):
        """_summary_"""
        win = vd.Plotter(title="PRE-PROCESS", sharecam=False, screensize=(1280, 720))
        mesh = (
            vd.Mesh(self.plotset["RENDER"]["filename"] + ".vtk").lineWidth(0.1).flat()
        )
        if self.plotset["QUALITY"]["lines"] == False:
            mesh = vd.Mesh(self.plotset["RENDER"]["filename"] + ".vtk")
        else:
            pass
        mesh.cmap("RdYlBu", on="cells", n=16).addScalarBar()
        text = vd.Text2D(
            "MYFEMPY " + get_version() + " < mesh quality > ",
            s=1,
            font="Arial",
            c="white",
        )
        mesh.addQuality(measure=self.plotset["QUALITY"]["method"]).cmap(
            "RdYlBu", on="cells"
        ).print()
        mesh.addScalarBar3D(
            c="white", title=("Meth. " + str(self.plotset["QUALITY"]["method"]))
        )
        win.show(text, mesh, viewup="y", bg="black", axes=4)
