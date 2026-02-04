#!/usr/bin/env python
import numpy as np
import vedo as vd

# from myfempy.utils.utils import get_version

class MeshProp:
    def __init__(self, plotset: dict, path):
        self.plotset = plotset
        self.filename = path + "/" + plotset["RENDER"]["filename"]

    def mesh_numbering(self):
        win = vd.Plotter(title="PRE-PROCESS", sharecam=False, screensize=(1280, 720))
        mesh = vd.UnstructuredGrid(self.filename + ".vtk")  # .lineWidth(0.1).flat()
        if self.plotset["LABELS"]["lines"] == False:
            mesh = vd.UnstructuredGrid(self.filename + ".vtk")
        else:
            pass
        # mesh.cmap("RdYlBu", on="cells")

        text = vd.Text2D(
            "MYFEMPY " + " < mesh numb. > ",
            s=1,
            font="Arial",
            c="white",
        )

        nodes = self.plotset["inci"][:, 4 : 4 + self.plotset["nodecon"]].reshape(
            (self.plotset["nnode"] * self.plotset["nodecon"],)
        )

        # noduni, idx = np.unique(nodes, return_index=True)

        labs0 = mesh.labels2d(
            content="id",  # nodes[np.sort(idx)].astype(int),  # 'id'
            on="points",
            scale=self.plotset["LABELS"]["scale"],
            # font="Arial",
            c="black",
        )

        labs1 = mesh.labels2d(
            content="cellid",  # self.plotset['inci'][:, 0].astype(int), # 'cellid'
            on="cells",
            scale=self.plotset["LABELS"]["scale"],
            # font="Arial",
            c="red",
        )

        win.show(text, mesh, labs0, labs1, viewup="y", bg="white", axes=2)

    def mesh_quality(self):
        """_summary_"""
        win = vd.Plotter(
            title="PRE-PROCESS",
            sharecam=False,
            screensize=(1280, 720),
            interactive=True,
        )
        mesh = (
            vd.UnstructuredGrid(self.plotset["RENDER"]["filename"] + ".vtk")
            .lineWidth(0.1)
            .flat()
        )
        # if self.plotset["QUALITY"]["lines"] == False:
        #     mesh = vd.Mesh(self.plotset["RENDER"]["filename"] + ".vtk")
        # else:
        #     pass
        mesh.cmap("RdYlBu", on="cells", n=16).addScalarBar()
        text = vd.Text2D(
            "MYFEMPY < mesh quality > ",
            s=1,
            font="Arial",
            c="white",
        )
        # https://vedo.embl.es/docs/vedo/mesh.html#Mesh.compute_quality
        mesh.compute_quality(metric=self.plotset["QUALITY"]["method"]).cmap(
            "RdYlBu", on="cells"
        ).print()
        mesh.addScalarBar3D(
            c="white", title=("Meth. " + str(self.plotset["QUALITY"]["method"]))
        )
        win.show(text, mesh, viewup="y", bg="black", axes=4)
