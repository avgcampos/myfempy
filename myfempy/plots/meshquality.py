#!/usr/bin/env python
import numpy as np
import vedo as vd

# from myfempy.utils.utils import get_version


__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""

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
