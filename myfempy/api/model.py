from __future__ import annotations

from myfempy.io.controllers import setPoints2NumericalIntegration
import numpy as np

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

class SetModel:
    """Model Class <ClassOrder>"""

    def __init__(self, Mesh, Element, Shape, Material, Geometry):
        self.mesh = Mesh
        self.shape = Shape
        self.material = Material
        self.geometry = Geometry
        self.element = Element

    # -----------------------------------------------
    def setElemList(self, modeldata):
        elemlist = SetModel.__elemlist(self, modeldata)
        self.elemlist = elemlist
        return elemlist

    def getElemList(self, modeldata):
        return SetModel.setElemList(self, modeldata)

    def setTabMat(self, modeldata):
        tabmat, mat_lib = SetModel.__tabmat(self, modeldata["MATERIAL"])
        # self.tabmat = tabmat
        self.mat_lib = mat_lib
        return tabmat

    def getTabMat(self, modeldata):
        return SetModel.setTabMat(self, modeldata)

    def setTabGeo(self, modeldata):
        tabgeo, geo_lib = SetModel.__tabgeo(self, modeldata["GEOMETRY"])
        # self.tabgeo = tabgeo
        self.geo_lib = geo_lib
        return tabgeo

    def getTabGeo(self, modeldata):
        return SetModel.setTabGeo(self, modeldata)

    def setInci(self, modeldata):
        SetModel.setElemList(self, modeldata)
        SetModel.setTabMat(self, modeldata)
        SetModel.setTabGeo(self, modeldata)

        inci, mesh_type_list = SetModel.__inci(
            self, self.elemlist, self.mat_lib, self.geo_lib
        )
        # self.inci = inci
        self.mesh_type_list = mesh_type_list
        return inci

    def getInci(self, modeldata):
        return SetModel.setInci(self, modeldata)

    def setCoord(self, modeldata):
        coordlist = SetModel.__coordlist(self, modeldata)
        coord = SetModel.__coord(self, coordlist)
        # self.coord = coord
        return coord

    def getCoord(self, modeldata):
        return SetModel.setCoord(self, modeldata)
    
    # def setIntGauss(self, element):
    #     intgauss = SetModel.__intgauss(element)
    #     self.intgauss = intgauss
    #     return intgauss

    # def getIntGauss(self, element):
    #     return SetModel.setIntGauss(self, element)

    # -----------------------------------------------
    # privates methods
    def __elemlist(self, modeldata):
        set_mesh = modeldata["MESH"]
        conec = self.mesh.getElementConection(set_mesh)
        elemset = self.element.getElementSet()
        shapset = self.shape.getShapeSet()
        meshset = int(f'{elemset["id"]}{shapset["id"]}')
        elemlist = self.mesh.getElementList(conec, meshset, modeldata)
        return elemlist

    def __coordlist(self, modeldata):
        # self.mesh(modeldata)
        set_mesh = modeldata["MESH"]
        coord = self.mesh.getNodesCoord(set_mesh)
        coordlist = [[None] * 4]
        nodes = [[None] * 4]
        for nn in range(len(coord)):
            nodes.append([int(coord[nn][0]), coord[nn][1], coord[nn][2], coord[nn][3]])
        nodes = nodes[1::][::]
        coordlist.extend(nodes)
        coordlist = coordlist[1::][::]
        return coordlist

    def __tabmat(self, matlist):
        """get material table"""
        nmat = len(matlist["PROPMAT"])
        mat_lib = dict()
        mat_prop = dict()
        key_mat_list = {
            "EXX": "EXX",
            "VXY": "VXY",
            "GXY": "GXY",
            "EYY": "EYY",
            "VYZ": "VYZ",
            "GYZ": "GYZ",
            "EZZ": "EZZ",
            "VZX": "VZX",
            "GZX": "GZX",
            "RHO": "RHO",
            "KXX": "KXX",
            "KYY": "KYY",
            "KZZ": "KZZ",
            "CTE": "CTE",
            "VIS": "VIS",
            "STIF": "STIF",
            "DAMP": "DAMP",
        }
        tabmat = [{}] * nmat 
        for mm in range(nmat):
            mat_lib[matlist["PROPMAT"][mm]["NAME"]] = mm + 1
            for pp in range(len(key_mat_list)):
                key = list(key_mat_list)[pp]
                if key in matlist["PROPMAT"][mm].keys():
                    mat_prop[key] = matlist["PROPMAT"][mm][key]
                else:
                    mat_prop[key] = "NULL"
            tabmat[mm] = {
                "EXX": mat_prop["EXX"],
                "VXY": mat_prop["VXY"],
                "GXY": mat_prop["GXY"],
                "EYY": mat_prop["EYY"],
                "VYZ": mat_prop["VYZ"],
                "GYZ": mat_prop["GYZ"],
                "EZZ": mat_prop["EZZ"],
                "VZX": mat_prop["VZX"],
                "GZX": mat_prop["GZX"],
                "RHO": mat_prop["RHO"],
                "KXX": mat_prop["KXX"],
                "KYY": mat_prop["KYY"],
                "KZZ": mat_prop["KZZ"],
                "CTE": mat_prop["CTE"],
                "VIS": mat_prop["VIS"],
                "STIF": mat_prop["STIF"],
                "DAMP": mat_prop["DAMP"],
            }
        return tabmat, mat_lib

    def __tabgeo(self, geolist):
        """get geometry table"""
        ngeo = len(geolist["PROPGEO"])
        geo_lib = dict()
        geo_prop = dict()
        key_geo_list = {
            "AREACS": "AREACS",
            "INERYY": "INERYY",
            "INERZZ": "INERZZ",
            "INERXX": "INERXX",
            "THICKN": "THICKN",
            "B": "b",
            "H": "h",
            "T": "t",
            "D": "d",
            "YMAX": "y_max",
            "YMIN": "y_min",
            "ZMAX": "z_max",
            "ZMIN": "z_min",
            "RMAX": "r_max",
            "ID": "ID",
        }
        tabgeo = [{}] * ngeo  # np.zeros((ngeo, len(key_geo_list) + 1))
        for gg in range(ngeo):
            geo_lib[geolist["PROPGEO"][gg]["NAME"]] = gg + 1

            if "DIM" in geolist["PROPGEO"][gg].keys():
                b = geolist["PROPGEO"][gg]["DIM"]["b"]
                h = geolist["PROPGEO"][gg]["DIM"]["h"]
                t = geolist["PROPGEO"][gg]["DIM"]["t"]
                d = geolist["PROPGEO"][gg]["DIM"]["d"]

                dim_sec = {
                    "b": b,
                    "h": h,
                    "t": t,
                    "d": d,
                }

                geoset = self.geometry.GeometrySet()
                idgeo = geoset["idgeo"]

                sect_prop = self.geometry.getSectionProp(dim_sec)

                tabgeo[gg] = {
                    "AREACS": sect_prop["areacs"],
                    "INERZZ": sect_prop["inerzz"],
                    "INERYY": sect_prop["ineryy"],
                    "INERXX": sect_prop["inerxx"],
                    "THICKN": sect_prop["thickn"],
                    "B": b,
                    "H": h,
                    "T": t,
                    "D": d,
                    "YMAX": 1.0,
                    "YMIN": -1.0,
                    "ZMAX": 1.0,
                    "ZMIN": -1.0,
                    "RMAX": 1.0,
                    "ID": idgeo,
                }

            else:
                for pp in range(len(key_geo_list)):
                    key = list(key_geo_list)[pp]
                    if key in geolist["PROPGEO"][gg].keys():
                        geo_prop[key] = geolist["PROPGEO"][gg][key]
                    else:
                        geo_prop[key] = "NULL"

                geoset = self.geometry.GeometrySet()
                idgeo = geoset["idgeo"]

                y_max, y_min, z_max, z_min, r_max = "NULL", "NULL", "NULL", "NULL", "NULL"

                if 'CG' in geolist["PROPGEO"][gg].keys():
                    y_max = geolist["PROPGEO"][gg]["CG"]["y_max"]
                    y_min = geolist["PROPGEO"][gg]["CG"]["y_min"]
                    z_max = geolist["PROPGEO"][gg]["CG"]["z_max"]
                    z_min = geolist["PROPGEO"][gg]["CG"]["z_min"]
                    r_max = geolist["PROPGEO"][gg]["CG"]["r_max"]
                else:
                    pass

                tabgeo[gg] = {
                    "AREACS": geo_prop["AREACS"],
                    "INERZZ": geo_prop["INERZZ"],
                    "INERYY": geo_prop["INERYY"],
                    "INERXX": geo_prop["INERXX"],
                    "THICKN": geo_prop["THICKN"],
                    "B": geo_prop["B"],
                    "H": geo_prop["H"],
                    "T": geo_prop["T"],
                    "D": geo_prop["D"],
                    "YMAX": y_max,
                    "YMIN": y_min,
                    "ZMAX": z_max,
                    "ZMIN": z_min,
                    "RMAX": r_max,
                    "ID": idgeo,
                }
        return tabgeo, geo_lib

    def __inci(self, elemlist, mat_lib, geo_lib):
        MAXCONECELM = int(20)  # Max 20-Nodes from hexa20
        inci = [[None] * (1 + 3 + MAXCONECELM)]
        nelem = len(elemlist)
        conec_elm = np.zeros((nelem, MAXCONECELM + 1), dtype=np.int32)
        prop_elm = np.zeros((nelem, 3), dtype=np.int32)
        mesh_type_list = dict()
        contelm = int(0)
        for kk in range(nelem):
            contelm += 1
            conec_elm[kk, 0] = contelm  # elemlist[dm][kk][0]
            nodes = np.array(elemlist[kk][4])
            if len(nodes) != MAXCONECELM:
                nodes = np.append(
                    nodes, np.zeros(int(MAXCONECELM - len(nodes))), axis=0
                )
            conec_elm[kk, 1:] = nodes

            elemeset = self.element.getElementSet()
            shapeset = self.shape.getShapeSet()
            keyelem = elemlist[kk][1]
            # elem = get_elemset(keyelem)
            # elemset = elem.elemset()
            mesh_type_list[keyelem] = [
                int(f'{elemeset["id"]}{shapeset["id"]}'),  # elemeset["id"],
                len(elemeset["dofs"]),
                len(shapeset["nodes"]),
                len(elemeset["tensor"]),
            ]
            prop_elm[kk, 0] = int(elemlist[kk][1])
            prop_elm[kk, 1] = mat_lib[elemlist[kk][2]]
            prop_elm[kk, 2] = geo_lib[elemlist[kk][3]]
        inci = np.concatenate(
            (conec_elm[:, 0][:, np.newaxis], prop_elm, conec_elm[:, 1:]), axis=1
        )
        return inci, mesh_type_list

    def __coord(self, coordlist):
        nnod = len(coordlist)
        coord = np.zeros((nnod, 4))
        for ii in range(0, nnod):
            coord[ii, :] = np.array(coordlist[ii][:])
        return coord

    # def __intgauss(element):
    #     if "INTGAUSS" in element.keys():
    #         intgauss = element["INTGAUSS"]
    #     else:
    #         intgauss = setPoints2NumericalIntegration(element["SHAPE"])
    #     return intgauss