from __future__ import annotations

import numpy as np


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
        self.tabmat = tabmat
        self.mat_lib = mat_lib
        return tabmat

    def getTabMat(self, modeldata):
        return SetModel.setTabMat(self, modeldata)

    def setTabGeo(self, modeldata):
        tabgeo, geo_lib = SetModel.__tabgeo(self, modeldata["GEOMETRY"])
        self.tabgeo = tabgeo
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
        self.inci = inci
        self.mesh_type_list = mesh_type_list
        return inci

    def getInci(self, modeldata):
        return SetModel.setInci(self, modeldata)

    def setCoord(self, modeldata):
        coordlist = SetModel.__coordlist(self, modeldata)
        coord = SetModel.__coord(self, coordlist)
        self.coord = coord
        return coord

    def getCoord(self, modeldata):
        return SetModel.setCoord(self, modeldata)

    def getIntGauss(self, modeldata):
        return SetModel.__intgauss(self, modeldata)

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
            "VXX": "VXX",
            "GXX": "GXX",
            "EYY": "EYY",
            "VYY": "VYY",
            "GYY": "GYY",
            "RHO": "RHO",
            "KXX": "KXX",
            "KYY": "KYY",
            "CTE": "CTE",
            "VIS": "VIS",
            "STIF": "STIF",
            "DAMP": "DAMP",
        }
        tabmat = [{}] * nmat  # np.zeros((nmat, len(key_mat_list)))
        for mm in range(nmat):
            mat_lib[matlist["PROPMAT"][mm]["NAME"]] = mm + 1
            for pp in range(len(key_mat_list)):
                key = list(key_mat_list)[pp]
                if key in matlist["PROPMAT"][mm].keys():
                    mat_prop[key] = matlist["PROPMAT"][mm][key]
                else:
                    mat_prop[key] = 0.0
            # matset = self.material.getMaterialSet()
            # idtyp = matset["idtyp"]  # mat_def(matlist[mm]["MAT"])
            # idmat = matset["idmat"]  # mat_beh(matlist[mm]["DEF"])
            # tabmat[mm] = {key:mat_prop[key]}
            tabmat[mm] = {
                "EXX": mat_prop["EXX"],
                "VXX": mat_prop["VXX"],
                "GXX": mat_prop["GXX"],
                "EYY": mat_prop["EYY"],
                "VYY": mat_prop["VYY"],
                "GYY": mat_prop["GYY"],
                "RHO": mat_prop["RHO"],
                "KXX": mat_prop["KXX"],
                "KYY": mat_prop["KYY"],
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
                    "ID": idgeo,
                }

            else:
                for pp in range(len(key_geo_list)):
                    key = list(key_geo_list)[pp]
                    if key in geolist["PROPGEO"][gg].keys():
                        geo_prop[key] = geolist["PROPGEO"][gg][key]
                    else:
                        geo_prop[key] = 0.0

                geoset = self.geometry.GeometrySet()
                idgeo = geoset["idgeo"]

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

    def __intgauss(self, modeldata):
        self.intgauss = modeldata["ELEMENT"]["INTGAUSS"]
        return self.intgauss
