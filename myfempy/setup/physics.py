from __future__ import annotations

import numpy as np

from myfempy.core.utilities import (get_elemen_from_nodelist,
                                    get_nodes_from_list)


class SetPhysics:
    """Set Physics Class <ClassOrder>"""

    def __init__(self, Model, Loads, BoundCond):
        self.loads = Loads
        self.boundcond = BoundCond
        self.model = Model
        # self.domain = Domain

    def getForceList(self, physicdata):
        return SetPhysics.setForceList(self, physicdata)

    def setForceList(self, physicdata):
        forcelist = SetPhysics.__forcelist(physicdata["LOAD"])
        self.forcelist = forcelist
        return forcelist

    def getBoundCondList(self, physicdata):
        return SetPhysics.setBoundCondList(self, physicdata)

    def setBoundCondList(self, physicdata):
        boundlist = SetPhysics.__boundcondlist(physicdata["BOUNDCOND"])
        self.boundcondlist = boundlist
        return boundlist

    def getLoadApply(self, physicdata):
        flist = SetPhysics.setForceList(self, physicdata["PHYSIC"])
        forcenodeaply = np.zeros((1, 4))
        for nforc in range(len(flist)):
            forcelist = flist[nforc]
            fapp = self.loads.getLoadApply(self.model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        forcenodeaply = forcenodeaply[1::][::]
        forcenodeaply[forcenodeaply[:, 3].argsort()]
        return forcenodeaply

    def getBoundCondApply(self, physicdata):
        bclist = SetPhysics.setBoundCondList(self, physicdata["PHYSIC"])
        boncdnodeaply = np.zeros((1, 4))
        for nbc in range(len(bclist)):
            listapply = bclist[nbc]
            bcapp = self.boundcond.getBCApply(self.model, listapply)
            boncdnodeaply = np.append(boncdnodeaply, bcapp, axis=0)
        boncdnodeaply = boncdnodeaply[1::][::]
        boncdnodeaply[boncdnodeaply[:, 3].argsort()]
        return boncdnodeaply

    def getLoadCoup(self, physicdata):
        forcenodeaply = np.zeros((1, 4))
        for nforc in range(len(physicdata['COUPLING']["POST"])):
            coupling = physicdata['COUPLING']["POST"][nforc]
            coupling["TYPE"] = physicdata['COUPLING']["TYPE"]
            coupling["STEP"] = int(nforc + 1)
            fapp = self.loads.getLoadApply(self.model, coupling)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        forcenodeaply = forcenodeaply[1::][::]
        forcenodeaply[forcenodeaply[:, 3].argsort()]
        return forcenodeaply

    def getUpdateMatrix(self, matrix, loadaply):
        return self.loads.getUpdateMatrix(self.model, matrix, loadaply)

    def getUpdateLoad(self):
        return None

    def getNodeList(self, domain_nodelist, coord, regions):
        return get_nodes_from_list(domain_nodelist, coord, regions)

    def getElementList(self, inci, nodelist):
        return get_elemen_from_nodelist(inci, nodelist)

    # -----------------------------------------------
    # privates methods
    def __forcelist(forcelist):
        """force set"""
        nforc = len(forcelist)
        flist = []  # np.zeros((1, 9))
        for fl in range(nforc):
            fap = forcelist[fl]
            for fs in range(len(fap["VAL"])):
                if "LOC" in fap.keys():

                    flist.append(
                        {
                            "TYPE": fap["TYPE"],
                            "DOF": fap["DOF"],
                            "VAL": fap["VAL"][fs],
                            "DIR": fap["DIR"],
                            "LOCX": fap["LOC"]["x"],
                            "LOCY": fap["LOC"]["y"],
                            "LOCZ": fap["LOC"]["z"],
                            "TAG": 0,
                            "STEP": int(fs + 1),
                        }
                    )

                elif "TAG" in fap.keys():

                    flist.append(
                        {
                            "TYPE": fap["TYPE"],
                            "DOF": fap["DOF"],
                            "VAL": fap["VAL"][fs],
                            "DIR": fap["DIR"],
                            "LOCX": 0.0,
                            "LOCY": 0.0,
                            "LOCZ": 0.0,
                            "TAG": fap["TAG"],
                            "STEP": int(fs + 1),
                        }
                    )

                else:
                    flist.append(
                        {
                            "TYPE": fap["TYPE"],
                            "DOF": 'fx',
                            "VAL": fap["VAL"][fs],
                            "DIR": 'node',
                            "LOCX": 0.0,
                            "LOCY": 0.0,
                            "LOCZ": 0.0,
                            "TAG": 0,
                            "STEP": int(fs + 1),
                        }
                    )
        return flist

    def __boundcondlist(boundcondlist):
        """boundary conditions set"""
        nbound = len(boundcondlist)
        blist = []  # np.zeros((1, 9))
        for bl in range(nbound):
            bap = boundcondlist[bl]
            if "VAL" in bap.keys():
                for bs in range(len(bap["VAL"])):
                    if "LOC" in bap.keys():
                        blist.append(
                            {
                                "TYPE": bap["TYPE"],
                                "DOF": bap["DOF"],
                                "DIR": bap["DIR"],
                                "LOCX": bap["LOC"]["x"],
                                "LOCY": bap["LOC"]["y"],
                                "LOCZ": bap["LOC"]["z"],
                                "TAG": 0,
                                "VAL": bap["VAL"][bs],
                                "STEP": int(bs + 1),
                            }
                        )

                    elif "TAG" in bap.keys():
                        blist.append(
                            {
                                "TYPE": bap["TYPE"],
                                "DOF": bap["DOF"],
                                "DIR": bap["DIR"],
                                "LOCX": 0.0,
                                "LOCY": 0.0,
                                "LOCZ": 0.0,
                                "TAG": bap["TAG"],
                                "VAL": bap["VAL"][bs],
                                "STEP": int(bs + 1),
                            }
                        )
                    else:
                        blist.append(
                            {
                                "TYPE": bap["TYPE"],
                                "DOF": bap["DOF"],
                                "DIR": bap["DIR"],
                                "LOCX": 0.0,
                                "LOCY": 0.0,
                                "LOCZ": 0.0,
                                "TAG": 0,
                                "VAL": bap["VAL"][bs],
                                "STEP": int(bs + 1),
                            }
                        )
            else:
                if "LOC" in bap.keys():
                    if "STEP" in bap.keys():
                        blist.append(
                            {
                                "TYPE": bap["TYPE"],
                                "DOF": bap["DOF"],
                                "DIR": bap["DIR"],
                                "LOCX": bap["LOC"]["x"],
                                "LOCY": bap["LOC"]["y"],
                                "LOCZ": bap["LOC"]["z"],
                                "TAG": 0,
                                "VAL": 0.0,
                                "STEP": bap["STEP"],
                            }
                        )
                    else:
                        blist.append(
                            {
                                "TYPE": bap["TYPE"],
                                "DOF": bap["DOF"],
                                "DIR": bap["DIR"],
                                "LOCX": bap["LOC"]["x"],
                                "LOCY": bap["LOC"]["y"],
                                "LOCZ": bap["LOC"]["z"],
                                "TAG": 0,
                                "VAL": 0.0,
                                "STEP": 0,
                            }
                        )
                    
                elif "TAG" in bap.keys():
                    if "STEP" in bap.keys():
                        blist.append(
                            {
                                "TYPE": bap["TYPE"],
                                "DOF": bap["DOF"],
                                "DIR": bap["DIR"],
                                "LOCX": 0.0,
                                "LOCY": 0.0,
                                "LOCZ": 0.0,
                                "TAG": bap["TAG"],
                                "VAL": 0.0,
                                "STEP": bap["STEP"],
                            }
                        )
                    else:
                        blist.append(
                            {
                                "TYPE": bap["TYPE"],
                                "DOF": bap["DOF"],
                                "DIR": bap["DIR"],
                                "LOCX": 0.0,
                                "LOCY": 0.0,
                                "LOCZ": 0.0,
                                "TAG": bap["TAG"],
                                "VAL": 0.0,
                                "STEP": 0,
                            }
                        )
                else:
                    blist.append(
                        {
                            "TYPE": bap["TYPE"],
                            "DOF": 'ux',
                            "DIR": 'node',
                            "LOCX": 0.0,
                            "LOCY": 0.0,
                            "LOCZ": 0.0,
                            "TAG": 0,
                            "VAL": 0.0,
                            "STEP": 0,
                        }
                    )
        return blist
